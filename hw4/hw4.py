import numpy as np
from Bio import SeqIO

def parse_score_matrix(score_path):
    """Parse the PAM score matrix from the file."""
    with open(score_path) as f:
        lines = f.readlines()

    matrix = {}
    headers = []
    for line in lines:
        # Skip comments and empty lines
        if line.startswith("#") or not line.strip():
            continue

        parts = line.split()

        # If headers are not yet set, initialize them
        if not headers:
            headers = parts
            continue

        # Process a row of the matrix
        row_char = parts[0]
        scores = list(map(int, parts[1:]))
        for col_char, score in zip(headers, scores):
            matrix[(row_char, col_char)] = score

    return matrix

def needleman_wunsch(seq1, seq2, score_matrix, gap_open, gap_extend):
    """Perform global alignment using the Needleman-Wunsch algorithm with three matrices."""
    m, n = len(seq1), len(seq2)

    # Initialize scoring matrices
    M = np.full((m + 1, n + 1), -np.inf)
    Ix = np.full((m + 1, n + 1), -np.inf)
    Iy = np.full((m + 1, n + 1), -np.inf)

    # Initialize traceback matrices
    traceback_M = np.zeros((m + 1, n + 1), dtype=int)
    traceback_Ix = np.zeros((m + 1, n + 1), dtype=int)
    traceback_Iy = np.zeros((m + 1, n + 1), dtype=int)

    # Base case initialization
    M[0, 0] = 0
    for i in range(1, m + 1):
        Ix[i, 0] = gap_open + (i - 1) * gap_extend
        traceback_Ix[i, 0] = 1  # Indicating a vertical gap
    for j in range(1, n + 1):
        Iy[0, j] = gap_open + (j - 1) * gap_extend
        traceback_Iy[0, j] = 2  # Indicating a horizontal gap

    # Fill scoring matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score_matrix.get((seq1[i - 1], seq2[j - 1]), -np.inf)

            # Update M
            M_scores = [
                M[i - 1, j - 1] + match,
                Ix[i - 1, j - 1] + match,
                Iy[i - 1, j - 1] + match,
            ]
            M[i, j] = max(M_scores)
            traceback_M[i, j] = np.argmax(M_scores)

            # Update Ix
            Ix_scores = [
                M[i - 1, j] + gap_open,
                Iy[i - 1, j] + gap_open,
                Ix[i - 1, j] + gap_extend,
            ]
            Ix[i, j] = max(Ix_scores)
            traceback_Ix[i, j] = np.argmax(Ix_scores)

            # Update Iy
            Iy_scores = [
                M[i, j - 1] + gap_open,
                Ix[i, j - 1] + gap_open,
                Iy[i, j - 1] + gap_extend,
            ]
            Iy[i, j] = max(Iy_scores)
            traceback_Iy[i, j] = np.argmax(Iy_scores)

    # Traceback to get aligned sequences
    aligned_seq1, aligned_seq2 = [], []
    i, j = m, n
    final_score = max(M[m, n], Ix[m, n], Iy[m, n])
    if M[m, n] >= Ix[m, n] and M[m, n] >= Iy[m, n]:
        current_matrix = 0  # Prefer M
    elif Ix[m, n] >= Iy[m, n]:
        current_matrix = 1
    else:
        current_matrix = 2

    while i > 0 or j > 0:
        if current_matrix == 0:  # M matrix (match/mismatch)
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            if traceback_M[i, j] == 0:  # Came from M
                current_matrix = 0
            elif traceback_M[i, j] == 1:  # Came from Ix
                current_matrix = 1
            elif traceback_M[i, j] == 2:  # Came from Iy
                current_matrix = 2
            i -= 1
            j -= 1

        elif current_matrix == 1:  # Ix matrix (gap in seq2)
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            if traceback_Ix[i, j] == 0:  # Came from M
                current_matrix = 0
            elif traceback_Ix[i, j] == 1:  # Came from Ix
                current_matrix = 1
            i -= 1

        elif current_matrix == 2:  # Iy matrix (gap in seq1)
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            if traceback_Iy[i, j] == 0:  # Came from M
                current_matrix = 0
            elif traceback_Iy[i, j] == 2:  # Came from Iy
                current_matrix = 2
            j -= 1

    print("Score b4 retrun:", final_score)
    print("M Score:", M[m, n], "Ix Score:", Ix[m, n], "Iy Score:", Iy[m, n])
    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), final_score

def smith_waterman(seq1, seq2, score_matrix, gap_open, gap_extend):
    """Perform local alignment using the Smith-Waterman algorithm."""
    m, n = len(seq1), len(seq2)
    
    # Initialize scoring and traceback matrices
    score = np.zeros((m + 1, n + 1))
    traceback = np.zeros((m + 1, n + 1), dtype=int)
    
    max_score = 0
    max_positions = []

    # Fill matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1, j - 1] + score_matrix.get((seq1[i - 1], seq2[j - 1]), -np.inf)
            delete = score[i - 1, j] + (gap_open if traceback[i - 1, j] != 1 else gap_extend)
            insert = score[i, j - 1] + (gap_open if traceback[i, j - 1] != 2 else gap_extend)
            score[i, j] = max(0, match, delete, insert)

            if score[i, j] == match:
                traceback[i, j] = 0
            elif score[i, j] == delete:
                traceback[i, j] = 1
            elif score[i, j] == insert:
                traceback[i, j] = 2

            if score[i, j] > max_score:
                max_score = score[i, j]
                max_positions = [(i, j)]
            elif score[i, j] == max_score:
                max_positions.append((i, j))

    # Traceback to get all aligned sequences with the maximum score
    best_alignments = []
    for max_pos in max_positions:
        aligned_seq1, aligned_seq2 = [], []
        i, j = max_pos
        while i > 0 and j > 0 and score[i, j] > 0:
            if traceback[i, j] == 0:  # diagonal
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                i -= 1
                j -= 1
            elif traceback[i, j] == 1:  # up
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append('-')
                i -= 1
            elif traceback[i, j] == 2:  # left
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j - 1])
                j -= 1

        best_alignments.append((''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), max_score))
    
    # Find the longest alignment among those with the max score
    longest_alignment = max(best_alignments, key=lambda x: len(x[0]))

    return longest_alignment[0], longest_alignment[1], longest_alignment[2]

def alignment(input_path, score_path, output_path, aln, gap_open, gap_extend):
    """Perform global or local sequence alignment."""
    # Parse sequences from input file
    records = list(SeqIO.parse(input_path, "fasta"))
    seq1, seq2 = str(records[0].seq), str(records[1].seq)
    print("sep1: ", seq1)
    print("sep2: ", seq2)

    # Parse scoring matrix
    score_matrix = parse_score_matrix(score_path)

    # Perform alignment
    if aln == "global":
        aligned_seq1, aligned_seq2, score = needleman_wunsch(seq1, seq2, score_matrix, gap_open, gap_extend)
    elif aln == "local":
        aligned_seq1, aligned_seq2, score = smith_waterman(seq1, seq2, score_matrix, gap_open, gap_extend)
    else:
        raise ValueError("Invalid alignment type. Use 'global' or 'local'.")

    # Write output
    with open(output_path, "w") as f:
        f.write(f">{records[0].id}\n{aligned_seq1}\n")
        f.write(f">{records[1].id}\n{aligned_seq2}\n")

    print("Aligned sequence 1:", aligned_seq1)
    print("Aligned sequence 2:", aligned_seq2)

    print("Score after return:", score)
    return score

# result_score = alignment("examples/test.fasta", "examples/pam250.txt", "examples/result_global.fasta", "global", -10, -2)
# threshold = 45

# if result_score >= threshold:
#     print("Alignment score meets the threshold!")
# else:
#     print("Alignment score is below the threshold.")

# result_score = alignment("examples/test.fasta", "examples/pam250.txt", "examples/result_local.fasta", "local", -10, -2)
# threshold = 59

# if result_score >= threshold:
#     print("Alignment score meets the threshold!")
# else:
#     print("Alignment score is below the threshold.")