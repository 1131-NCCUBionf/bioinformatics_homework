from Bio import SeqIO

from Bio import SeqIO

def read_fasta(input_path):
    """Read a FASTA file and return a list of (header, sequence)."""
    sequences = []
    with open(input_path, 'r') as file:
        for record in SeqIO.parse(file, 'fasta'):
            sequences.append((record.id, str(record.seq)))
    return sequences


def read_scoring_matrix(score_path):
    """Parses a PAM matrix file into a dictionary format."""
    matrix = {}
    amino_acids = []
    with open(score_path, 'r') as file:
        for line in file:
            # Skip comment lines and empty lines
            if line.startswith("#") or not line.strip():
                continue
            # First non-comment line contains amino acid labels
            if not amino_acids:
                amino_acids = line.strip().split()
                continue
            # Process the matrix rows
            parts = line.strip().split()
            row_label = parts[0]
            scores = list(map(int, parts[1:]))
            matrix[row_label] = dict(zip(amino_acids, scores))
    return matrix


def initialize_matrices(seq1, seq2, gap_open, gap_extend, local):
    """Initialize score and traceback matrices for alignment."""
    m, n = len(seq1), len(seq2)
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    traceback_matrix = [[None] * (n + 1) for _ in range(m + 1)]

    if not local:  # Global alignment: initialize first row/column with gap penalties
        for i in range(1, m + 1):
            score_matrix[i][0] = gap_open + (i - 1) * gap_extend
            traceback_matrix[i][0] = 'U'
        for j in range(1, n + 1):
            score_matrix[0][j] = gap_open + (j - 1) * gap_extend
            traceback_matrix[0][j] = 'L'
    return score_matrix, traceback_matrix


def align_sequences(seq1, seq2, scoring_matrix, gap, local):
    """Align two sequences using either global or local alignment."""
    len1, len2 = len(seq1), len(seq2)

    # Initialize the scoring (dp) and traceback matrices
    dp = [[0] * (len2 + 1) for _ in range(len1 + 1)]
    traceback = [[''] * (len2 + 1) for _ in range(len1 + 1)]

    # Initialize matrices for global alignment
    if not local:
        dp[0][0] = 0
        for i in range(1, len1 + 1):
            dp[i][0] = dp[i - 1][0] + gap  # Cumulative gap penalty
            traceback[i][0] = 'U'  # Gap in seq2
        for j in range(1, len2 + 1):
            dp[0][j] = dp[0][j - 1] + gap  # Cumulative gap penalty
            traceback[0][j] = 'L'  # Gap in seq1

    # Fill the matrices
    max_score, start_pos = float('-inf'), (0, 0)
    for i in range(1, len1 + 1):
        for j in range(1, len2 + 1):
            match = dp[i - 1][j - 1] + scoring_matrix[seq1[i - 1]][seq2[j - 1]]
            delete = dp[i - 1][j] + gap  # Gap in seq2
            insert = dp[i][j - 1] + gap  # Gap in seq1

            if local:
                dp[i][j] = max(0, match, delete, insert)
            else:
                dp[i][j] = max(match, delete, insert)

            # Record traceback
            if dp[i][j] == match:
                traceback[i][j] = 'D'
            elif dp[i][j] == delete:
                traceback[i][j] = 'U'
            elif dp[i][j] == insert:
                traceback[i][j] = 'L'
            elif local:
                traceback[i][j] = '0'

            # Track the max score for local alignment
            if local and dp[i][j] > max_score:
                max_score = dp[i][j]
                start_pos = (i, j)

    # For global alignment, start traceback from the bottom-right corner
    if not local:
        start_pos = (len1, len2)

    return traceback, start_pos


def traceback_alignment(traceback_matrix, seq1, seq2, start_pos, local=False):
    """Traceback to get the aligned sequences."""
    aligned_seq1, aligned_seq2 = [], []
    i, j = start_pos

    while i > 0 or j > 0:
        if traceback_matrix[i][j] == 'D':  # Diagonal (match/mismatch)
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 'U':  # Up (gap in seq2)
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            i -= 1
        elif traceback_matrix[i][j] == 'L':  # Left (gap in seq1)
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            j -= 1
        else:
            break

    # Reverse the alignments to construct the final sequences
    return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2))


def alignment(input_path, score_path, output_path, aln, gap):
    sequences = read_fasta(input_path)
    scoring_matrix = read_scoring_matrix(score_path)
    seq1, seq2 = sequences[0][1], sequences[1][1]
    local = (aln == 'local')

    # Perform alignment
    traceback_matrix, start_pos = align_sequences(seq1, seq2, scoring_matrix, gap, local)
    aligned_seq1, aligned_seq2 = traceback_alignment(traceback_matrix, seq1, seq2, start_pos, local)

    # Write results to file
    with open(output_path, 'w') as f:
        f.write(f">{sequences[0][0]}\n{aligned_seq1}\n")
        f.write(f">{sequences[1][0]}\n{aligned_seq2}\n")



# alignment("test_global.fasta", "pam250.txt", "result_global.fasta", "global", -10)
# alignment("test_local.fasta", "pam100.txt", "result_local.fasta", "local", -10)