from itertools import combinations

def read_sequences(input_path):
    sequences = {}
    with open(input_path, 'r') as f:
        seq_id = None
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                seq_id = line[1:]  # Sequence ID without '>'
                sequences[seq_id] = ""
            else:
                sequences[seq_id] += line  # Add sequence line by line
    return sequences

def read_pam_matrix(score_path):
    pam_matrix = {}
    with open(score_path, 'r') as f:
        # Skip header lines until we find the matrix header row
        for line in f:
            line = line.strip()
            if line and not line.startswith('#') and line[0].isalpha():
                headers = line.split()  # Get amino acid labels from this line
                break

        # Read the matrix lines
        for line in f:
            if line.strip() and not line.startswith('#'):
                parts = line.split()
                row_header = parts[0]  # The amino acid label for this row
                scores = list(map(int, parts[1:]))  # Convert scores to integers
                pam_matrix[row_header] = dict(zip(headers, scores))
                
    return pam_matrix


def calculate_pairwise_score(seq1, seq2, pam_matrix, gopen, gextend):
    score = 0
    gap_open = False  # Tracks if we're currently in a gap
    
    for a, b in zip(seq1, seq2):
        if a == '-' or b == '-':  # If there is a gap
            if gap_open:
                score += gextend  # Extend the gap penalty
            else:
                score += gopen  # Open a new gap penalty
                gap_open = True
        else:
            score += pam_matrix[a][b]  # Match/mismatch score from PAM
            gap_open = False  # Reset gap status

    return score

def calculate_SoP(input_path, score_path, gopen, gextend):
    sequences = read_sequences(input_path)
    pam_matrix = read_pam_matrix(score_path)
    total_score = 0
    
    # Calculate scores for all unique pairs of sequences
    for (seq_id1, seq1), (seq_id2, seq2) in combinations(sequences.items(), 2):
        total_score += calculate_pairwise_score(seq1, seq2, pam_matrix, gopen, gextend)
    
    return total_score


# total_score1 = calculate_SoP("examples/test1.fasta", "examples/pam250.txt", -10, -2) #score=1047
# total_score2 = calculate_SoP("examples/test2.fasta", "examples/pam100.txt", -8, -2) #score=606

# print("Total score of alignment:", total_score1)
# print("Total score of alignment:", total_score2)