import pandas as pd
import numpy as np

def generate_pam(x, input_path, output_path):
    mutation = pd.read_csv(input_path, index_col=0, header=0, sep='\s+', skiprows=1)
    # print(mutation)
    mutation_matrix = np.array(mutation, dtype=float) / 10000
    # print(mutation_matrix)

    pam_matrix = np.linalg.matrix_power(mutation_matrix, x)
    # print(pam_matrix * 100)

    target_frequencies = [0.087, 0.041, 0.040, 0.047, 0.033, 
                          0.038, 0.050, 0.089, 0.034, 0.037,
                          0.085, 0.081, 0.015, 0.040, 0.051,
                          0.070, 0.058, 0.010, 0.030, 0.065]
    # print(sum(target_frequencies))

    pam_log = np.zeros_like(pam_matrix)
    for i in range(pam_matrix.shape[0]):
        for j in range(pam_matrix.shape[0]):
            pam_log[i][j] = 10 * np.log10(pam_matrix[i][j]/target_frequencies[i])

    rounded_pam_log = np.round(pam_log).astype(int)

    # # Create a DataFrame with appropriate headers
    amino_acids = list(mutation.columns)  # Get the amino acid labels
    pam_df = pd.DataFrame(rounded_pam_log, index=amino_acids, columns=amino_acids)

    # Save the output in the desired format
    with open(output_path, 'w') as f:
        # Write the header
        f.write('   ' + ' '.join(amino_acids) + '\n')
        # Write each row of the PAM matrix
        for idx, (aa, row) in enumerate(zip(amino_acids, pam_df.values)):
            if idx == len(amino_acids) - 1:
                f.write(f"{aa} " + ' '.join(map(str, row)))
            else:
                f.write(f"{aa} " + ' '.join(map(str, row)) + '\n')

# generate_pam(250, input_path="example/mut.txt", output_path="example/pamx.txt")