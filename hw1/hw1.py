import pandas as pd
import numpy as np

input_path = "example/mut.txt"
output_path = "pamx.txt"

def generate_pam(x, input_path, output_path):
    mutation = pd.read_csv(input_path, index_col=0, header=0, delim_whitespace=True)
    # print(mutation)
    mutation_matrix = np.array(mutation, dtype=float)
    # print(mutation_matrix)
    mutation_matrix = mutation_matrix / 10000
    # print(mutation_matrix)

    pam_matrix = np.linalg.matrix_power(mutation_matrix, x)
    # print(pam_matrix)

    target_frequencies = pam_matrix.sum(axis=1)
    frequency_sum = target_frequencies.sum()
    target_frequencies = target_frequencies / frequency_sum
    # print(target_frequencies, frequency_sum)

    pam_log = np.zeros_like(pam_matrix)
    for i in range(pam_matrix.shape[0]):
        for j in range(pam_matrix.shape[0]):
            R = pam_matrix[i][j]/target_frequencies[i]
            if R <= 0:
                pam_log[i][j] = -5
            else: pam_log[i][j] = 10 * np.log10(R)

    print(np.round(pam_log).astype(int))

generate_pam(250, input_path=input_path, output_path=output_path)