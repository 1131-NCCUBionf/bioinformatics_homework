import pandas as pd
import numpy as np

input_path = "example/mut.txt"
output_path = "pamx.csv"  # Change output path to .csv

def generate_pam(x, input_path, output_path):
    # Load the mutation matrix
    mutation = pd.read_csv(input_path, index_col=0, header=0, delim_whitespace=True)
    mutation_matrix = np.array(mutation, dtype=float)

    # Normalize the mutation matrix
    mutation_matrix = mutation_matrix / 10000

    # Calculate the PAM matrix (Mx)
    pam_matrix = np.linalg.matrix_power(mutation_matrix, x)

    # Round the PAM matrix and convert to integers
    rounded_pam_matrix = np.round(pam_matrix * 100).astype(int)

    # Print the rounded PAM matrix
    print(rounded_pam_matrix)

    # Create a DataFrame to save to CSV
    pam_df = pd.DataFrame(rounded_pam_matrix, index=mutation.index, columns=mutation.index)

    # Save the PAM matrix to a CSV file
    pam_df.to_csv(output_path)

    print(f"PAM matrix saved to {output_path}.")

# Example usage
generate_pam(250, input_path=input_path, output_path=output_path)
