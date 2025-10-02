import os
from Bio import SeqIO

def multifasta_to_single_fasta(input_file, output_folder):
    """
    Convert a multifasta file to multiple single FASTA files.
    
    Parameters:
    input_file (str): Path to the input multifasta file.
    output_folder (str): Path to the folder where single fasta files will be saved.
    """
    # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Open the multifasta file
    try:
        with open(input_file, "r") as multifasta:
            # Parse the multifasta file and iterate through each sequence
            for record in SeqIO.parse(multifasta, "fasta"):
                # Generate the output file name based on the sequence ID
                output_filename = os.path.join(output_folder, f"{record.id}.fasta")
                
                # Write the individual sequence to a separate FASTA file
                with open(output_filename, "w") as single_fasta:
                    SeqIO.write(record, single_fasta, "fasta")
                    
                print(f"Written: {output_filename}")
    
    except Exception as e:
        print(f"Error processing the file: {e}")

# Example usage:
input_file = "/Users/kirstyn.brunker/GitHub/RABV_Philippines_2025/redcap/rabv_4b_2025_r10_run7/rabv_4b_2025_r10_run7.fasta"  # Path to the input multifasta file
output_folder = "/Users/kirstyn.brunker/GitHub/RABV_Philippines_2025/redcap/rabv_4b_2025_r10_run7/split_fasta"  # Path to the output folder where single fasta files will be saved

# Call the function
multifasta_to_single_fasta(input_file, output_folder)
