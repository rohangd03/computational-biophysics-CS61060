README: Protein Structure Analysis Program
Author: ROHAN GHOSH DASTIDAR	22CH30028

Problem 1: Protein Sequence and Structure Analysis

============================================================================================================
### OVERVIEW ###
This Python program analyzes protein structures from PDB files by:
- Extracting chain sequences from PDB & FASTA files.
- Comparing PDB and FASTA sequences for consistency.
- Identifying chain breaks using Needleman-Wunsch alignment.
- Calculating molecular weight and solvent-accessible surface area (ASA) using NACCESS.
- Generating a comprehensive output report in a text file.
============================================================================================================

### DISCLAIMER ###
This code has been tailored to PARTICULARLY run on Windows systems and it excutes NACCESS through WSL (Ubuntu 24.04.2 LTS). Running on any other OS will not lead to desired results.

============================================================================================================
### HOW TO RUN THE PROGRAM ###
1. Install Dependencies
Make sure you have Python 3.11.4 or newer installed along with the required modules:

RUN --> pip install numpy requests

2. Modify File Paths
Before running the script, update the file paths in main():

if __name__ == "__main__":
    root_folder = "D:\\______6th_SEMESTER_FILES\\COMPUTATIONAL-BIOPHYSICS\\PROBLEM-1"
    naccess_path = "/mnt/d/______6th_SEMESTER_FILES/COMPUTATIONAL-BIOPHYSICS/naccess/naccess"
    wsl_root_folder = "/mnt/d/______6th_SEMESTER_FILES/COMPUTATIONAL-BIOPHYSICS/PROBLEM-1"

    main(root_folder, naccess_path, wsl_root_folder)

root_folder --> Folder where you want the PDB, FASTA, and NACCESS files to be stored
wsl_root_folder --> Contains the same location as {root_folder} but the path is modified to be run by WSL Ubuntu
naccess_path --> Contains the location of the naccess file that calculates Solvent Accessible Surface Area

3. Run the Program
Go to the folder containing the PYTHON program --> Open in Terminal

RUN --> python A1_22CH30028.py

4. View Results
The program will generate a report in the {OUTPUT} folder present inside the {root_folder} FOLDER.
============================================================================================================


============================================================================================================
### FUNCTION DESCRIPTIONS ###

FUNCTION NAME							DESCRIPTION

download_pdb(pdb_id, save_path)					Downloads the PDB file from RCSB.
download_fasta(pdb_id, save_path)				Downloads the FASTA sequence from RCSB.
residue_code_to_single_letter(residue_code, atom_count)		Converts 3-letter amino acid codes to single-letter codes.
heteromer__fasta__sequence(fasta_file_path, target_chain)	Extracts sequence of a specific chain from a FASTA file.
homomer__fasta__sequence(fasta_file_path)			Extracts sequence from a homomeric complex.
protein_sequence__pdb1(pdb_file_path, chain_id)			Extracts amino acid sequence from a PDB file.
compare_sequences(pdb_sequence, fasta_sequence)			Compares PDB and FASTA sequences for consistency.
chain_breaks(seq1, seq2)					Implements Needleman-Wunsch alignment to identify chain breaks.
check_chain_breaks(pdb_sequence, fasta_sequence)		Identifies and prints chain breaks in protein sequences.
get__no__of__chains(pdb_file)					Extracts the number and names of chains from a PDB file.
extract__chains(pdb_file, chain_id, output_file)		Extracts individual chains from a PDB file into separate PDB files.
run_naccess(file_path, naccess_path)				Runs NACCESS to calculate solvent-accessible surface area (ASA).
calculate__solvent__accessible__area(file_path)			Extracts ASA values from NACCESS output files.
molecular_weight__calculator(seq)				Calculates the molecular weight of a protein sequence.
============================================================================================================


============================================================================================================
### OUTPUT REPORT ###

The output file would look something like:

Chain Analysis for PDB ID: 8wix
The complex is heteromeric. 
Chain breaks detected in the complex. 
======================================================================================================================================================
Chain ID  Residues       FASTA Mol. Weight (Da)        PDB Mol. Weight (Da)          Solvent Area (Ų)   
------------------------------------------------------------------------------------------------------------------------------------------------------
A         141            19223.30                      18257.00                      7628.09             
B         146            19223.30                      19223.30                      7997.68             
======================================================================================================================================================


============================================================================================================
### NOTES ###
- Ensure NACCESS is installed and executable in WSL.
- If the PDB or FASTA sequence is missing, the program will exit with an error.
- The program will auto-remove old output files before running new analyses.




