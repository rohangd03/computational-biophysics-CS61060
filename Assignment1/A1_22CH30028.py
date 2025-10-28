import subprocess
import os
import requests
from pathlib import Path
import numpy as np
import sys
import time
import math

### Step 1: DOWNLOAD PDB FILE
def download_pdb(pdb_id, save_path):
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    save_path_PDB = save_path
    
    response = requests.get(url)
    if response.status_code == 200:
        # Request was successful
        # Save the file

        with open(save_path_PDB, 'w') as file:
            file.write(response.text)
        print(f'PDB file {pdb_id} downloaded successfully')
        print(f"PDB file save to {save_path_PDB}")
    else:
        print("Failed to download PDB file. Please check the PDB ID.")
        time.sleep(3)
        sys.exit()
    file.close()

### STEP 2: DOWNLOAD FASTA FILE
def download_fasta(pdb_id, save_path):
    fasta_url = f'https://www.rcsb.org/fasta/entry/{pdb_id}'
    response = requests.get(fasta_url)
    save_path_FASTA = save_path
    
    if(response.status_code == 200):
        with open(save_path_FASTA, 'w') as file:
            file.write(response.text)
        print(f"FASTA file downloaded successfully")
    else:
        print("Failed to download FASTA file. Please check the PDB ID.")
    file.close()

### STEP 3: CONVERT THREE LETTER CODE TO SINGLE LETTER CODE
def residue_code_to_single_letter(residue_code, atom_count):
    three_aa = [
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL"]

    single_aa = [
        "A", "R", "N", "D", "C",
        "Q", "E", "G", "H", "I",
        "L", "K", "M", "F", "P",
        "S", "T", "W", "Y", "V"]
    
    index = three_aa.index(residue_code)
    try:
        return single_aa[index]
    except Exception as e:
        print(f"Error raised for {residue_code} at {atom_count}: {e}")

### STEP 4(a): EXTRACT PROTEIN SEQUENCE FROM FASTA
def heteromer__fasta__sequence(fasta_file_path, target_chain):

    with open(fasta_file_path, "r") as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip()

            if line.startswith(">"):  # Sequence header
                # Save the sequence that has been generated so far . . . 
                if target_chain in line:
                    found_chain = True
                    sequence = "" # Reset sequence
                else:
                    found_chain = False
            
            elif found_chain:
                sequence += line

    return sequence

### STEP 4(b): EXTRACT PROTEIN SEQUENCE FROM FASTA
def homomer__fasta__sequence(fasta_file_path):

    with open(fasta_file_path, "r") as file:
            lines = file.readlines()
            for line in lines:
                line = line.strip()

                if line.startswith(">"):  # Sequence header 
                    sequence = ""
                else:
                    sequence += line

    return sequence

### STEP 5: EXTRACT PROTEIN SEQUENCE FROM PDB FILE
def protein_sequence__pdb1(pdb_file_path, i):
    with open(pdb_file_path, "r") as file:
        
        sequence = ""
        lines = file.readlines()
        for line in lines:
            if line.startswith("ATOM"):
                entry = line.split()
                atom_count = entry[1]

                if entry[4] == i:
                    if entry[2] == 'CA':
                        residue_3_letter = entry[3]
                        try:
                            singleLetter = residue_code_to_single_letter(residue_3_letter, atom_count)
                            sequence += singleLetter
                        except Exception as e:
                            print(f"Error raised for {residue_3_letter} at {atom_count}: {e}")
                            continue
        
        return sequence

### STEP 6: COMPARE PDB AND FASTA SEQUENCE
def compare_sequences(pdb_sequence, fasta_sequence):
    for chain_id in pdb_sequence:
        if chain_id in fasta_sequence:
            if pdb_sequence[chain_id] == fasta_sequence[chain_id]:
                print(f"Sequence for chain {chain_id} is same in PDB and FASTA")
                return 1
            else:
                print(f"Sequence for chain {chain_id} is different in PDB and FASTA")
                return 0
        else:
            return Exception("Chain not found in FASTA file")

### STEP 7(a): FIND THE POSITIONS OF CHAIN BREAKS
### THIS IS IMPLEMENTED USING THE [[ NEEDLEMAN-WUNSCH ALGORITHM ]]
def chain_breaks(seq1, seq2):
    match = 1
    mismatch = -1
    gap = -2

    m, n = len(seq1), len(seq2)
    
    # SCORING MATRIX
    score_matrix = np.zeros((m + 1, n + 1))
    
    for i in range(m + 1):
        score_matrix[i][0] = i * gap
    for j in range(n + 1):
        score_matrix[0][j] = j * gap

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match_mismatch = score_matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)
            delete = score_matrix[i - 1][j] + gap
            insert = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(match_mismatch, delete, insert)

    aligned_seq1, aligned_seq2 = "", ""
    i, j = m, n
    
    while i > 0 or j > 0:
        current_score = score_matrix[i][j]

        if i > 0 and j > 0 and (current_score == score_matrix[i - 1][j - 1] + (match if seq1[i - 1] == seq2[j - 1] else mismatch)):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and (current_score == score_matrix[i - 1][j] + gap):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    return aligned_seq1, aligned_seq2

### STEP 7(b): CHECK IF THERE ARE ANY CHAIN BREAKS
def check_chain_breaks(pdb_sequence, fasta_sequence):
    for chain_id in pdb_sequence:
        if chain_id in fasta_sequence:
            pdb_seq = pdb_sequence[chain_id]
            fasta_seq = fasta_sequence[chain_id]
            if pdb_seq != fasta_seq:

                print(f"Chain breaks detected in chain {chain_id}:")
                aligned_seq1, aligned_seq2 = chain_breaks(pdb_seq, fasta_seq)
                
                print(f"Aligned PDB sequence: {'-' * 50} \n {aligned_seq1}")                
                print(f"{aligned_seq2} \n Aligned FASTA sequence: {'-' * 50}")

                location = []
                for i in range(len(aligned_seq1)):
                    if aligned_seq1[i] == "-":
                        location.append(i)
                print(f"PDB chain breaks at positions: {[i for i in location]}")
                return 1


            else:
                print(f"No chain breaks in chain {chain_id}.")
                return 0
        else:
            print(f"Chain {chain_id} not found in FASTA file.")

### STEP 8(a): GET NO. OF CHAINS IN THE COMPLEX
def get__no__of__chains(pdb_file):
    list_of_chain_id = []
    
    with open(pdb_file, "r") as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("ATOM"):
                entry = line.split()
                if entry[4] not in list_of_chain_id:
                    list_of_chain_id.append(entry[4])

    return list_of_chain_id

### STEP 8(b): EXTRACT CHAINS TO SEPARATE PDB FILES
def extract__chains(pdb_file, chain_id, output_file):
    with open(output_file, "w") as f:
        with open(pdb_file, "r") as file:
            for line in file:
                if line.startswith("ATOM"):
                    parts = line.split()
                    if len(parts) > 4 and parts[4] == chain_id:
                        f.write(line)

    new_output_file = output_file.replace(".txt", ".pdb")
    
    if os.path.exists(new_output_file):
        os.remove(new_output_file)

    os.rename(output_file, new_output_file)
    
    print(f"Chain {chain_id} extracted to {new_output_file}")

### STEP 9: CALCULATE SURFACE AREA WITH NACCESS
def run_naccess(file_path, naccess_path):
    
    # file_path - WSL path
    # naccess_path - WSL path
    # file_name - windows path
    # root_path - windows path

    try:
        subprocess.run(["wsl", "tcsh", naccess_path, file_path])
        print(f"NACCESS executed successfully for {file_path}!")
    except Exception as e:
        print(f"Unexpected error: {e}")

### STEP 10: CALCULATE SOLVENT ACCESSIBLE AREA
def calculate__solvent__accessible__area(file_path):
    area = 0.0
    with open (file_path, "r") as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("ATOM"):
                entry = line.split()
                area = area + float(entry[-2])
    return area

### STEP 11: CALCULATE MOLECULAR WEIGHT OF EACH CHAIN
def molecular_weight__calculator(seq):
    amino_acid_weights = {
    'A': 89.1,   # Alanine
    'R': 174.2,  # Arginine
    'N': 132.1,  # Asparagine
    'D': 133.1,  # Aspartate
    'C': 121.2,  # Cysteine
    'E': 147.1,  # Glutamate
    'Q': 146.2,  # Glutamine
    'G': 75.1,   # Glycine
    'H': 155.2,  # Histidine
    'I': 131.2,  # Isoleucine
    'L': 131.2,  # Leucine
    'K': 146.2,  # Lysine
    'M': 149.2,  # Methionine
    'F': 165.2,  # Phenylalanine
    'P': 115.1,  # Proline
    'S': 105.1,  # Serine
    'T': 119.1,  # Threonine
    'W': 204.2,  # Tryptophan
    'Y': 181.2,  # Tyrosine
    'V': 117.1   # Valine
    }

    molecular_weight = 0.0 
    for res in seq:
        molecular_weight = molecular_weight + amino_acid_weights[res]
    
    return molecular_weight


### ---------------------------------------------------------------------- ###
### ---------------------------------------------------------------------- ###

def main(root_folder, naccess_path, wsl_root_folder):
    pdb_id = input("Enter PDB ID: ")
    os.makedirs(root_folder, exist_ok=True)

    pdb_file_path = f"{root_folder}\\{pdb_id}.pdb"
    fasta_file_path = f"{root_folder}\\{pdb_id}.fasta"
    
    download_pdb(pdb_id, pdb_file_path)
    download_fasta(pdb_id, fasta_file_path)

    list_ = get__no__of__chains(pdb_file_path)
    # This gives the NUMBER and NAMES of chains in the complex
    # We check if the complex is homomeric
    with open(fasta_file_path, "r") as file:
        content = file.read()
        frequency = content.count(">")
    

    pdb_sequence = {}
    fasta_sequence = {}

    # Add chain sequences into respective dictionaries
    for i in list_:
        pdb_sequence[i] = protein_sequence__pdb1(pdb_file_path, i)

    # COMPLEX == HOMOMERIC ??? OR HETEROMERIC ???
    if frequency == len(list_):
        typeComplex = "The complex is heteromeric."
        for i in list_:
            fasta_sequence[i] = heteromer__fasta__sequence(fasta_file_path, i)
    else:
        typeComplex = "The complex is homomeric."
        for i in list_:
            fasta_sequence[i] = homomer__fasta__sequence(fasta_file_path)
    print(typeComplex)


    ### COMPARE SEQUENCES - PDB vs FASTA
    output = compare_sequences(pdb_sequence, fasta_sequence)
    if (output == 0):
        check_chain_breaks(pdb_sequence, fasta_sequence)

    for i in list_:
        chain_id = i
        output_file = f"{root_folder}\\{pdb_id}_chain_{i}.txt"
        pdb_file = f"{root_folder}\\{pdb_id}.pdb"
        extract__chains(pdb_file, chain_id, output_file)
    
    
    solvent_accessible_area = {}
    FASTA_molecular_weight = {}
    PDB_molecular_weight = {}
    no_of_residues = {}
    mol_wt_difference_info = {}
    

    for i in list_:    
        FILE_NAME = f"{pdb_id}_chain_{i}"
        path_to_naccess = f"{wsl_root_folder}/{FILE_NAME}.pdb"
        
        if os.path.exists(f"{root_folder}\\{FILE_NAME}.pdb"): 
            print(f"Calculating solvent accessible area for {FILE_NAME}.pdb \n")
            run_naccess(path_to_naccess, naccess_path)
        else:
            print(f"{FILE_NAME}.pdb NOT FOUND . . . \n")
            continue

        
        naccess_file = f"{root_folder}\\{FILE_NAME}.asa"

        if os.path.exists(naccess_file):
            print(f"Naccess output file found --> {FILE_NAME}.asa \n")
            area = calculate__solvent__accessible__area(naccess_file)
            solvent_accessible_area[i] = area
        else:
            print(f"{FILE_NAME}.asa NOT FOUND . . . \n")


        # CALCULATING MOLECULAR WEIGHT USING FASTA SEQUENCE
        FASTA_molecular_weight[i] = molecular_weight__calculator(fasta_sequence[i])
        no_of_residues[i] = len(fasta_sequence[i])

        # CALCULATING MOLECULAR WEIGHT USING FASTA SEQUENCE
        PDB_molecular_weight[i] = molecular_weight__calculator(pdb_sequence[i])
        no_of_residues[i] = len(pdb_sequence[i])

        # COMPARING MOLECULAR WEIGHTS
        if FASTA_molecular_weight[i] == PDB_molecular_weight[i]:
            print(f"Chain {i} - Molecular weights are same")
            mol_wt_difference_info[i] = 0
        else:
            difference = ((FASTA_molecular_weight[i] - PDB_molecular_weight[i])/FASTA_molecular_weight[i])*100
            difference = math.fabs(difference)
            print(f"Chain {i} - PDB Molecular weight differs from FASTA by {difference:.2f} %")
            mol_wt_difference_info[i] = difference


    # WRITING ALL THE OUTPUTS INTO A TEXT FILEEEEEEEEEEEEEE 
    output_folder = f"{root_folder}\\OUTPUT"
    output_file = os.path.join(output_folder, f"{pdb_id}_OUTPUT.txt")
    
    if os.path.exists(output_file):
        os.remove(output_file)

    with open(output_file, "w") as file:
        file.write(f"Chain Analysis for PDB ID: {pdb_id}\n")
        
        if(len(list_) > 1):
            file.write(f"{typeComplex} \n")
            if(output == 0):
                file.write("Chain breaks detected in the complex. \n")
                
                for i in list_:
                    if (mol_wt_difference_info[i] != 0):
                        file.write(f"Chain {i} - PDB Molecular weight differs from FASTA by {mol_wt_difference_info[i]:.2f} % \n")
            
            elif(output == 1):
                file.write("No chain breaks detected in the complex. \n")
        else:
            file.write("This protein contains 1 chain. \n")

        file.write("=" * 150 + "\n")
        
        file.write(f"{'Chain ID':<10}{'Residues':<15}{'FASTA Mol. Weight (Da)':<30}{'PDB Mol. Weight (Da)':<30}{'Solvent Area (A^2)':<20}\n")
        
        file.write("-" * 150 + "\n")
        
        for chain_id in list_:
            file.write(f"{chain_id:<10}{no_of_residues[chain_id]:<15}{FASTA_molecular_weight[chain_id]:<30.2f}{PDB_molecular_weight[chain_id]:<30.2f}{solvent_accessible_area[chain_id]:<20.2f}\n")
        
        file.write("=" * 150 + "\n")
    
    print(f"Chain analysis report saved to {output_file} \n")


if __name__ == "__main__":
    root_folder = "D:\\______6th_SEMESTER_FILES\\COMPUTATIONAL-BIOPHYSICS\\PROBLEM-1"
    naccess_path = "/mnt/d/______6th_SEMESTER_FILES/COMPUTATIONAL-BIOPHYSICS/naccess/naccess"
    wsl_root_folder = "/mnt/d/______6th_SEMESTER_FILES/COMPUTATIONAL-BIOPHYSICS/PROBLEM-1"

    main(root_folder, naccess_path, wsl_root_folder)