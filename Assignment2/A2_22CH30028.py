import subprocess
import os
import requests
from pathlib import Path


### Step 1: DOWNLOAD PDB FILE
def download_pdb(pdb_id, save_path):
    url = f'https://files.rcsb.org/download/{pdb_id}.pdb'
    save_path_PDB = save_path + f"{pdb_id}.pdb"

    response = requests.get(url)
    if response.status_code == 200:
        save_path = save_path_PDB
        with open(save_path, 'w') as file:
            file.write(response.text)
        print(f'PDB file {pdb_id} downloaded successfully')
        print(f'PDB file saved to {save_path_PDB}')
    else:
        print(f'Failed to download PDB file {pdb_id}. Please check the PDB ID.')
    file.close()


### STEP 2: GET NO. OF CHAINS IN THE COMPLEX
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


### STEP 3: EXTRACT CHAINS INTO SEPARATE PDB FILES
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


### STEP 4: RUN NACCESS OF THE CHAINS & COMPLEX
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


### STEP 5: CALCULATE SOLVENT ACCESSIBLE AREA
def calculate__solvent__accessible__area(file_path):
    area = 0.0
    with open (file_path, "r") as file:
        lines = file.readlines()
        for line in lines:
            entry = line.split()
            area = area + float(entry[-2])
    return area
    

### STEP 6: FIND INTERFACE ATOMS (ATOMS WITH ASA REDUCTION >= 0.1 A^2)
def find_interface_atoms(chain_file_path, complex_file_path):
    interface_atoms = []

    chain_asa_dict = {}

    with open(chain_file_path, "r") as chain_file:
        for line in chain_file:
            entry = line.split()
            if len(entry) > 5:
                # Identification of each atom in the chain based on its NAME & RESIDUE ID
                atom_key = (entry[2], entry[5])
                try:
                    chain_asa_dict[atom_key] = float(entry[-2])
                except ValueError:
                    continue

    with open(complex_file_path, "r") as complex_file:
        for line in complex_file:
            parts = line.split()
            if len(parts) > 5:
                # Seeing if its the same atom in the chain file and complex file
                atom_key = (parts[2], parts[5])

                try:
                    complex_asa = float(parts[-2])
                    if atom_key in chain_asa_dict:
                        asa_reduction = chain_asa_dict[atom_key] - complex_asa
                        if asa_reduction > 0.1:
                            atom_name = atom_key[0]
                            residue_ = parts[3]
                            residue_no = atom_key[1]

                            interface_atoms.append((atom_name, residue_, residue_no, asa_reduction))  
                            # (Atom Name, Residue Name, Residue ID, ASA Reduction)
                except ValueError:
                    continue

    return interface_atoms


# STEP 7: GET SOLVATION PARAMETERS
def get_solvation_parameter(atom_name, residue_name):
    
    solvation_params = {
        "C": 16,  # Carbon
        "S": 21,  # Sulfur
        "Charged_N": -50,  # Charged Nitrogen (Arg, Lys)
        "Charged_O": -24,  # Charged Oxygen (Asp, Glu)
        "Other_NO": -6  # Rest Nitrogen & Oxygen
    }

    if atom_name.startswith("C"):  
        return solvation_params["C"]
    elif atom_name.startswith("S"):  
        return solvation_params["S"]
    elif atom_name in ["NE", "NH1", "NH2", "NZ"] and residue_name in ["ARG", "LYS"]:
        return solvation_params["Charged_N"]
    elif atom_name in ["OD1", "OD2", "OE1", "OE2"] and residue_name in ["ASP", "GLU"]:
        return solvation_params["Charged_O"]
    elif atom_name.startswith("N") or atom_name.startswith("O"):  
        return solvation_params["Other_NO"]
    
    return None


# STEP 8: COMPUTE SOLVATION ENERGY
def compute_solvation_energy(interface_atoms):
    
    total_solvation_energy = 0.0

    for atom in interface_atoms:
        atom_name, residue_name, residue_id, asa_reduction = atom
        solvation_param = get_solvation_parameter(atom_name, residue_name)
        solvation_energy = solvation_param * asa_reduction
        
        total_solvation_energy = total_solvation_energy + solvation_energy
    
    return total_solvation_energy


# STEP 9: CALCULATE LRMSD
def find_LRMSD_coordinates(file_path):
    coordinates = {}
        
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("ATOM"):
                entry = line.split()
                if entry[4] == 'B':  # Chain 'B' refers to LIGAND
                    atom_no = entry[1]  # Atom number [Maintain unique atom recoginition]
                    x = float(entry[6])
                    y = float(entry[7])
                    z = float(entry[8])
                    coords = [x, y, z]

                    coordinates[atom_no] = coords

    # Returns coordinates of all the atoms in CHAIN "B" of the target/decoy in a dictionary format
    # Atom_no: [x, y, z] 
    return coordinates  


# STEP 10: CALCULATE iRMSD (Extracts Interface Residue Coordinates)
def find_iRMSD_coordinates(file_path, list_of_interface_atoms):
    coordinates = {}

    # Remove ASA reduction from tuples
    interface_atoms_filtered = [(x[0], x[1], x[2]) for x in list_of_interface_atoms]

    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("ATOM"):
                entry = line.split()
                
                # Ensure valid line format
                if len(entry) < 9:
                    continue  

                atom_no = entry[1]  # Atom number
                atom_name = entry[2]  # Atom name
                residue_name = entry[3]  # Residue name
                residue_no = entry[5]  # Residue sequence number
                x = float(entry[6])
                y = float(entry[7])
                z = float(entry[8])
                coords = [x, y, z]

                # Check if the atom is in the interface list
                if (atom_name, residue_name, residue_no) in interface_atoms_filtered:
                    coordinates[atom_no] = coords

    return coordinates


# STEP 11: CALCULATE FNAT SCORE
def calculate_FNAT_score(interface_residue_no, decoy_interface_residue_no):
    no_of_target_residues = len(interface_residue_no)
    no_of_decoy_residues = len(decoy_interface_residue_no)
    fnat = no_of_decoy_residues/no_of_target_residues
    return fnat

# ---------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------- #


def main(naccess_path, root_folder, path_to_PDB, target_location, target_root):
    Interface__Area__MASTER = []
    Solvation__Energy__MASTER = []
    LRMSD__MASTER = []
    IRMSD__MASTER = []
    Fnat_score_MASTER = []

    for file_name in os.listdir(path_to_PDB):
        if file_name.startswith("complex.") and file_name.endswith(".pdb"):
            new_name = file_name.replace(".", "_", 1)  # Replace the first '.' with '_'
            old_path = os.path.join(path_to_PDB, file_name)
            new_path = os.path.join(path_to_PDB, new_name)

            os.rename(old_path, new_path)  # Rename the file
            print(f"Renamed: {file_name} → {new_name}")

    print("All PDB files renamed successfully!")

    list_of_decoys = []

    for decoy in os.listdir(path_to_PDB):
        print(decoy)
        list_of_decoys.append(decoy)

    score = {}
    for decoy in list_of_decoys:
        score[decoy] = {
            "Interface area": None,
            "Solvation Energy": None,
            "LRMSD": None,
            "IRMSD": None,
            "Fnat Score": None
        }


    print(f"\n \n \n \n \nCALCULATING INTERFACE AREA WITH NACCESS . . . . . Please wait while the program runs . . . \n \n \n \n")

    for decoy in list_of_decoys:
        decoy_name = os.path.splitext(decoy)[0]

        decoy_file_path = f"{path_to_PDB}\\{decoy}"
        
        # Since all the given complexes are dimers with chains labelled A and B, we are not manaully extracting the number and # names of chains from the PDB file
        list_ = ['A', 'B']

        for i in list_:
            chain_id = i
            output_file = f"{path_to_PDB}\\{decoy_name}_chain_{chain_id}.txt"
            extract__chains(decoy_file_path, chain_id, output_file)

        for i in list_:    
            FILE_NAME = f"{decoy_name}_chain_{i}.pdb"
            path_to_naccess = f"{root_folder}/{FILE_NAME}"
            run_naccess(path_to_naccess, naccess_path)

        path_to_naccess = f"{root_folder}/{decoy}"
        run_naccess(path_to_naccess, naccess_path)
        

        solvent__accessible__area = []
        solvent__accessible__area.append(calculate__solvent__accessible__area(f"{path_to_PDB}\\{decoy_name}.asa"))

        for i in  list_:
            chain_id = i
            path = f"{path_to_PDB}\\{decoy_name}_chain_{chain_id}.asa"
            solvent__accessible__area.append(calculate__solvent__accessible__area(path))

        if len(solvent__accessible__area) < 3:
            print(f"Error: Not enough solvent accessibility data for {decoy_name} !!!!!")
            continue 

        complex = solvent__accessible__area[0]
        receptor = solvent__accessible__area[1]
        ligand = solvent__accessible__area[2]
        interface_area = (receptor + ligand - complex)/2

        Interface__Area__MASTER.append(interface_area)
        score[decoy]["Interface area"] = interface_area
        print(f"Interface area for {decoy_name}.pdb = {interface_area} A^2")

    for decoy in list_of_decoys:
        try:
            print(score[decoy].get("Interface area"))
        except Exception as e:
            print(f"Error: {e}")

    
    print(f"\n \n \n \n \nCALCULATING SOLVATION ENERGY FOR DECOYS . . . . . Please wait while the program runs . . . \n \n \n \n")

    for decoy in list_of_decoys:
        decoy_name = os.path.splitext(decoy)[0]
        total_solvation_energy = 0.0

        list_of_interface_atoms = []

        for i in list_:
            chain_file_path = f"{path_to_PDB}\\{decoy_name}_chain_{i}.asa"
            complex_file_path = f"{path_to_PDB}\\{decoy_name}.asa"
            chain_interface_atoms = find_interface_atoms(chain_file_path, complex_file_path)
            
            list_of_interface_atoms.extend(chain_interface_atoms)
            
            total_solvation_energy = total_solvation_energy + compute_solvation_energy(chain_interface_atoms)

        total_solvation_energy = (0.5 * total_solvation_energy)/1000
        print(f"Solvation energy for {decoy_name}.pdb = {total_solvation_energy} kcal/mol")
        Solvation__Energy__MASTER.append(total_solvation_energy)
        score[decoy]["Solvation Energy"] = total_solvation_energy

    for decoy in list_of_decoys:
        try:
            print(score[decoy].get("Solvation Energy"))
        except Exception as e:
            print(f"Error: {e}")



    print(f"\n \n \n \n \nCALCULATING LRMSD BETWEEN DECOYS AND TARGET . . . . . Please wait while the program runs . . . \n \n \n \n")


    for decoy in list_of_decoys:

        list_of_interface_atoms = []

        for i in list_:
            chain_file_path = f"{path_to_PDB}\\{decoy_name}_chain_{i}.asa"
            complex_file_path = f"{path_to_PDB}\\{decoy_name}.asa"
            chain_interface_atoms = find_interface_atoms(chain_file_path, complex_file_path)
            
            list_of_interface_atoms.extend(chain_interface_atoms)

        target_path = target_location
        decoy_path = f"{path_to_PDB}\\{decoy}"

        target_coordinates = find_LRMSD_coordinates(target_path)
        decoy_coordinates = find_LRMSD_coordinates(decoy_path)
        
        RMSD = 0.0
        count = 0

        for atom_no, target_coords in target_coordinates.items():
            if atom_no in decoy_coordinates:
                decoy_coords = decoy_coordinates[atom_no]
                distance_squared = sum((t - d) ** 2 for t, d in zip(target_coords, decoy_coords))
                RMSD += distance_squared
                count += 1

        if count > 0:
            RMSD = (RMSD / count) ** 0.5
        else:
            RMSD = None

        print(f"LRMSD for {decoy} = {RMSD} Å \n")
        LRMSD__MASTER.append(RMSD)
        score[decoy]["LRMSD"] = RMSD


    print(f"\n \n \n \n \nCALCULATING IRMSD BETWEEN DECOYS AND TARGET . . . . . Please wait while the program runs . . . \n \n \n \n")


    for decoy in list_of_decoys:
        target_path = target_location
        decoy_path = f"{path_to_PDB}\\{decoy}"

        target_coordinates = find_iRMSD_coordinates(target_path, list_of_interface_atoms)
        decoy_coordinates = find_iRMSD_coordinates(decoy_path, list_of_interface_atoms)
        
        RMSD = 0.0
        count = 0

        for atom_no, target_coords in target_coordinates.items():
            if atom_no in decoy_coordinates:
                decoy_coords = decoy_coordinates[atom_no]
            
                distance_squared = sum((t - d) ** 2 for t, d in zip(target_coords, decoy_coords))
                RMSD += distance_squared
                count += 1

        if count > 0:
            RMSD = (RMSD / count) ** 0.5
        else:
            RMSD = None

        print(f"iRMSD for {decoy} = {RMSD} Å \n")
        IRMSD__MASTER.append(RMSD)
        score[decoy]["IRMSD"] = RMSD


    print(f"\n \n \n \n \nCALCULATING FNAT SCORE FOR THE DECOYS AND TARGET . . . . . Please wait while the program runs . . . \n \n \n")

    # CREATE SEPARATE CHAINS FOR THE TARGET.pdb
    for i in list_:
        chain_id = i
        output_file = f"{target_root}\\target_chain_{chain_id}.txt"
        extract__chains(target_path, chain_id, output_file)

    # Interface atoms in TARGET . . . 
    # "target_interface_atoms" -- array -- contains TUPLE  
    target_interface_atoms = []
    for i in list_:
        chain_file_path = f"{path_to_PDB}\\{decoy_name}_chain_{i}.asa"
        complex_file_path = f"{path_to_PDB}\\{decoy_name}.asa"
        chain_interface_atoms = find_interface_atoms(chain_file_path, complex_file_path)
        
        target_interface_atoms.extend(chain_interface_atoms)

    interface_residue_no = []
    interface_residues = []

    for item in target_interface_atoms:
        residue_no = item[2]
        residue_name = item[1]

        if residue_no not in interface_residue_no:
            interface_residue_no.append(residue_no)
            interface_residues.append(residue_name)
        else:
            continue

    print(f"No. of interface residues in target.pdb = {len(interface_residue_no)} \n")

    # find_interface_atoms ---> returns --> (atom_name, residue_, residue_no, asa_reduction)

    # Interface atoms in DECOYS . . . 
    for decoy in list_of_decoys:
        decoy_interface_atoms = []

        for i in list_:
            chain_file_path = f"{path_to_PDB}\\{decoy_name}_chain_{i}.asa"
            complex_file_path = f"{path_to_PDB}\\{decoy_name}.asa"
            chain_interface_atoms = find_interface_atoms(chain_file_path, complex_file_path)
            # chain_interface_atoms -- array -- (atom_name, residue_, residue_no, asa_reduction)

            decoy_interface_atoms.extend(chain_interface_atoms)

        # decoy_interface_residue_no -- array -- Contains (atom_name, residue_, residue_no, asa_reduction) for interface atoms
        decoy_interface_residue_no = []

        for item in decoy_interface_atoms:
            residue_no = item[2]
            residue_name = item[1]

            if residue_no not in decoy_interface_residue_no:
                decoy_interface_residue_no.append(residue_no)
            else:
                continue
        
        print(f"No. of interface residues in {decoy} = {len(decoy_interface_residue_no)}")

        # Calculate Fnat score . . . 
        no_of_target_residues = len(interface_residue_no)
        no_of_decoy_residues = len(decoy_interface_residue_no)
        fnat_score = no_of_decoy_residues/no_of_target_residues
        score[decoy]["Fnat Score"] = fnat_score


    return score

if __name__ == "__main__":
    # Paths for WSL - UBUNTU Compatible
    naccess_path = "/mnt/d/______6th_SEMESTER_FILES/COMPUTATIONAL-BIOPHYSICS/naccess/naccess"
    root_folder = "/mnt/d/______6th_SEMESTER_FILES/COMPUTATIONAL-BIOPHYSICS/Problem-2/upload/Decoys8/"

    # Paths for File Handing - WINDOWS Compatible
    path_to_PDB = "Decoys"
    target_location = "D:\\______6th_SEMESTER_FILES\\COMPUTATIONAL-BIOPHYSICS\\Problem-2\\upload\\target.pdb"
    target_root = "D:\\______6th_SEMESTER_FILES\\COMPUTATIONAL-BIOPHYSICS\\Problem-2\\upload"

    score = main(naccess_path, root_folder, path_to_PDB, target_location, target_root)
    
    with open("score.txt", "w") as file:
        
        # HEADER Row . . . 
        file.write("Filename".ljust(20) + "| Interface Area (A^2)".ljust(20) + "| Solvation Energy (kcal/mol)".ljust(20) +
                "| LRMSD (A)".ljust(10) + "| IRMSD (A)".ljust(10) + "| Fnat Score".ljust(10) + "\n")
        
        file.write("=" * 100 + "\n")

        # Write data for each decoy
        for decoy, values in score.items():
            file.write(decoy.ljust(20) + "| " + 
                    str(values["Interface area"]).ljust(20) + "| " +
                    str(values["Solvation Energy"]).ljust(20) + "| " +
                    str(values["LRMSD"]).ljust(10) + "| " +
                    str(values["IRMSD"]).ljust(10) + "| " +
                    str(values["Fnat Score"]).ljust(10) + "\n")

    print(f"✅ Score data successfully written to score.txt")

