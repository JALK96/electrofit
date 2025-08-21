import glob
import logging
import os
import shutil
from openbabel import openbabel
import json
import sys
import re
import numpy as np

def find_project_root(current_dir, project_name="electrofit"):
    root = None
    while True:
        parent_dir = os.path.dirname(current_dir)
        if os.path.basename(current_dir) == project_name:
            root = current_dir  # Set root to the current project_name directory
        if parent_dir == current_dir:
            # We've reached the filesystem root
            if root is None:
                raise FileNotFoundError(f"Project root directory '{project_name}' not found.")
            return root  # Return the outermost match found
        current_dir = parent_dir

script_dir = os.path.dirname(os.path.abspath(__file__))
project_path = find_project_root(current_dir=script_dir)

def replace_posres_in_file(file_path):
    """Opens a topology file and replaces"POSRES_LIG" with "POSRES"

    Args:
        file_path (file): Topology file to eddit.
    """
    try:
        # Open the file and read the content
        with open(file_path, "r") as file:
            content = file.read()

        # Check if POSRES_LIG is present, then replace it
        if "POSRES_LIG" in content:
            updated_content = content.replace("POSRES_LIG", "POSRES")
        else:
            logging.error("No 'POSRES_LIG' found in the file.")
            return

        # Write the updated content back to the file
        with open(file_path, "w") as file:
            file.write(updated_content)

        logging.info(f"Successfully replaced 'POSRES_LIG' with 'POSRES' in {file_path}")

    except FileNotFoundError:
        logging.error(f"The file {file_path} does not exist.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")


def include_tip3p(file_path):
    try:
        # Open the file and read the content
        with open(file_path, "r") as file:
            content = file.read()

        # Check if '[ system ]' is present
        if "[ system ]" in content:
            # Find where to insert the tip3p topology (itp file)
            updated_content = content.replace(
                "[ system ]",
                """
; Include water topology
#include "amber14sb.ff/tip3p.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

[ system ]
                                              """,
            )
        else:
            logging.error(
                "No '[ system ]' section found in the file to include tip3p topology (tip3p.itp)."
            )
            return

        # Write the updated content back to the file
        with open(file_path, "w") as file:
            file.write(updated_content)

        logging.info(
            f"Successfully included tip3p.itp and position restraints in {file_path}"
        )

    except FileNotFoundError:
        logging.error(f"The file {file_path} does not exist.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")


def include_ions(file_path):
    try:
        # Open the file and read the content
        with open(file_path, "r") as file:
            content = file.read()

        # Check if '[ system ]' is present
        if "[ system ]" in content:
            # Find where to insert the topology (itp file)
            updated_content = content.replace(
                "[ system ]",
                """
; Include ion topology
#include "amber14sb.ff/ions.itp"

[ system ]
                                              """,
            )
        else:
            logging.error(
                "No '[ system ]' section found in the file to include ions.itp)."
            )
            return

        # Write the updated content back to the file
        with open(file_path, "w") as file:
            file.write(updated_content)

        logging.info(f"Successfully included ions.itp in {file_path}")

    except FileNotFoundError:
        logging.error(f"The file {file_path} does not exist.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")


def remove_defaults_section_lines(file_path):
    """
    Removes the '[ defaults ]' section and the next three lines from a GROMACS topology file.

    Parameters:
    - file_path (str): Path to the GROMACS topology (.top) file.

    Returns:
    - None

    Raises:
    - FileNotFoundError: If the specified file does not exist.
    - Exception: For any other unforeseen errors during file processing.
    """
    try:
        # Check if the file exists
        if not os.path.isfile(file_path):
            raise FileNotFoundError(f"The file '{file_path}' does not exist.")

        # Read the original content of the file
        with open(file_path, "r") as file:
            lines = file.readlines()

        # Initialize variables to track the removal process
        new_lines = []
        skip_next = 0  # Counter to skip lines after [ defaults ]

        for index, line in enumerate(lines):
            stripped_line = line.strip()

            # Detect the start of the '[ defaults ]' section
            if stripped_line == "[ defaults ]" and skip_next == 0:
                logging.info(
                    "Found '[ defaults ]' section. Initiating removal of this section and the next three lines."
                )
                skip_next = 3  # [ defaults ] + 2 subsequent lines
                continue  # Skip the '[ defaults ]' line

            if skip_next > 0:
                logging.debug(f"Skipping line {index + 1}: {line.strip()}")
                skip_next -= 1
                continue  # Skip the lines within the [ defaults ] section

            # Retain all other lines
            new_lines.append(line)

        # Write the updated content back to the file
        with open(file_path, "w") as file:
            file.writelines(new_lines)

        logging.info(
            f"Successfully removed '[ defaults ]' section and the following three lines from '{file_path}'."
        )

    except FileNotFoundError as fnf_error:
        logging.error(fnf_error)
        raise
    except Exception as e:
        logging.error(f"An error occurred while removing '[ defaults ]' section: {e}")
        raise


def include_ff(file_path):
    try:
        # Open the file and read all lines
        with open(file_path, "r") as file:
            lines = file.readlines()

        # Define the lines to insert
        include_comment = "; Include forcefield\n"
        include_line = '#include "amber14sb.ff/forcefield.itp"\n'

        # Check if the include line already exists to prevent duplication
        if include_line in lines:
            logging.info(
                f'"amber14sb.ff/forcefield.itp" is already included in {file_path}.'
            )
            return

        # Insert the include lines as the second and third lines
        if len(lines) >= 1:
            lines.insert(1, include_comment)
            lines.insert(2, include_line)
        else:
            # If the file is empty, add the include lines at the beginning
            lines.append(include_comment)
            lines.append(include_line)

        # Write the updated content back to the file
        with open(file_path, "w") as file:
            file.writelines(lines)

        logging.info(
            f'Successfully included "amber14sb.ff/forcefield.itp" in {file_path}.'
        )

    except FileNotFoundError:
        logging.error(f"The file {file_path} does not exist.")
    except Exception as e:
        logging.error(f"An error occurred: {e}")


# Dictionary containing the binary to IP mapping
binary_to_ip = {
    "010101": "IP21",
    "101010": "IP42",
    "101100": "IP44",
    "111000": "IP56",
    "010111": "IP23",
    "101101": "IP45",
    "111001": "IP57",
    "011111": "IP31",
    "111011": "IP59",
    "101111": "IP47",
    "111101": "IP61"
}

def copy_and_rename_folders(source, destination, bash_script_source=f"{project_path}/electrofit/bash/pis.sh", python_script_source=f"{project_path}/electrofit/execution/check_symmetry_and_write.py", nested_folder="run_gau_create_gmx_in"):
    """
    Copy folders from the source directory to the destination with an additional nested folder.
    Rename .mol2 files within each folder according to a binary-to-IP mapping.
    
    Parameters:
    - source (str): Path to the source directory.
    - destination (str): Path to the destination directory.
    - nested_folder (str): Name of the intermediate folder to nest files (default is "run_gau_create_gmx_in").
    """

    # Ensure the source directory exists
    if not os.path.isdir(source):
        print(f"Source directory '{source}' does not exist.")
        return

    # Create the destination directory if it does not exist
    os.makedirs(destination, exist_ok=True)

    # Walk through each folder in the source directory
    for folder_name in os.listdir(source):
        folder_path = os.path.join(source, folder_name)

        # Only process directories
        if os.path.isdir(folder_path):

            # Define the new nested path in the destination
            nested_dest_path = os.path.join(destination, folder_name, nested_folder)

            # Create the nested directory structure
            os.makedirs(nested_dest_path, exist_ok=True)

            # Copy all files from the source folder to the nested destination folder
            for item in os.listdir(folder_path):
                item_source_path = os.path.join(folder_path, item)
                item_dest_path = os.path.join(nested_dest_path, item)

                # Copy file or directory
                if os.path.isfile(item_source_path):
                    shutil.copy2(item_source_path, item_dest_path)
                elif os.path.isdir(item_source_path):
                    shutil.copytree(item_source_path, item_dest_path, dirs_exist_ok=True)

            print(f"Copied '{folder_path}' to '{nested_dest_path}'")

            # Copy bash script (to process initial structure) to destination directory
            bash_script_name= os.path.basename(bash_script_source)
            bash_dest_path = os.path.join(nested_dest_path, bash_script_name)
            try:
                shutil.copy2(bash_script_source, bash_dest_path)
                print(f"Copied {bash_script_source} to {bash_dest_path}")
            except Exception as e:
                print(f"Error copying file: {e}")

def rename_mol2_binary(base_dir, binary):
    """
    Rename .mol2 files in the specified folder according to a binary-to-IP mapping.
    
    Parameters:
    - base_dir (str): Base directory containing the folder with binary name to be renamed.
    - binary (str): The binary string corresponding to this folder.
    """

    # Get the IP name corresponding to the binary
    ip_name = binary_to_ip.get(binary)
    if ip_name is None:
        print(f"Binary '{binary}' not found in mapping.")
        return

    # Construct the full path to the folder
    folder_path = os.path.join(base_dir, "run_gau_create_gmx_in")
    
    # Check if the folder exists
    if os.path.isdir(folder_path):
        # Find the mol2 file in the folder (assuming there's only one mol2 file)
        for file_name in os.listdir(folder_path):
            if file_name.endswith(".mol2"):
                # Construct the old file path
                old_file_path = os.path.join(folder_path, file_name)
                
                # Construct the new file name based on the IP name
                new_file_name = f"{ip_name}.mol2"
                
                # Construct the new file path
                new_file_path = os.path.join(folder_path, new_file_name)
                
                # Rename the mol2 file
                os.rename(old_file_path, new_file_path)
                print(f"Renamed '{old_file_path}' to '{new_file_path}'")
    else:
        print(f"Folder '{folder_path}' does not exist or does not match the binary name in mapping.")

def find_file_with_extension(extension):

    cwd = os.getcwd()
    search_pattern = os.path.join(cwd, f"*.{extension}")
    files = glob.glob(search_pattern)

    if files:
        file_name = os.path.basename(files[0])  # Take the first found result
        print(f"Found file: {file_name}")
        return file_name
    else:
        print(f"No files with .{extension} extension found in {cwd}")
        return None

def get_parent_folder_name():
    current_dir = os.getcwd()
    parent_dir = os.path.dirname(current_dir)
    parent_folder_name = os.path.basename(parent_dir)
    return parent_folder_name

def get_parent_parent_folder_name():
    current_dir = os.getcwd()
    parent_dir = os.path.dirname(current_dir)
    parent_parent_dir = os.path.dirname(parent_dir)
    parent_parent_folder_name = os.path.basename(parent_parent_dir)
    return parent_parent_folder_name


def strip_extension(file_name):
    # Split the file name into the name and the extension
    name, extension = os.path.splitext(file_name)

    print(f"File name: {name}")
    return name

def mol2_to_pdb_and_back(input_file, output_file, residue_name, cwd=None):

    if cwd:
        os.chdir(cwd)
        
    obConversion = openbabel.OBConversion()
    
    # Convert MOL2 to PDB
    obConversion.SetInAndOutFormats("mol2", "pdb")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_file)
    
    # Write to PDB format (in memory)
    pdb_string = obConversion.WriteString(mol)
    
    # Read PDB string back into OBMol
    obConversion.SetInAndOutFormats("pdb", "mol2")
    mol2 = openbabel.OBMol()
    obConversion.ReadString(mol2, pdb_string)
    
    # Set the residue name
    for residue in openbabel.OBResidueIter(mol2):
        residue.SetName(residue_name)
    
    # Write back to MOL2 format
    obConversion.WriteFile(mol2, output_file)

def pdb_to_mol2(input_file, output_file, residue_name, cwd=None):

    if cwd:
        os.chdir(cwd)
    
    obConversion = openbabel.OBConversion()
    
    # Convert PDB to MOL2
    obConversion.SetInAndOutFormats("pdb", "mol2")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_file)
    print("Converting ", input_file, " to ", output_file, "...")
    
    # Set the residue name
    for residue in openbabel.OBResidueIter(mol):
        residue.SetName(residue_name)
    print("Changed Residue Name to:", residue_name)
    
    # Write to MOL2 format
    obConversion.WriteFile(mol, output_file)

def mol2_to_pdb_with_bonds(input_file, existing_pdb_file, cwd=None):
    """
    Convert a MOL2 file to a PDB format, extract the bond information, 
    and insert it into an existing PDB file. Inserts bond info before ENDMDL if present;
    otherwise, inserts before END. Adds a REMARK line documenting the bond information addition.

    Args:
        input_file (str): Path to the input MOL2 file.
        existing_pdb_file (str): Path to the existing PDB file to append bond information.
        cwd (str, optional): Working directory to change to before running. Defaults to None.
    """
    if cwd:
        os.chdir(cwd)

    obConversion = openbabel.OBConversion()

    # Convert MOL2 to PDB
    obConversion.SetInAndOutFormats("mol2", "pdb")
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, input_file)

    # Write to a temporary PDB file to hold bond information
    temp_pdb_path = "temp_with_bonds.pdb"
    obConversion.WriteFile(mol, temp_pdb_path)

    # Extract bond information from temporary PDB file
    bond_info = []
    with open(temp_pdb_path, 'r') as pdb_temp:
        for line in pdb_temp:
            if line.startswith("CONECT"):
                bond_info.append(line)

    # Clean up temporary file
    os.remove(temp_pdb_path)

    # Create the REMARK line to document the added bond information
    remark_line = f"REMARK   1 EDITED BY electrofit: Added bond information from {os.path.basename(input_file)}\n"

    # Read and modify the existing PDB content
    modified_pdb_content = [remark_line]  # Add the REMARK line at the beginning
    endmdl_found = False
    with open(existing_pdb_file, 'r') as existing_pdb:
        for line in existing_pdb:
            # Insert bond info right before ENDMDL or END, based on the first found
            if not endmdl_found and line.strip() == "ENDMDL":
                modified_pdb_content.extend(bond_info)
                endmdl_found = True
            elif line.strip() == "END" and not endmdl_found:
                modified_pdb_content.extend(bond_info)
                endmdl_found = True
            modified_pdb_content.append(line)  # Add current line

    # Write the modified content back to the PDB file
    with open(existing_pdb_file, 'w') as existing_pdb:
        existing_pdb.writelines(modified_pdb_content)

    print(f"Bond information inserted before ENDMDL or END in {existing_pdb_file}")
    print(f"REMARK added: {remark_line.strip()}")

def parse_mol2(mol2_file):
    """
    Parses the MOL2 file to extract atoms and bonds.

    Parameters:
    - mol2_file: Path to the MOL2 file.

    Returns:
    - atoms: Dictionary mapping atom_id to atom properties.
    - bonds: List of bonds with origin_atom_id, target_atom_id, and bond_type.
    """
    atoms = {}
    bonds = []
    with open(mol2_file, 'r') as f:
        lines = f.readlines()

    section = None
    for line in lines:
        line = line.strip()
        if line.startswith('@<TRIPOS>ATOM'):
            section = 'ATOM'
            continue
        elif line.startswith('@<TRIPOS>BOND'):
            section = 'BOND'
            continue
        elif line.startswith('@<TRIPOS>'):
            section = None
            continue

        if section == 'ATOM':
            if not line:
                continue
            parts = line.split()
            if len(parts) < 9:
                continue  # Incomplete atom line
            atom_id = int(parts[0])
            atom_name = parts[1]
            atom_type = parts[5]
            atoms[atom_id] = {
                'atom_name': atom_name,
                'atom_type': atom_type,
                'connections': []
            }
        elif section == 'BOND':
            if not line:
                continue
            parts = line.split()
            if len(parts) < 4:
                continue  # Incomplete bond line
            origin = int(parts[1])
            target = int(parts[2])
            bond_type = parts[3]
            bonds.append((origin, target, bond_type))
    return atoms, bonds

def build_connections(atoms, bonds):
    """
    Establishes connections between atoms based on bonds.

    Parameters:
    - atoms: Dictionary of atoms.
    - bonds: List of bonds.

    Returns:
    - None (modifies atoms in place).
    """
    for origin, target, bond_type in bonds:
        atoms[origin]['connections'].append({'atom_id': target, 'bond_type': bond_type})
        atoms[target]['connections'].append({'atom_id': origin, 'bond_type': bond_type})

def create_equivalence_groups(atoms):
    """
    Creates equivalence groups for oxygen atoms connected exclusively to the same phosphorus atom.

    Parameters:
    - atoms: Dictionary of atoms with connection information.

    Returns:
    - equiv_groups: Dictionary for equiv_groups.json.
    """
    equiv_groups = {}
    # Collect all oxygen atom_ids
    oxygen_atom_ids = [atom_id for atom_id, atom in atoms.items() if atom['atom_type'].startswith('O')]
    if not oxygen_atom_ids:
        print("No oxygen atoms found in the MOL2 file.")
        return {}
    min_oxygen_id = min(oxygen_atom_ids)

    # Create label map: atom_id -> label (e.g., 19 -> "O7" if min_oxygen_id=13)
    label_map = {atom_id: f"O{atom_id - min_oxygen_id +1}" for atom_id in oxygen_atom_ids}
    
    # Debug: Print label map
    print("\nLabel Map:")
    for atom_id in sorted(label_map.keys()):
        print(f"  Atom ID {atom_id}: {label_map[atom_id]}")
    print()

    for atom_id, atom in atoms.items():
        # Identify phosphorus atoms
        if not atom['atom_type'].startswith('P'):
            continue

        # Find connected oxygen atoms via single or double bonds
        connected_oxys = []
        for conn in atom['connections']:
            conn_atom = atoms.get(conn['atom_id'])
            if conn_atom and conn_atom['atom_type'].startswith('O'):
                # Check if this oxygen is bonded only to this phosphorus atom
                if len(conn_atom['connections']) == 1:
                    connected_oxys.append({
                        'atom_id': conn['atom_id'],
                        'bond_type': conn['bond_type']
                    })
                elif len(conn_atom['connections']) == 2:
                    # If oxygen is bonded to phosphorus and one more atom, exclude it
                    bonded_atoms = conn_atom['connections']
                    bonded_to_p = any(
                        other['atom_id'] == atom_id for other in bonded_atoms
                    )
                    bonded_to_others = any(
                        other['atom_id'] != atom_id for other in bonded_atoms
                    )
                    if bonded_to_p and not bonded_to_others:
                        connected_oxys.append({
                            'atom_id': conn['atom_id'],
                            'bond_type': conn['bond_type']
                        })
                else:
                    # More than two connections, exclude
                    continue

        if len(connected_oxys) < 2:
            continue  # Need at least two oxygens to form equivalence

        # Select central oxygen with the lowest atom_id
        central_oxy = min(connected_oxys, key=lambda oxy: oxy['atom_id'])

        # Generate unique labels based on atom_id
        central_label = label_map[central_oxy['atom_id']]
        equiv_labels = [
            label_map[oxy['atom_id']]
            for oxy in connected_oxys
            if oxy['atom_id'] != central_oxy['atom_id']
        ]

        # Debugging Output
        print(f"Phosphorus Atom ID: {atom_id}")
        print(f"  Central Oxygen: {central_label}")
        print(f"  Equivalent Oxygens: {equiv_labels}\n")

        # Ensure no duplicates and valid labels
        if equiv_labels:
            equiv_groups[central_label] = equiv_labels

    return equiv_groups

def write_equivalence_groups(equiv_groups, output_file):
    """
    Writes the equivalence groups to a JSON file.

    Parameters:
    - equiv_groups: Dictionary of equivalence groups.
    - output_file: Path to the output JSON file.

    Returns:
    - None
    """
    with open(output_file, 'w') as f:
        json.dump(equiv_groups, f, indent=4)
    print(f"Equivalence groups successfully written to '{output_file}'.")

def generate_atom_labels(atom_numbers):
    element_symbols = {
        1: 'H', 6: 'C', 7: 'N', 8: 'O', 15: 'P',
        16: 'S', 17: 'Cl', 35: 'Br', 53: 'I',
    }
    label_counters = {}
    atom_labels = []
    for z in atom_numbers:
        symbol = element_symbols.get(z, f'Z{z}')
        count = label_counters.get(symbol, 0) + 1
        label_counters[symbol] = count
        atom_labels.append(f"{symbol}{count}")
    return atom_labels

def load_equivalence_groups(equiv_groups_file):
    try:
        with open(equiv_groups_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading equivalence groups: {e}")
        sys.exit(1)

def edit_resp_input(input_file, equiv_groups_file, output_file, ignore_sym=False):
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error reading RESP input file: {e}")
        sys.exit(1)

    charge_mult_pattern = re.compile(r'^\s*-?\d+\s+(\d+)\s*$')
    atom_data_start = None
    num_atoms = None

    for idx, line in enumerate(lines):
        if 'Resp charges for organic molecule' in line and '1.0' in lines[idx-1]:
            charge_mult_line = lines[idx + 1].strip()
            match = charge_mult_pattern.match(charge_mult_line)
            if match:
                num_atoms = int(match.group(1))
                atom_data_start = idx + 2
            else:
                print(f"Error parsing charge and multiplicity line: '{charge_mult_line}'")
                sys.exit(1)
            break
    else:
        print("Error locating atom data in RESP input file.")
        sys.exit(1)

    atom_numbers = []
    for i in range(num_atoms):
        line = lines[atom_data_start + i].strip()
        parts = line.split()
        try:
            atom_numbers.append(int(parts[0]))
        except (IndexError, ValueError):
            atom_numbers.append(0)

    atom_labels = generate_atom_labels(atom_numbers)
    label_to_index = {label: idx + 1 for idx, label in enumerate(atom_labels)}

    equivalence_groups_labels = load_equivalence_groups(equiv_groups_file)
    
    if ignore_sym:
        equivalence_groups_labels = {}

    equivalence_groups = {}
    for central_label, equiv_labels in equivalence_groups_labels.items():
        central_idx = label_to_index.get(central_label)
        if not central_idx:
            print(f"Warning: Central atom '{central_label}' not found.")
            continue
        equiv_indices = [label_to_index.get(lbl) for lbl in equiv_labels if label_to_index.get(lbl)]
        if equiv_indices:
            equivalence_groups[central_idx] = equiv_indices

    for i in range(num_atoms):
        line_number = atom_data_start + i
        parts = lines[line_number].strip().split()
        try:
            z_int = int(parts[0])
            lines[line_number] = f"  {z_int:>3}    0\n"
        except (IndexError, ValueError):
            continue

    for central_idx, equiv_indices in equivalence_groups.items():
        line_number = atom_data_start + (central_idx - 1)
        try:
            z_int = int(lines[line_number].strip().split()[0])
            lines[line_number] = f"  {z_int:>3}    0\n"
        except (IndexError, ValueError):
            continue
        for equiv_idx in equiv_indices:
            line_number = atom_data_start + (equiv_idx - 1)
            try:
                z_int = int(lines[line_number].strip().split()[0])
                lines[line_number] = f"  {z_int:>3}  {central_idx:>3}\n"
            except (IndexError, ValueError):
                continue

    try:
        with open(output_file, 'w') as f:
            f.writelines(lines)
        print(f"Modified RESP input file saved as '{output_file}'.")
    except Exception as e:
        print(f"Error writing output file: {e}")
        sys.exit(1)

def modify_gaussian_input(file_path):
    """
    Modifies the Gaussian input file to remove 'opt' and set charge/multiplicity.
    
    Parameters:
    - file_path (str): Path to the Gaussian input file.
    """
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
        
        with open(file_path, 'w') as file:
            for line in lines:
                # Remove 'opt' from the route section
                if ' opt' in line:
                    line = line.replace(' opt', '')
                    logging.info("Removed 'opt' from the route section.")
                
                file.write(line)
                
    except Exception as e:
        logging.error(f"Error modifying Gaussian input file: {e}")
        raise

def parse_charges_from_mol2(mol2_file):
    """
    Parses the MOL2 file to extract atom names and charges.
    
    Parameters:
    - mol2_file: Path to the MOL2 file.
    
    Returns:
    - atoms: Dictionary mapping atom names to a list of charges.
    """
    atoms = {}
    with open(mol2_file, 'r') as f:
        lines = f.readlines()

    section = None
    for line in lines:
        line = line.strip()
        if line.startswith('@<TRIPOS>ATOM'):
            section = 'ATOM'
            continue
        elif line.startswith('@<TRIPOS>'):
            section = None
            continue

        if section == 'ATOM':
            if not line:
                continue
            parts = line.split()
            if len(parts) < 9:
                continue  # Incomplete atom line
            atom_name = parts[1]
            atom_charge = float(parts[8]) 
            if atom_name not in atoms:
                atoms[atom_name] = {'charges': []}
            atoms[atom_name]['charges'].append(atom_charge)
        else:
            continue

    return atoms

def adjust_atom_names(atoms_dict):
    """
    Adjusts atom names in the provided dictionary by appending a count to each unique element symbol.
    
    Parameters:
        atoms_dict (dict): Dictionary where keys are atom names and values are properties associated with each atom.
    
    Returns:
        dict: New dictionary with adjusted atom names.
    """
    
    # Get the list of atom names from the dictionary keys
    atom_names = list(atoms_dict.keys())

    # Initialize a dictionary to keep track of counts for each element
    counts = {}
    adjusted_atom_names = []
    
    for name in atom_names:
        # Extract the element symbol (handles one or two-letter symbols)
        match = re.match(r'^([A-Z][a-z]?)(\d*)', name)
        if match:
            element = match.group(1)
        else:
            # If the name doesn't match the pattern, keep it unchanged
            adjusted_atom_names.append(name)
            continue

        # Update the count for the element
        counts.setdefault(element, 0)
        counts[element] += 1

        # Adjust the name with the element symbol and its count
        adjusted_name = f"{element}{counts[element]}"
        adjusted_atom_names.append(adjusted_name)

    # Create a mapping from old names to adjusted names
    name_mapping = dict(zip(atom_names, adjusted_atom_names))

    # Update atoms_dict keys with adjusted names
    adjusted_atoms_dict = {new_name: atoms_dict[old_name] for old_name, new_name in name_mapping.items()}
    
    return adjusted_atoms_dict

def extract_charges_from_subdirectories(base_dir, results_dir):
    """
    Walk through subdirectories of the base_dir and extract charges from mol2 files.
    
    Parameters:
    - base_dir: Path to the directory containing subdirectories with mol2 files.
    
    Returns:
    - adjusted_atoms_dict: Dictionary of atoms with adjusted names and charges collected from all subdirectories.
    """
    atoms_dict = {}
    # Traverse each subdirectory in 'base_dir'
    subdirs = [f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))]
    
    for subdir in subdirs:
        subdir_path = os.path.join(base_dir, subdir)
        for file_name in os.listdir(subdir_path):
            if file_name.endswith("_resp.mol2"):
                mol2_path = os.path.join(subdir_path, file_name)
                atoms = parse_charges_from_mol2(mol2_path)
                
                # If it's the first subdirectory, initialize the atom names
                if not atoms_dict:
                    # Initialize atoms_dict with atom names and empty charge lists
                    for atom_name, atom_data in atoms.items():
                        atoms_dict[atom_name] = {'charges': []}

                # Collect the charges for the atoms in the atoms_dict
                for atom_name, atom_data in atoms.items():
                    atoms_dict[atom_name]['charges'].extend(atom_data['charges'])

    # Calculate the mean charge for each atom
    for atom_name, atom_data in atoms_dict.items():
        charges = atom_data['charges']
        atom_data['average_charge'] = np.mean(charges) if charges else 0
    
    adjusted_atoms_dict = adjust_atom_names(atoms_dict)

    # Write the average charges to the output file
    output_file = os.path.join(results_dir, "average_charges.txt")
    try:
        with open(output_file, 'w') as f:
            f.write("#Atom_Name\tAverage_Charge\n")
            for atom_name, atom_data in adjusted_atoms_dict.items():
                f.write(f"{atom_name}\t{atom_data['average_charge']:.4f}\n")
        print(f"Average charges successfully written to {output_file}")
    except Exception as e:
        print(f"An error occurred while writing to {output_file}: {e}")

    return adjusted_atoms_dict

# Load symmetry groups from a JSON file
def load_symmetry_groups(json_path):
    with open(json_path, 'r') as file:
        symmetry_groups = json.load(file)
    return symmetry_groups


def replace_charge_in_ac_file(file_path, new_charge_float, cwd=None):
    """
    Replaces all CHARGE lines in the file with the new charge value.
    """
    if cwd:
        home = os.getcwd()
        os.chdir(cwd)

    new_charge_int = int(round(new_charge_float))
    new_charge_line = f"CHARGE     {new_charge_float:.2f} ( {new_charge_int} )\n"

    with open(file_path, 'r') as file:
        lines = file.readlines()

    charge_replaced = False
    updated_lines = []
    for line in lines:
        if line.startswith('CHARGE'):
            updated_lines.append(new_charge_line)
            charge_replaced = True
            print(f"Replaced line:\nOld: {line.strip()}\nNew: {new_charge_line.strip()}")
        else:
            updated_lines.append(line)

    if not charge_replaced:
        print("No CHARGE line found in the file.")
        return

    with open(file_path, 'w') as file:
        file.writelines(updated_lines)
    
    if cwd:
        os.chdir(home)

    print(f"All CHARGE entries have been updated to {new_charge_float:.2f} ({new_charge_int}).")