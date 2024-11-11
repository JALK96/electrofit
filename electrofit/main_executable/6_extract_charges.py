import os
import numpy as np
import matplotlib.pyplot as plt
import sys

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


sys.path.append(project_path)

from electrofit.helper.file_manipulation import find_file_with_extension, mol2_to_pdb_with_bonds
from electrofit.helper.config_parser import ConfigParser


def parse_mol2(mol2_file):
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

def extract_charges_from_subdirectories(base_dir):
    """
    Walk through subdirectories of the base_dir and extract charges from mol2 files.
    
    Parameters:
    - base_dir: Path to the directory containing subdirectories with mol2 files.
    
    Returns:
    - atoms_dict: Dictionary of atoms with charges collected from all subdirectories.
    """
    atoms_dict = {}
    # Traverse each subdirectory in 'base_dir'
    subdirs = [f for f in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, f))]
    
    for subdir in subdirs:
        subdir_path = os.path.join(base_dir, subdir)
        for file_name in os.listdir(subdir_path):
            if file_name.endswith("_resp.mol2"):
                mol2_path = os.path.join(subdir_path, file_name)
                atoms = parse_mol2(mol2_path)
                
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

    return atoms_dict

def plot_charges(atoms_dict, base_dir):
    """
    Plot the charges of atoms along with their average charges.

    Parameters:
    - atoms_dict: Dictionary of atoms and their charges.
    """
    atom_names = list(atoms_dict.keys())
    min_charges = [min(atoms_dict[atom]['charges']) for atom in atom_names]
    max_charges = [max(atoms_dict[atom]['charges']) for atom in atom_names]
    avg_charges = [atoms_dict[atom]['average_charge'] for atom in atom_names]
    avg_sum = sum(avg_charges)



    # Plotting
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Draw vertical lines showing charge range for each atom
    for i, atom_name in enumerate(atom_names):
        ax.plot([i, i], [min_charges[i], max_charges[i]], color="darkblue", linewidth=2)
    
    # Add red dots for the average charge
    ax.scatter(range(len(atom_names)), avg_charges, color="darkblue", marker='^', label="Average Charge", zorder=5)

    # Labeling
    ax.set_xticks(range(len(atom_names)))
    ax.set_xticklabels(atom_names, rotation=90)
    ax.set_xlabel('Atom Names')
    ax.set_ylabel('Atomic Partial Charge (e)')
    ax.set_title(f'Charge Range and Mean for Each Atom. Average Sum {round(avg_sum,0)}')
    ax.legend(["Charge Range", "Average Charge"])

    # Save and show the plot
    plt.tight_layout()
    figure_path = os.path.join(base_dir, "charges.pdf")
    plt.savefig(figure_path)
    #plt.show()





process_dir = os.path.join(project_path, "process")


# Iterate over each subdirectory in the process directory
for sub_dir in os.listdir(process_dir):

    # Base directory containing the subdirectories with mol2 files
    base_dir = os.path.join(process_dir, sub_dir, "extracted_conforms")

    # Extract charges from all subdirectories
    atoms_dict = extract_charges_from_subdirectories(base_dir)

    # Plot the charges and average charges
    plot_charges(atoms_dict, base_dir)