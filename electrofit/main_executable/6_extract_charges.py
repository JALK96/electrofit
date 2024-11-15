import os
import numpy as np
import matplotlib.pyplot as plt
import sys
import re
import shutil
import glob
import fnmatch
import matplotlib.colors as mcolors
import json


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

# Example usage:
# atoms_dict = {"H1": {}, "O1": {}, "H2": {}, "O2": {}}
# adjusted_atoms_dict = adjust_atom_names(atoms_dict)
# print(adjusted_atoms_dict)

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

def plot_charges_old(atoms_dict, initial_charges_dict, base_dir):
    """
    Plot the charges of atoms along with their average charges.

    Parameters:
    - atoms_dict: Dictionary of atoms and their charges.
    """
    atom_names = list(atoms_dict.keys())
    min_charges = [min(atoms_dict[atom]['charges']) for atom in atom_names]
    max_charges = [max(atoms_dict[atom]['charges']) for atom in atom_names]
    avg_charges = [atoms_dict[atom]['average_charge'] for atom in atom_names]
    init_charges = [initial_charges_dict[atom]['charges'] for atom in atom_names]
    print(init_charges)
    avg_sum = sum(avg_charges)



    # Plotting
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Draw vertical lines showing charge range for each atom
    for i, atom_name in enumerate(atom_names):
        ax.plot([i, i], [min_charges[i], max_charges[i]], color="darkblue", linewidth=2)
    
    # Add blue triangels for the average charge
    ax.scatter(range(len(atom_names)), avg_charges, color="darkblue", marker='^', label="Average Charge", zorder=5)

    # Add red dots for the initial charge
    ax.scatter(range(len(atom_names)), init_charges, color="darkred", marker='o', label="Average Charge", zorder=5)

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
    charges_path = os.path.join(base_dir, "average_charges.chg")
    with open(charges_path, "w") as output:
        for i in avg_charges:
            output.write(str(round(i,4)) + '\n')
    #plt.show()


def plot_charges_old2(atoms_dict, initial_charges_dict, base_dir):
    """
    Plot the charge distributions of atoms along with their average and initial charges.

    Parameters:
    - atoms_dict: Dictionary of atoms and their charges.
    - initial_charges_dict: Dictionary of initial charges for each atom.
    - base_dir: Directory to save the plot and charges.
    """
    atom_names = list(atoms_dict.keys())
    charges_data = [atoms_dict[atom]['charges'] for atom in atom_names]
    avg_charges = [atoms_dict[atom]['average_charge'] for atom in atom_names]
    init_charges = [initial_charges_dict[atom]['charges'] for atom in atom_names]
    avg_sum = sum(avg_charges)

    # Plotting
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot violin plots for each atom's charge distribution
    vp = ax.violinplot(charges_data, positions=range(len(atom_names)), widths=0.6, showmeans=True, showmedians=False,)
    color= "blue" # Hydrogen (all atoms with O in its name)
    color= "red" # Oxygen 
    color= "orange" # Phosphor atoms 
    color= "grey" # Carbon
    for body in vp['bodies']:
        body.set_facecolor(f'dark{color}')
        body.set_edgecolor(f'dark{color}')
        body.set_alpha(0.7)
    vp['cmeans'].set_edgecolor(f'dark{color}')
    vp['cmeans'].set_linewidth(2)

    # Set the color of the vertical lines (sticks) and min/max markers (caps) to dark blue
    for partname in ('cbars', 'cmins', 'cmaxes'):
        if partname in vp:
            vp[partname].set_edgecolor(f'dark{color}')
            vp[partname].set_linewidth(1.5)
            if partname == 'cbars':
                vp[partname].set_alpha(0.7)
    

    
    # Add blue triangles for the average charge
    #ax.scatter(range(len(atom_names)), avg_charges, color="darkblue", marker='^', label="Average Charge", zorder=5, alpha=0.5)

    # Add red dots for the initial charge
    ax.scatter(range(len(atom_names)), init_charges, color="red", marker='o', label="Initial Charge", zorder=5, alpha=1, s=5)

    # Labeling
    ax.set_xticks(range(len(atom_names)))
    ax.set_xticklabels(atom_names, rotation=90)
    ax.set_xlabel('Atom Names')
    ax.set_ylabel('Atomic Partial Charge (e)')
    ax.set_title(f'Charge Distribution for Each Atom (Average Sum: {round(avg_sum, 4)})')
    ax.legend()

    # Save and close the plot
    plt.tight_layout()
    figure_path = os.path.join(base_dir, "charges.pdf")
    plt.savefig(figure_path)
    plt.close(fig)  # Close the figure to free up memory

    # Save average charges to a file
    charges_path = os.path.join(base_dir, "average_charges.chg")
    with open(charges_path, "w") as output:
        for i in avg_charges:
            output.write(str(round(i,4)) + '\n')


def plot_charges_by_atom(atoms_dict, initial_charges_dict, base_dir):
    """
    Plot the charge distributions of atoms along with their average and initial charges,
    coloring the violin plots based on atom types.

    Parameters:
    - atoms_dict: Dictionary of atoms and their charges.
    - initial_charges_dict: Dictionary of initial charges for each atom.
    - base_dir: Directory to save the plot and charges.
    """
    atom_names = list(atoms_dict.keys())
    charges_data = [atoms_dict[atom]['charges'] for atom in atom_names]
    avg_charges = [atoms_dict[atom]['average_charge'] for atom in atom_names]
    init_charges = [initial_charges_dict[atom]['charges'] for atom in atom_names]
    avg_sum = sum(avg_charges)

    # Define a color mapping for each atom type
    color_map = {
        'H': 'royalblue',     # Hydrogen
        'O': 'darkred',      # Oxygen
        'P': 'darkorange',   # Phosphorus
        'C': 'darkblue',     # Carbon
        # Add more elements and colors as needed
    }

    # Plotting
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot violin plots for each atom's charge distribution
    vp = ax.violinplot(
        charges_data,
        positions=range(len(atom_names)),
        widths=0.6,
        showmeans=True,
        showmedians=False,
        showextrema=True,
    )

    # Iterate over each violin plot and set colors based on atom type
    for i, body in enumerate(vp['bodies']):
        atom_name = atom_names[i]

        # Extract the element symbol from the atom name
        match = re.match(r'^([A-Z][a-z]?)(\d*)', atom_name)
        if match:
            element = match.group(1)
        else:
            element = 'Unknown'  # Default value if no match is found

        # Get the color for this element from the color map
        color = color_map.get(element, 'black')  # Default to black if element not in color_map

        # Set the face and edge color of the violin body
        body.set_facecolor(color)
        body.set_edgecolor(color)

    # Customize the mean lines
    if 'cmeans' in vp:
        means = vp['cmeans']
        # Prepare a list of colors for the mean lines
        mean_colors = []
        for i in range(len(atom_names)):
            atom_name = atom_names[i]
            match = re.match(r'^([A-Z][a-z]?)(\d*)', atom_name)
            if match:
                element = match.group(1)
            else:
                element = 'Unknown'
            color = color_map.get(element, 'black')
            mean_colors.append(color)
        # Set the colors and linewidths for the mean lines
        vp['cmeans'].set_color(mean_colors)
        vp['cmeans'].set_linewidth(2)

    # Set the color of the vertical lines and caps
    for partname in ('cbars', 'cmins', 'cmaxes'):
        if partname in vp:
            items = vp[partname]
            if isinstance(items, list):
                # For older versions of Matplotlib where items are lists
                for i, item in enumerate(items):
                    atom_name = atom_names[i]
                    match = re.match(r'^([A-Z][a-z]?)(\d*)', atom_name)
                    if match:
                        element = match.group(1)
                    else:
                        element = 'Unknown'
                    color = color_map.get(element, 'black')
                    item.set_edgecolor(color)
                    item.set_linewidth(1.5)
                    if partname == 'cbars':
                        item.set_alpha(0.7)
            else:
                # For newer versions where items is a LineCollection
                # Prepare a list of colors
                line_colors = []
                for i in range(len(atom_names)):
                    atom_name = atom_names[i]
                    match = re.match(r'^([A-Z][a-z]?)(\d*)', atom_name)
                    if match:
                        element = match.group(1)
                    else:
                        element = 'Unknown'
                    color = color_map.get(element, 'black')
                    line_colors.append(color)
                items.set_color(line_colors)
                items.set_linewidth(1.5)

    # Add scatter points for average and initial charges
    #ax.scatter(range(len(atom_names)), avg_charges, color="black", marker='^', label="Average Charge", zorder=5)
    ax.scatter(range(len(atom_names)), init_charges, color="black", marker='o', label="Initial Charge", zorder=5, s=5)

    # Labeling
    ax.set_xticks(range(len(atom_names)))
    ax.set_xticklabels(atom_names, rotation=90)
    ax.set_xlabel('Atom Names')
    ax.set_ylabel('Atomic Partial Charge (e)')
    ax.set_title(f'Charge Distribution for Each Atom (Average Sum: {round(avg_sum, 4)})')

    # Create custom legend
    handles = [
        plt.Line2D([], [], color='royalblue', marker='None', linestyle='-', markersize=10, label='Hydrogen'),
        plt.Line2D([], [], color='darkred', marker='None', linestyle='-', markersize=10, label='Oxygen'),
        plt.Line2D([], [], color='darkorange', marker='None', linestyle='-', markersize=10, label='Phosphorus'),
        plt.Line2D([], [], color='darkblue', marker='None', linestyle='-', markersize=10, label='Carbon'),

        plt.Line2D([], [], color='black', marker='o', linestyle='None', markersize=5, label='Initial Charge'),
    ]
    ax.legend(handles=handles, title='Average Charges', frameon=False)

    # Save and close the plot
    plt.tight_layout()
    figure_path = os.path.join(base_dir, "charges.pdf")
    plt.savefig(figure_path)
    plt.close(fig)  # Close the figure to free up memory

    # Save average charges to a file
    charges_path = os.path.join(base_dir, "average_charges.chg")
    with open(charges_path, "w") as output:
        for i in avg_charges:
            output.write(str(round(i,4)) + '\n')



def create_atom_color_mapping(atom_names, symmetry_groups):
    # List of colors to use for the groups
    group_colors_list = ['darkred', 'darkgreen', 'darkorange', 'purple', 'royalblue', 'lightcoral', 'deepskyblue', 'mediumvioletred']
    group_colors = {}

    # Map each group to a color
    for i, (group_representative, group_atoms) in enumerate(symmetry_groups.items()):
        color = group_colors_list[i % len(group_colors_list)]  # Cycle through colors if needed
        # Include the representative atom in the group
        group = [group_representative] + group_atoms
        for atom in group:
            group_colors[atom] = color

    # Assign 'darkblue' to atoms not in any group
    atom_to_color = {}
    for atom in atom_names:
        color = group_colors.get(atom, 'darkblue')
        atom_to_color[atom] = color

    return atom_to_color

def plot_charges_by_symmetry(atoms_dict, initial_charges_dict, base_dir, symmetry_groups):
    """
    Plot the charge distributions of atoms, coloring the violin plots based on symmetry groups.

    Parameters:
    - atoms_dict: Dictionary of atoms and their charges.
    - initial_charges_dict: Dictionary of initial charges for each atom.
    - base_dir: Directory to save the plot and charges.
    - symmetry_groups: Dictionary mapping representative atoms to lists of equivalent atoms.
    """
    import os
    import matplotlib.pyplot as plt

    atom_names = list(atoms_dict.keys())
    charges_data = [atoms_dict[atom]['charges'] for atom in atom_names]
    avg_charges = [atoms_dict[atom]['average_charge'] for atom in atom_names]
    init_charges = [initial_charges_dict[atom]['charges'] for atom in atom_names]
    avg_sum = sum(avg_charges)

    # Create atom-to-color mapping
    atom_to_color = create_atom_color_mapping(atom_names, symmetry_groups)

    # Plotting
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot violin plots for each atom's charge distribution
    vp = ax.violinplot(
        charges_data,
        positions=range(len(atom_names)),
        widths=0.6,
        showmeans=True,
        showmedians=False,
        showextrema=True,
    )

    # Customize the violin plots
    for i, body in enumerate(vp['bodies']):
        atom_name = atom_names[i]
        color = atom_to_color[atom_name]
        body.set_facecolor(color)
        body.set_edgecolor(color)
        #body.set_alpha(0.7)

    # Customize the mean lines
    if 'cmeans' in vp:
        mean_colors = [atom_to_color[atom_names[i]] for i in range(len(atom_names))]
        vp['cmeans'].set_color(mean_colors)
        vp['cmeans'].set_linewidth(2)

    # Set the color of the vertical lines and min/max markers (caps)
    for partname in ('cbars', 'cmins', 'cmaxes'):
        if partname in vp:
            items = vp[partname]
            colors = [atom_to_color[atom_names[i]] for i in range(len(atom_names))]
            items.set_color(colors)
            items.set_linewidth(1.5)

    # Add red dots for the initial charge
    ax.scatter(
        range(len(atom_names)),
        init_charges,
        color="black",
        marker='o',
        label="Initial Charge",
        zorder=5,
        alpha=1,
        s=5
    )

    # Labeling
    ax.set_xticks(range(len(atom_names)))
    ax.set_xticklabels(atom_names, rotation=90)
    ax.set_xlabel('Atom Names')
    ax.set_ylabel('Atomic Partial Charge (e)')
    ax.set_title(f'Charge Distribution for Each Atom (Average Sum: {round(avg_sum, 4)})')
    ax.legend(frameon=False)

    # Save and close the plot
    plt.tight_layout()
    figure_path = os.path.join(base_dir, "charges_by_symmetry.pdf")
    plt.savefig(figure_path)
    plt.close(fig)  # Close the figure to free up memory

    # Save average charges to a file
    charges_path = os.path.join(base_dir, "average_charges.chg")
    with open(charges_path, "w") as output:
        for i in avg_charges:
            output.write(str(round(i, 4)) + '\n')



# Load symmetry groups from a JSON file
def load_symmetry_groups(json_path):
    with open(json_path, 'r') as file:
        symmetry_groups = json.load(file)
    return symmetry_groups



script_dir = os.path.dirname(os.path.abspath(__file__))
project_path = find_project_root(current_dir=script_dir)
sys.path.append(project_path)
process_dir = os.path.join(project_path, "process")


from electrofit.commands.run_commands import run_command, run_acpype
from electrofit.helper.file_manipulation import find_file_with_extension
from electrofit.helper.config_parser import ConfigParser
from electrofit.helper.set_logging import setup_logging, reset_logging

# Iterate over each subdirectory in the process directory
for sub_dir in os.listdir(process_dir):

    # Base directory containing the subdirectories with mol2 files
    base_dir = os.path.join(process_dir, sub_dir, "extracted_conforms")

    # Initial processing directory 
    pis_dir = os.path.join(process_dir, sub_dir, "run_gau_create_gmx_in")

    # Acpype directory
    abstract_ac_dir = os.path.join(pis_dir, "*.acpype")


    # Use glob to find all matching files or directories
    ac_files = glob.glob(abstract_ac_dir)

    ac_dir = ac_files[0]

    results_dir = os.path.join(process_dir, sub_dir, "results")
    os.makedirs(results_dir, exist_ok=True)

    log_file_path = os.path.join(results_dir, "results.log")

    

    # Change to base directory containig configuration file (.ef) and copy it into the results directory
    os.chdir(base_dir)
    config_file_path = os.path.join(base_dir, find_file_with_extension("ef"))
    shutil.copy2(config_file_path, results_dir)


    # Define molecule parameters
    config = ConfigParser(config_file_path)

    charge = config.Charge
    atom_type = config.AtomType
    molecule_name = config.MoleculeName
    adjust_sym = config.AdjustSymmetry
    print(adjust_sym)

    updated_mol2_file = os.path.join(results_dir, f"averaged_{molecule_name}.mol2")

    
    mol2_file_pattern = f"*{atom_type}.mol2"


    for file_name in os.listdir(ac_dir):
        if fnmatch.fnmatch(file_name, mol2_file_pattern):
            mol2_source_file_path = os.path.join(ac_dir, file_name)
            mol2_file_name = file_name
            # Copy the file to the destination directory
            #shutil.copy(source_file_path, results_dir)
            #print(f"Copied {file_name} to {results_dir}")

    # Extract the initial charges
    initial_charges_dict = adjust_atom_names(parse_mol2(mol2_source_file_path))
  

    # Extract charges from all subdirectories
    atoms_dict = extract_charges_from_subdirectories(base_dir, results_dir)

    if adjust_sym == True:
        os.chdir(pis_dir)
        equiv_group_path = os.path.join(pis_dir, find_file_with_extension("json"))
        equiv_group = load_symmetry_groups(equiv_group_path)

        # Plot the charges and average charges
        plot_charges_by_symmetry(atoms_dict, initial_charges_dict, results_dir, equiv_group)

    plot_charges_by_atom(atoms_dict, initial_charges_dict, results_dir)
    setup_logging(log_file_path)
    average_charges = os.path.join(results_dir, "average_charges.chg")
    update_mol2 = os.path.join(project_path, "electrofit/execution/update_mol2.py")
    run_command(f"python {update_mol2} {mol2_source_file_path} {average_charges} {updated_mol2_file}")
    run_acpype(updated_mol2_file, charge, results_dir, atom_type)

    reset_logging()




