import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import re

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

def shift_atom_number(name):
    match = re.match(r"([A-Za-z]+)(\d+)?", name)
    if match:
        base = match.group(1)
        num = match.group(2)
        if num is not None:
            new_num = int(num) + 1
            return f"{base}{new_num}"
        else:
            # For names without numbers, add '1'
            return f"{base}1"
    else:
        # Return the name as is if it doesn't match the pattern
        return name
    
script_dir = os.path.dirname(os.path.abspath(__file__))
project_path = find_project_root(current_dir=script_dir)

sys.path.append(project_path)

# Define the base process directory
process_dir = os.path.join(project_path, "process.nobackup")

# Loop through each subdirectory in the process directory
for folder_name in os.listdir(process_dir):
    folder_path = os.path.join(process_dir, folder_name)
    
    # Check if it's a directory
    if os.path.isdir(folder_path):
        # Define the 'run_final_gmx_simulation' directory within this folder
        run_final_sim_dir = os.path.join(folder_path, 'run_final_gmx_simulation')
        
        # Check if 'run_final_gmx_simulation' exists
        if os.path.isdir(run_final_sim_dir):
            # Define the destination directory 'analyze_final_sim' within the same working directory
            dest_dir = os.path.join(folder_path, 'analyze_final_sim')
            os.makedirs(dest_dir, exist_ok=True)
            os.chdir(dest_dir)

            gro_file_path = os.path.join(run_final_sim_dir, 'md.gro')
            xtc_file_path = os.path.join(run_final_sim_dir, 'md_center.xtc')

            # Load Universe
            u = mda.Universe(gro_file_path, xtc_file_path)

            # Define dihedral groups
            dihedral_groups = [
                {
                    'name': 'Ring Dihedrals',
                    'dihedrals': [
                        ("C", "C1", "C2", "C3"),   # Dihedral 1
                        ("C1", "C2", "C3", "C4"),  # Dihedral 2
                        ("C2", "C3", "C4", "C5"),  # Dihedral 3
                        ("C3", "C4", "C5", "C"),   # Dihedral 4
                        ("C4", "C5", "C", "C1"),   # Dihedral 5
                        ("C5", "C", "C1", "C2"),   # Dihedral 6
                    ]
                },
                {
                    'name': 'Dihedrals for C1',
                    'dihedrals': [
                        ("C", "O", "P", "O6"),       # C1-O1-P1-O7
                        ("C", "O", "P", "O7"),       # C1-O1-P1-O8
                        ("C", "O", "P", "O8"),       # C1-O1-P1-O9
                    ]
                },
                {
                    'name': 'Dihedrals for C2',
                    'dihedrals': [
                        ("C1", "O1", "P1", "O9"),    # C2-O2-P2-O10
                        ("C1", "O1", "P1", "O10"),   # C2-O2-P2-O11
                        ("C1", "O1", "P1", "O11"),   # C2-O2-P2-O12
                    ]
                },
                {
                    'name': 'Dihedrals for C3',
                    'dihedrals': [
                        ("C2", "O2", "P2", "O12"),   # C3-O3-P3-O13
                        ("C2", "O2", "P2", "O13"),   # C3-O3-P3-O14
                        ("C2", "O2", "P2", "O14"),   # C3-O3-P3-O15
                    ]
                },
                {
                    'name': 'Dihedrals for C4',
                    'dihedrals': [
                        ("C3", "O3", "P3", "O15"),   # C4-O4-P4-O16
                        ("C3", "O3", "P3", "O16"),   # C4-O4-P4-O17
                        ("C3", "O3", "P3", "O17"),   # C4-O4-P4-O18
                    ]
                },
                {
                    'name': 'Dihedrals for C5',
                    'dihedrals': [
                        ("C4", "O4", "P4", "O18"),   # C5-O5-P5-O19
                        ("C4", "O4", "P4", "O19"),   # C5-O5-P5-O20
                        ("C4", "O4", "P4", "O20"),   # C5-O5-P5-O21
                    ]
                },
                {
                    'name': 'Dihedrals for C6',
                    'dihedrals': [
                        ("C5", "O5", "P5", "O21"),   # C6-O6-P6-O22
                        ("C5", "O5", "P5", "O22"),   # C6-O6-P6-O23
                        ("C5", "O5", "P5", "O23"),   # C6-O6-P6-O24
                    ]
                },
            ]

            # Collect all dihedral definitions and group indices
            dihedral_definitions = []
            dihedral_group_indices = []
            group_names = []

            for group_idx, group in enumerate(dihedral_groups):
                group_names.append(group['name'])
                for dihedral in group['dihedrals']:
                    dihedral_definitions.append(dihedral)
                    dihedral_group_indices.append(group_idx)

            # Collect atom indices for each dihedral
            dihedral_indices = []
            for dihedral in dihedral_definitions:
                indices = []
                for name in dihedral:
                    sel = u.select_atoms(f"name {name}")
                    if len(sel) != 1:
                        raise ValueError(f"Atom selection for name '{name}' did not return exactly one atom.")
                    indices.append(sel[0].index)
                dihedral_indices.append(indices)

            print("Dihedral Indices:", dihedral_indices)

            # Initialize lists to store times and dihedral angles
            times = []
            dihedral_angles = []

            # Loop over trajectory frames
            for ts in u.trajectory:
                positions = u.atoms.positions
                box = u.dimensions
                dihedrals_this_frame = []
                for indices in dihedral_indices:
                    A_pos = positions[indices[0]]
                    B_pos = positions[indices[1]]
                    C_pos = positions[indices[2]]
                    D_pos = positions[indices[3]]
                    dihedral = calc_dihedrals(
                        A_pos[np.newaxis, :], B_pos[np.newaxis, :],
                        C_pos[np.newaxis, :], D_pos[np.newaxis, :],
                        box=box
                    )
                    dihedrals_this_frame.append(dihedral[0])
                dihedral_angles.append(dihedrals_this_frame)
                times.append(ts.time)

            # Convert lists to numpy arrays
            times = np.array(times) 
            times = times / 1000 # in ns
            dihedral_angles = np.array(dihedral_angles)  # Shape: (num_frames, num_dihedrals)

            # Convert radians to degrees
            angles_deg = np.rad2deg(dihedral_angles)

            # Organize dihedral indices per group
            num_groups = len(dihedral_groups)
            group_dihedral_indices = [[] for _ in range(num_groups)]

            for dihedral_idx, group_idx in enumerate(dihedral_group_indices):
                group_dihedral_indices[group_idx].append(dihedral_idx)

            # Plot dihedrals per group
            for group_idx, group in enumerate(dihedral_groups):
                dihedral_idxs = group_dihedral_indices[group_idx]
                num_dihedrals_in_group = len(dihedral_idxs)

                if group['name'].startswith('Dihedrals for C'):
                    # Create subplots for groups 'Dihedrals for C1' to 'Dihedrals for C6'
                    num_subplots = num_dihedrals_in_group
                    fig, axes = plt.subplots(nrows=num_subplots, ncols=1, figsize=(8, 2*num_subplots), sharex=True)
                    # Ensure axes is iterable
                    if num_subplots == 1:
                        axes = [axes]
                    for idx_in_group, dihedral_idx in enumerate(dihedral_idxs):
                        ax = axes[idx_in_group]
                        # Get the color from the colors list
                        colors_list = ['darkblue', 'darkgreen', 'darkred']
                        color = colors_list[idx_in_group % len(colors_list)]
                        # Create the label from the dihedral definition
                        dihedral_atoms = dihedral_definitions[dihedral_idx]
                        dihedral_atoms_shifted = [shift_atom_number(atom) for atom in dihedral_atoms]
                        label = '-'.join(dihedral_atoms_shifted)
                        label = label.replace(' ', '_')  # Replace spaces if any
                        ax.plot(times, angles_deg[:, dihedral_idx], label=label, color=color)
                        if idx_in_group == 1: # Second subplot
                            ax.set_ylabel('Dihedral Angle (degrees)')
                        #ax.set_title(f'{group["name"]}: {label}')
                        ax.legend(loc='best')
                    axes[-1].set_xlabel('Time (ns)')
                    plt.tight_layout()
                    # Save the plot
                    plot_filename = f'dihedrals_{group_idx+1}_{group["name"].replace(" ", "_")}.pdf'
                    plt.savefig(plot_filename)
                    plt.close()
                else:
                    # For other groups, keep the existing plotting style
                    plt.figure(figsize=(8, 4))
                    colors_list = ['darkblue', 'royalblue', 'darkred', 'green', 'darkorange', 'gold']
                    for idx_in_group, dihedral_idx in enumerate(dihedral_idxs):
                        plt.plot(times, angles_deg[:, dihedral_idx], label=f'τ {idx_in_group+1}', color=colors_list[idx_in_group], alpha=1-idx_in_group/10)
                    plt.xlabel('Time (ns)')
                    plt.xlim((min(times), max(times)))
                    plt.ylabel('Dihedral Angle (degrees)')
                    plt.title(f'Dihedral Angles Over Time - {group["name"]}')
                    plt.legend(loc='best')
                    plt.tight_layout()
                    # Save the plot
                    plot_filename = f'dihedrals_{group_idx+1}_{group["name"].replace(" ", "_")}.pdf'
                    plt.savefig(plot_filename)
                    plt.close()