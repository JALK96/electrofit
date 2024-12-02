import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np
import matplotlib.pyplot as plt
import os
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

            # Define the list of dihedrals with atom names
            dihedral_list = [
                ("C", "C1", "C2", "C3"),   # Dihedral 1
                ("C1", "C2", "C3", "C4"),  # Dihedral 2
                ("C2", "C3", "C4", "C5"),  # Dihedral 3
                ("C3", "C4", "C5", "C"),   # Dihedral 4
                ("C4", "C5", "C", "C1"),   # Dihedral 5
                ("C5", "C", "C1", "C2"),   # Dihedral 6
            ]

            # Collect atom indices for each dihedral
            dihedral_indices = []
            for dihedral in dihedral_list:
                indices = []
                for name in dihedral:
                    sel = u.select_atoms(f"name {name}")
                    if len(sel) != 1:
                        raise ValueError(f"Atom selection for name '{name}' did not return exactly one atom.")
                    indices.append(sel[0].index)
                dihedral_indices.append(indices)

            print(dihedral_indices)

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
            dihedral_angles = np.array(dihedral_angles)  # Shape: (num_frames, num_dihedrals)

            # Convert radians to degrees
            angles_deg = np.rad2deg(dihedral_angles)
            colors_list = [
                'darkblue', 'royalblue', 'darkred',
                'green', 'darkorange', 'gold'
            ]
            # Plot all dihedrals in one figure
            plt.figure(figsize=(8, 4))
            for idx in range(angles_deg.shape[1]):
                plt.plot(times, angles_deg[:, idx], label=f'τ {idx+1}', color=colors_list[idx], alpha=1-idx/10)
            plt.xlabel('Time (ps)')
            plt.xlim((min(times), max(times)))
            plt.ylabel('Dihedral Angle (degrees)')
            plt.title('Dihedral Angles Over Time')
            plt.legend(loc='center right')
            plt.savefig('dihedrals_C_backbone.pdf')
            plt.close()