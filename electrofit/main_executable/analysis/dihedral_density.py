#!/usr/bin/env python

import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import re
import seaborn as sns
from scipy.stats import gaussian_kde

sns.set_context("talk")
sns.set_style({'font.family':'serif', 'font.serif':'Times New Roman'})

################################################################################
# HACKY CONFIG OBJECT:
# We define a small object "cf" to mimic your usage of cf.color_map and cf.dark_gray
################################################################################
import types
cf = types.SimpleNamespace()
cf.color_map = 'viridis'
cf.dark_gray = 'dimgray'

################################################################################
# FUNCTIONS
################################################################################

def find_project_root(current_dir, project_name="electrofit"):
    """
    Walk upwards in directories until we find 'project_name' or reach the filesystem root.
    """
    root = None
    while True:
        parent_dir = os.path.dirname(current_dir)
        if os.path.basename(current_dir) == project_name:
            root = current_dir  # Found the project_name directory
        if parent_dir == current_dir:
            # Reached filesystem root
            if root is None:
                raise FileNotFoundError(f"Project root directory '{project_name}' not found.")
            return root
        current_dir = parent_dir


def shift_atom_number(name):
    """
    If 'name' ends with a number, shift that number by 1. 
    Example: C1 -> C2, O10 -> O11, etc.
    """
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
        # Return the name as-is if it doesn't match the pattern
        return name

################################################################################
# MAIN SCRIPT
################################################################################

if __name__ == "__main__":
    # Get the current script directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Try to find the project root
    project_path = find_project_root(current_dir=script_dir)

    sys.path.append(project_path)

    # Define the base process directory
    process_dir = os.path.join(project_path, "process.nobackup")

    # Loop through each subdirectory in the process directory
    for folder_name in os.listdir(process_dir):
        folder_path = os.path.join(process_dir, folder_name)
        species_id = folder_name.replace("IP_", "")

        # Check if it's a directory
        if os.path.isdir(folder_path):
            # Define the 'run_final_gmx_simulation' directory within this folder
            run_final_sim_dir = os.path.join(folder_path, 'run_final_gmx_simulation')

            # Check if 'run_final_gmx_simulation' exists
            if os.path.isdir(run_final_sim_dir):
                # Define the destination directory 'analyze_final_sim'
                dest_dir = os.path.join(folder_path, 'analyze_final_sim')
                dest_dir = os.path.join(dest_dir, 'dihedral')
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
                            ("C", "O", "P", "O6"),      
                            ("C", "O", "P", "O7"),
                            ("C", "O", "P", "O8"),
                        ]
                    },
                    {
                        'name': 'Dihedrals for C2',
                        'dihedrals': [
                            ("C1", "O1", "P1", "O9"),
                            ("C1", "O1", "P1", "O10"),
                            ("C1", "O1", "P1", "O11"),
                        ]
                    },
                    {
                        'name': 'Dihedrals for C3',
                        'dihedrals': [
                            ("C2", "O2", "P2", "O12"),
                            ("C2", "O2", "P2", "O13"),
                            ("C2", "O2", "P2", "O14"),
                        ]
                    },
                    {
                        'name': 'Dihedrals for C4',
                        'dihedrals': [
                            ("C3", "O3", "P3", "O15"),
                            ("C3", "O3", "P3", "O16"),
                            ("C3", "O3", "P3", "O17"),
                        ]
                    },
                    {
                        'name': 'Dihedrals for C5',
                        'dihedrals': [
                            ("C4", "O4", "P4", "O18"),
                            ("C4", "O4", "P4", "O19"),
                            ("C4", "O4", "P4", "O20"),
                        ]
                    },
                    {
                        'name': 'Dihedrals for C6',
                        'dihedrals': [
                            ("C5", "O5", "P5", "O21"),
                            ("C5", "O5", "P5", "O22"),
                            ("C5", "O5", "P5", "O23"),
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
                times = times / 1000.0  # convert ps to ns
                dihedral_angles = np.array(dihedral_angles)  # Shape: (num_frames, num_dihedrals)

                # Convert radians to degrees
                angles_deg = np.rad2deg(dihedral_angles)

                # Organize dihedral indices per group
                num_groups = len(dihedral_groups)
                group_dihedral_indices = [[] for _ in range(num_groups)]
                for dihedral_idx, group_idx in enumerate(dihedral_group_indices):
                    group_dihedral_indices[group_idx].append(dihedral_idx)

                ############################################################################
                #  PART 1: TIME-SERIES PLOTS
                ############################################################################
                # Filter only groups whose name starts with 'Dihedrals for C'
                groups_to_plot = [group for group in dihedral_groups if group['name'].startswith('Dihedrals for C')]
                # Also record the corresponding original indices
                groups_to_plot_indices = [i for i, group in enumerate(dihedral_groups) if group['name'].startswith('Dihedrals for C')]
                num_groups = len(groups_to_plot)

                # Define a color map for however many groups we have
                colors_summary = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(groups_to_plot)))


                # Create one subplot per filtered group in one overall figure.
                fig, axes = plt.subplots(nrows=num_groups, ncols=1, figsize=(5, 2 * num_groups), sharex=True)
                if num_groups == 1:
                    axes = [axes]

                # Loop over each filtered group and plot only the first dihedral for that group.
                for idx, group in enumerate(groups_to_plot):
                    # Retrieve the original group index to access the correct dihedral indices.
                    original_group_idx = groups_to_plot_indices[idx]
                    dihedral_idxs = group_dihedral_indices[original_group_idx]
                    if len(dihedral_idxs) == 0:
                        continue

                    # Always choose the first dihedral from the group.
                    first_dihedral_idx = dihedral_idxs[0]
                    
                    # Create the label from the dihedral definition.
                    dihedral_atoms = dihedral_definitions[first_dihedral_idx]
                    dihedral_atoms_shifted = [shift_atom_number(atom) for atom in dihedral_atoms]
                    label = '-'.join(dihedral_atoms_shifted).replace(' ', '_')
                    
                    # Plot the trace for this dihedral using dark blue.
                    axes[idx].plot(times, angles_deg[:, first_dihedral_idx],
                                label=label, color=colors_summary[idx], linestyle="", marker="o")
                    
                    axes[idx].set_ylabel('Angle (°)', fontsize=18)
                    #axes[idx].legend(loc='best', fontsize=10)

                axes[-1].set_xlabel('Time (ns)', fontsize=18)
                plt.tight_layout()
                plt.savefig("dihedrals_summary.pdf")
                plt.close()

                # Individual time-series plots for each group independent
                for group_idx, group in enumerate(dihedral_groups):
                    dihedral_idxs = group_dihedral_indices[group_idx]
                    num_dihedrals_in_group = len(dihedral_idxs)

                    if group['name'].startswith('Dihedrals for C'):
                        # Create subplots for groups 'Dihedrals for C1' to 'Dihedrals for C6'
                        num_subplots = num_dihedrals_in_group
                        fig, axes = plt.subplots(
                            nrows=num_subplots, 
                            ncols=1, 
                            figsize=(7, 2 * num_subplots), 
                            sharex=True
                        )
                        # Ensure axes is iterable
                        if num_subplots == 1:
                            axes = [axes]

                        # We'll choose some distinct colors to differentiate the lines
                        colors_list = ['darkblue', 'darkgreen', 'darkred']
                        
                        for idx_in_group, dihedral_idx in enumerate(dihedral_idxs):
                            ax = axes[idx_in_group]
                            color = colors_list[idx_in_group % len(colors_list)]
                            # Create the label from the dihedral definition
                            dihedral_atoms = dihedral_definitions[dihedral_idx]
                            dihedral_atoms_shifted = [shift_atom_number(atom) for atom in dihedral_atoms]
                            label = '-'.join(dihedral_atoms_shifted).replace(' ', '_')
                            ax.plot(times, angles_deg[:, dihedral_idx], label=label, color=color, linestyle="",marker="o")
                            if idx_in_group == 1:  # Middle subplot gets a label
                                ax.set_ylabel('Dihedral Angle (°)', fontsize=14)
                            ax.legend(loc='best', fontsize=10)

                        axes[-1].set_xlabel('Time (ns)', fontsize=14)
                        plt.tight_layout()
                        # Save the plot
                        plot_filename = f'dihedrals_{group_idx+1}_{group["name"].replace(" ", "_")}.pdf'
                        plt.savefig(plot_filename)
                        plt.close()

                    else:
                        # For 'Ring Dihedrals' or any group not "Dihedrals for C*"
                        plt.figure(figsize=(7, 4))
                        colors_carbon = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, 6))
                        for idx_in_group, dihedral_idx in enumerate(dihedral_idxs):
                            plt.plot(times, 
                                    angles_deg[:, dihedral_idx], 
                                    label=f'τ {idx_in_group+1}', 
                                    color=colors_carbon[idx_in_group % len(colors_carbon)], 
                                    alpha=1 - idx_in_group / 10)
                        plt.xlabel('Time (ns)', fontsize=14)
                        plt.xlim((min(times), max(times)))
                        plt.ylabel('Dihedral Angle (°)', fontsize=14)
                        plt.title(f'Dihedral Angles Over Time - {group["name"]}', fontsize=16)
                        plt.legend(loc='best', fontsize=10)
                        plt.tight_layout()
                        # Save the plot
                        plot_filename = f'dihedrals_{group_idx+1}_{group["name"].replace(" ", "_")}.pdf'
                        plt.savefig(plot_filename)
                        plt.close()

                ############################################################################
                #  PART 2: SEPARATED Dihedral PLOTS (Histogram Approach with vertical offset)
                #
                #   - All "carbon" dihedrals in ONE figure
                #   - "oxygen" dihedrals in triple-wise chunks, each chunk in its own figure
                #
                ############################################################################

                import numpy as np

                # Define a helper function to split lists into chunks
                def chunks(lst, n):
                    """Yield successive n-sized chunks from lst."""
                    for i in range(0, len(lst), n):
                        yield lst[i : i + n]

                # Let's separate the dihedrals into carbon vs. oxygen.
                carbon_data = []
                oxygen_data = []

                for group_idx, group in enumerate(dihedral_groups):
                    dihedral_idxs = group_dihedral_indices[group_idx]
                    if group['name'].startswith('Ring'):
                        dihedral_type = "carbon"
                    else:
                        dihedral_type = "oxygen"

                    for dihedral_idx in dihedral_idxs:
                        angles_1d = angles_deg[:, dihedral_idx]  # in degrees
                        dihedral_atoms = dihedral_definitions[dihedral_idx]
                        dihedral_atoms_shifted = [shift_atom_number(atom) for atom in dihedral_atoms]
                        label = '-'.join(dihedral_atoms_shifted).replace(' ', '_')

                        if dihedral_type == "carbon":
                            carbon_data.append((dihedral_idx, angles_1d, label))
                        else:
                            oxygen_data.append((dihedral_idx, angles_1d, label))

                ############################################################################
                # Plot CARBON in one figure, offset each histogram, reversed order, dynamic linewidth
                ############################################################################
                if len(carbon_data) > 0:
                    fig_carbon, ax_carbon = plt.subplots(figsize=(5, 4), constrained_layout=True)
                    colors_carbon = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(carbon_data)))

                    # Reverse the order of the data so the largest offset is plotted first (in the background)
                    base_linewidth = 2
                    linewidth_decrement = 0.25  # Adjust as needed

                    for i_idx, (d_idx, angles_deg_array, label) in reversed(list(enumerate(carbon_data))):
                        # 1) Compute the histogram in memory
                        counts, bin_edges = np.histogram(
                            angles_deg_array,
                            bins=360,
                            range=(-180, 180),
                            density=True
                        )
                        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

                        # 2) Offset the counts, e.g., 0.02 per dihedral
                        offset = 0.02 * i_idx
                        offset_counts = counts + offset

                        # 3) Adjust linewidth: thicker for foreground (higher i_idx)
                        linewidth = base_linewidth - (linewidth_decrement * i_idx)
                        linewidth = max(linewidth, 1.0)  # Ensure linewidth doesn't get too small

                        # 4) Plot as a step line
                        ax_carbon.plot(
                            bin_centers,
                            offset_counts,
                            drawstyle='steps-mid',
                            label=label,
                            color=colors_carbon[i_idx],
                            linewidth=linewidth
                        )

                    # We don't want the y-axis to pretend it's real "density" anymore
                    ax_carbon.set_yticks([])
                    ax_carbon.set_ylabel("Density", fontsize=14)

                    ax_carbon.set_xlabel("Angle (°)", fontsize=14)
                    ax_carbon.set_xlim(-180, 180)
                    ax_carbon.set_title("Carbon (Ring) Dihedral Hist with Vertical Offset", fontsize=16)
                    ax_carbon.legend(fontsize=10)
                    ax_carbon.grid(False)

                    plt.savefig("dihedral_hist_carbon_offset.pdf")
                    plt.close()

                if len(carbon_data) > 0:
                    fig_carbon, ax_carbon = plt.subplots(figsize=(5, 4), constrained_layout=True)
                    colors_carbon = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(carbon_data)))
                    
                    # Define a grid for the smooth KDE curves
                    x_grid = np.linspace(-180, 180, 361)
                    
                    base_linewidth = 2
                    linewidth_decrement = 0.25  # Adjust as needed

                    # Reverse the order so that the highest offset (background) is plotted first.
                    for i_idx, (d_idx, angles_deg_array, label) in reversed(list(enumerate(carbon_data))):
                        # Compute the smooth KDE using gaussian_kde
                        kde = gaussian_kde(angles_deg_array)
                        density = kde(x_grid)
                        
                        # Apply an offset to separate overlapping curves
                        offset = 0.02 * i_idx
                        density_offset = density + offset
                        
                        # Adjust linewidth: thicker for curves with a higher index (foreground)
                        linewidth = base_linewidth - (linewidth_decrement * i_idx)
                        linewidth = max(linewidth, 1.0)  # Prevent the linewidth from becoming too thin
                        
                        # Plot the smooth KDE curve with the offset
                        ax_carbon.plot(
                            x_grid,
                            density_offset,
                            label=label,
                            color=colors_carbon[i_idx],
                            linewidth=1.5
                        )

                    ax_carbon.set_yticks([])
                    ax_carbon.set_ylabel("Density", fontsize=18)
                    ax_carbon.set_xlabel("Angle (°)", fontsize=18)
                    ax_carbon.set_xlim(-180, 180)
                    ax_carbon.set_title(f"({species_id})", fontsize=18)
                    ax_carbon.legend(fontsize=10)
                    ax_carbon.grid(False)

                    plt.savefig("dihedral_kde_carbon_offset.pdf")
                    plt.close()

                ############################################################################
                # Plot OXYGEN dihedrals in TRIPLE-WISE fashion, offset each histogram, reversed order, dynamic linewidth
                ############################################################################
                if len(oxygen_data) > 0:
                    chunk_size = 3
                    oxygen_chunks = list(chunks(oxygen_data, chunk_size))

                    for chunk_idx, chunk_data in enumerate(oxygen_chunks):
                        fig_oxy, ax_oxy = plt.subplots(figsize=(7, 4), constrained_layout=True)
                        colors_oxy = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(chunk_data)))

                        # Reverse the order within the chunk to plot background first
                        base_linewidth = 2
                        linewidth_decrement = 0.5  # Adjust as needed

                        for i_idx, (d_idx, angles_deg_array, label) in reversed(list(enumerate(chunk_data))):
                            # 1) Compute the histogram
                            counts, bin_edges = np.histogram(
                                angles_deg_array,
                                bins=360,
                                range=(-180, 180),
                                density=True
                            )
                            bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

                            # 2) Offset
                            offset = 0.01 * i_idx  # Adjust offset as needed
                            offset_counts = counts + offset

                            # 3) Adjust linewidth and color: thicker for foreground (higher i_idx)
                            linewidth = base_linewidth - (linewidth_decrement * i_idx)
                            linewidth = max(linewidth, 1.0)  # Ensure linewidth doesn't get too small
                            colors_list = ['darkblue', 'darkgreen', 'darkred']

                            # 4) Plot
                            ax_oxy.plot(
                                bin_centers,
                                offset_counts,
                                drawstyle='steps-mid',
                                label=label,
                                color=colors_list[i_idx],
                                linewidth=linewidth
                            )

                            
                            # 5) Add a horizontal dashed line at the offset level
                            ax_oxy.hlines(
                                offset,
                                xmin=-180,
                                xmax=180,
                                linestyles="--",
                                colors="black",
                                linewidth=1,
                                alpha=1  
                            )

                        ax_oxy.set_yticks([])
                        ax_oxy.set_ylabel("Density", fontsize=14)

                        ax_oxy.set_xlabel("Angle (°)", fontsize=14)
                        ax_oxy.set_xlim(-180, 180)
                        ax_oxy.set_title(f"Oxygen Dihedral Hist Offset - Group {chunk_idx+1}", fontsize=16)
                        ax_oxy.legend(fontsize=10)
                        ax_oxy.grid(False)

                        outname = f"dihedral_hist_oxygen_group{chunk_idx+1}_offset.pdf"
                        plt.savefig(outname)
                        plt.close()

                ############################################################################
                # FINAL SUMMARY PLOT: One histogram per oxygen chunk
                # We'll pick the first dihedral from each chunk_data (chunk_data[0])
                ############################################################################

                if len(oxygen_data) > 0:
                    # We'll reuse the same chunk_size and oxygen_chunks from above
                    # chunk_size = 3
                    # oxygen_chunks = list(chunks(oxygen_data, chunk_size))

                    fig_summary, ax_summary = plt.subplots(figsize=(7, 4), constrained_layout=True)

                    # Define a color map for however many chunks we have
                    colors_summary = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(oxygen_chunks)))

                    for chunk_idx, chunk_data in enumerate(oxygen_chunks):
                        # Just pick the first dihedral in this chunk.
                        # chunk_data could have 1–3 elements, so let's do a quick check:
                        if len(chunk_data) == 0:
                            continue
                        d_idx, angles_deg_array, label = chunk_data[0]

                        # Compute the histogram with no offset
                        counts, bin_edges = np.histogram(
                            angles_deg_array,
                            bins=360,
                            range=(-180, 180),
                            density=True
                        )
                        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
                        offset = 0.01 * chunk_idx
                        counts = counts + offset

                        # Plot with steps
                        ax_summary.plot(
                            bin_centers,
                            counts,
                            drawstyle='steps-mid',
                            label=f"{label}",
                            color=colors_summary[chunk_idx],
                            linewidth=1.5
                        )
                    ax_summary.set_yticks([])
                    ax_summary.set_xlabel("Angle (°)", fontsize=14)
                    ax_summary.set_ylabel("Density", fontsize=14)
                    ax_summary.set_xlim(-180, 180)
                    ax_summary.set_title("Summary: Oxygen Dihedrals", fontsize=16)
                    ax_summary.legend(fontsize=10)
                    ax_summary.grid(False)

                    outname = "dihedral_hist_oxygen_summary.pdf"
                    plt.savefig(outname)
                    plt.close()

                if len(oxygen_data) > 0:
                    fig_summary, ax_summary = plt.subplots(figsize=(7, 4), constrained_layout=True)

                    # Define a color map for however many chunks we have
                    colors_summary = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(oxygen_chunks)))

                    for chunk_idx, chunk_data in enumerate(oxygen_chunks):
                        # Just pick the first dihedral in this chunk.
                        if len(chunk_data) == 0:
                            continue
                        d_idx, angles_deg_array, label = chunk_data[0]

                        # Plot the KDE for the dihedral angle data.
                        # Using clip=(-180,180) to restrict the density estimation to the desired range.
                        sns.kdeplot(
                            x=angles_deg_array,
                            bw_adjust=1,  # Adjust the bandwidth if needed.
                            label=label,
                            color=colors_summary[chunk_idx],
                            linewidth=1.5,
                            ax=ax_summary,
                            clip=(-180, 180)
                        )

                    ax_summary.set_yticks([])
                    ax_summary.set_xlabel("Angle (°)", fontsize=14)
                    ax_summary.set_ylabel("Density", fontsize=14)
                    ax_summary.set_xlim(-180, 180)
                    ax_summary.set_title("Summary: Oxygen Dihedrals (KDE)", fontsize=16)
                    ax_summary.legend(fontsize=10)
                    ax_summary.grid(False)

                    outname = "dihedral_kde_oxygen_summary.pdf"
                    plt.savefig(outname)
                    plt.close()

                # Define a grid over the range of angles
                x_grid = np.linspace(-180, 180, 361)

                fig_summary, ax_summary = plt.subplots(figsize=(4,4), constrained_layout=True)

                # Define a color map for however many chunks we have
                colors_summary = plt.cm.get_cmap(cf.color_map)(np.linspace(0, 1, len(oxygen_chunks)))

                for chunk_idx, chunk_data in enumerate(oxygen_chunks):
                    if len(chunk_data) == 0:
                        continue
                    d_idx, angles_deg_array, label = chunk_data[0]

                    # Compute the smooth KDE using gaussian_kde
                    kde = gaussian_kde(dataset=angles_deg_array, bw_method="silverman")
                    density = kde(x_grid)
                    
                    # Compute an offset, e.g., 0.01 * chunk_idx
                    offset = 0.01 * chunk_idx
                    density_offset = density + offset

                    # Plot the smooth KDE curve with the offset (no 'steps-mid' drawstyle)
                    ax_summary.plot(x_grid, density_offset, label=label, 
                                    color=colors_summary[chunk_idx], linewidth=1.5)

                ax_summary.set_yticks([])
                ax_summary.set_xlabel("Angle (°)", fontsize=18)
                ax_summary.set_ylabel("Density", fontsize=18)
                ax_summary.set_xlim(-180, 180)
                ax_summary.set_title(f"({species_id})", fontsize=18)
                ax_summary.legend(fontsize=10)
                ax_summary.grid(False)

                outname = "dihedral_kde_offset_oxygen_summary.pdf"
                plt.savefig(outname)
                plt.close()