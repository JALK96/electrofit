import os
import logging
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# Set the style for seaborn
sns.set_context("talk")

def find_project_root(current_dir, project_name="electrofit"):
    """
    Walks up directories until it finds `project_name` or reaches root.
    Returns the outermost match of `project_name`.
    """
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

# Import your custom modules
from electrofit.commands.run_commands import run_command
from electrofit.helper.set_logging import setup_logging

def create_index_file(structure_file, index_file):
    """
    Creates an index file with groups P, P1, P2, P3, P4, P5 using gmx make_ndx.
    """
    logging.info("Starting index file creation with gmx make_ndx...")

    # Define the commands to be sent to gmx make_ndx
    make_ndx_commands = (
        "a P\n"
        "a P1\n"
        "a P2\n"
        "a P3\n"
        "a P4\n"
        "a P5\n"
        "q\n"
    )

    # Construct the shell command to pipe the commands into gmx make_ndx
    make_ndx_command = f'printf "{make_ndx_commands}" | gmx make_ndx -f {structure_file} -o {index_file}'
    run_command(make_ndx_command)
    logging.info("Index file created successfully.")

def run_pairdist_commands(trajectory, topology, index_file, groups, selection_group, output_prefix):
    """
    Runs gmx pairdist for each group in groups against the selection_group.
    Produces .xvg files with min distance vs time.
    """
    logging.info("Starting pair distance calculations with gmx pairdist...")

    for i, group in enumerate(groups, start=1):
        output_file = f"{output_prefix}{i}.xvg"
        pairdist_command = (
            f'gmx pairdist -f {trajectory} -s {topology} -n {index_file} '
            f'-ref "{group}" -sel "{selection_group}" -pbc no -rmpbc yes -o {output_file}'
        )
        logging.info(f"Running gmx pairdist for group '{group}' -> '{output_file}'...")
        run_command(pairdist_command)
        logging.info(f"Distance calculation for group '{group}' completed.")

    logging.info("All gmx pairdist commands executed successfully.")

def plot_all_distances_subplots(output_prefix, num_groups, plot_filename='all_distances_subplots.pdf'):
    """
    Plots all distance files (.xvg) as subplots (distance vs. time).
    """
    logging.info("Starting to plot all distance data as subplots...")
    # Temporarily suppress logging for cleaner output
    logging.disable(logging.CRITICAL)

    nrows, ncols = 2, 3  # 6 subplots in a 2x3 grid
    fig, axes = plt.subplots(nrows, ncols, figsize=(8, 6), sharex=True, sharey=True)
    axes = axes.flatten()

    for i in range(1, num_groups + 1):
        filename = f"{output_prefix}{i}.xvg"
        ax = axes[i-1]
        try:
            data = np.loadtxt(filename, comments=('#', '@'))
            time = data[:, 0] / 1000.0  # Convert ps to ns
            distance = data[:, 1]
            group_label = f"P{i}"

            ax.plot(time, distance, color='black', linestyle='-', linewidth=0.5)
            mean_distance = np.mean(distance)
            # Add a red dashed line for the mean
            ax.axhline(y=mean_distance, color='red', linestyle='--', linewidth=1.5)
            ax.text(0.95, 0.95, f"Mean: {mean_distance:.3f} nm",
                    ha='right', va='top', transform=ax.transAxes, color='red',
                    fontsize=12, bbox=dict(facecolor='white', alpha=0, edgecolor='none'))

            ax.set_title(f"{group_label} - Na$^+$", fontsize=14)
            ax.set_xlabel('Time (ns)')
            if i % ncols == 1:
                ax.set_ylabel('Distance (nm)')
            ax.grid(False)
        except Exception as e:
            ax.text(0.5, 0.5, 'Error loading data', ha='center',
                    va='center', transform=ax.transAxes, color='red')
            ax.set_title(f"Distance: P{i} - NA (Error)", fontsize=14)

    plt.tight_layout()
    #fig.suptitle('Minimum Distances Between Phosphorus Groups and Na-Ions', fontsize=18, y=1.05)
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    # plt.show()  # Uncomment if you want to display
    logging.disable(logging.NOTSET)
    logging.info("Subplot plotting completed successfully.")

# --------------------------------------------------------------------
# NEW FEATURE: Counting how many ions are within a certain cutoff
# --------------------------------------------------------------------

def run_ion_count_commands(trajectory, topology, index_file, groups, selection_group,
                           cutoff=0.5, output_prefix="ion_count_"):
    """
    Uses gmx select to count how many ions (selection_group) are within a given
    distance cutoff of each phosphate group. Produces .xvg files with counts over time.
    """
    logging.info(f"Starting ion count within {cutoff} nm for each group in {groups}...")

    for group in groups:
        # We'll name each output file based on the group, e.g. ion_count_P.xvg, ion_count_P1.xvg, etc.
        output_file = f"{output_prefix}{group}.xvg"
        # Construct the 'gmx select' expression
        select_expr = f'group "{selection_group}" and within {cutoff} of group "{group}"'
        command = (
            f'gmx select -f {trajectory} -s {topology} -n {index_file} '
            f'-select \'{select_expr}\' '
            f'-os {output_file} '        # .xvg with # of selected atoms vs time
            f'-on dummy_index.ndx '      # optional index file (not used further, but needed by gmx)
            f'2>&1'
        )
        logging.info(f"Running gmx select for group '{group}', output -> '{output_file}'...")
        run_command(command)
    
    logging.info("All ion count commands executed successfully.")

def plot_ion_counts_subplots(output_prefix, groups, plot_filename='ion_counts_subplots.pdf'):
    """
    Plots the number of ions within a cutoff for each group as subplots (count vs. time).
    """
    logging.info("Starting to plot ion counts as subplots...")
    logging.disable(logging.CRITICAL)

    num_groups = len(groups)
    nrows, ncols = 2, 3  # Adjust for 6 groups
    fig, axes = plt.subplots(nrows, ncols, figsize=(8, 6), sharex=True, sharey=True)
    axes = axes.flatten()

    for i, group in enumerate(groups):
        filename = f"{output_prefix}{group}.xvg"
        ax = axes[i]
        try:
            data = np.loadtxt(filename, comments=('#', '@'))
            time = data[:, 0] / 1000.0  # time in ns if .xvg is in ps
            ion_count = data[:, 1]
            mean_count = np.mean(ion_count)

            ax.plot(time, ion_count, color='darkblue', linestyle='-', linewidth=0.5)
            ax.axhline(y=mean_count, color='red', linestyle='--', linewidth=1.5)
            ax.text(0.95, 0.95, f"Mean: {mean_count:.1f}",
                    ha='right', va='top', transform=ax.transAxes, color='red',
                    fontsize=12, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

            ax.set_title(f"{group} Ion Count", fontsize=14)
            ax.set_xlabel('Time (ns)')
            if i % ncols == 0:
                ax.set_ylabel('# of Na+ ions')
            ax.grid(False)
        except Exception as e:
            ax.text(0.5, 0.5, 'Error loading data', ha='center',
                    va='center', transform=ax.transAxes, color='red')
            ax.set_title(f"{group} (Error)", fontsize=14)

    # If we have fewer than nrows*ncols subplots, hide the empty ones
    for j in range(num_groups, nrows*ncols):
        fig.delaxes(axes[j])

    plt.tight_layout()
    #fig.suptitle('Number of Na+ Ions Within Cutoff Distance', fontsize=18, y=1.05)
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    # plt.show()
    logging.disable(logging.NOTSET)
    logging.info("Ion count subplot plotting completed successfully.")

# Optional: Count how many ions are near the *entire molecule* (e.g., 'MOL' group)
def run_ion_count_whole_molecule(trajectory, topology, index_file, whole_group="Other",
                                 selection_group="NA", cutoff=0.5, output_xvg="ion_count_MOL.xvg"):
    """
    If you have an index group for the entire molecule (named e.g. 'MOL'),
    this function will count how many 'selection_group' atoms are within
    `cutoff` nm of that entire molecule.
    """
    logging.info(f"Starting ion count for entire molecule group '{whole_group}'...")
    select_expr = f'group "{selection_group}" and within {cutoff} of group "{whole_group}"'
    command = (
        f'gmx select -f {trajectory} -s {topology} -n {index_file} '
        f'-select \'{select_expr}\' '
        f'-os {output_xvg} '    
        f'-on dummy_mol_index.ndx '
        f'2>&1'
    )
    run_command(command)
    logging.info("Ion count for entire molecule completed successfully.")

def plot_whole_molecule_ion_count(xvg_file, plot_filename='ion_count_whole_mol.pdf'):
    """
    Plots a single line: # of NA ions within cutoff for the entire molecule vs time.
    """
    logging.info("Plotting ion count for the entire molecule...")
    try:
        data = np.loadtxt(xvg_file, comments=('#','@'))
        time = data[:, 0] / 1000.0
        ion_count = data[:, 1]
        mean_count = np.mean(ion_count)

        plt.figure(figsize=(6,4))
        plt.plot(time, ion_count, color='darkblue', linestyle='-', linewidth=0.5)
        plt.axhline(y=mean_count, color='red', linestyle='--', linewidth=1.5)
        plt.text(0.95, 0.1, f"Mean: {mean_count:.1f}",
                 ha='right', va='top', transform=plt.gca().transAxes, color='red',
                 fontsize=12, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        #plt.title('Na+ Ion Count Near Entire Molecule')
        plt.xlabel('Time (ns)')
        plt.ylabel('# of Na+ ions')
        plt.tight_layout()
        plt.savefig(plot_filename, dpi=300)
        # plt.show()
        logging.info("Ion count for entire molecule plotted successfully.")
    except Exception as e:
        logging.error(f"Error plotting whole molecule ion count: {e}")

# ----------------------------
# Main Execution Flow
# ----------------------------
def main():
    # Define the base process directory
    process_dir = os.path.join(project_path, "process")

    # Loop through each subdirectory in the process directory
    for folder_name in os.listdir(process_dir):
        
        folder_path = os.path.join(process_dir, folder_name)
        
        if os.path.isdir(folder_path):
            # Define the 'run_final_gmx_simulation' directory within this folder
            run_final_sim_dir = os.path.join(folder_path, 'run_final_gmx_simulation')
            
            if os.path.isdir(run_final_sim_dir):
                # Define the destination directory 'analyze_final_sim'
                dest_dir = os.path.join(folder_path, 'analyze_final_sim')
                dest_dir = os.path.join(dest_dir, 'NaP_dist_count')

                os.makedirs(dest_dir, exist_ok=True)
                os.chdir(dest_dir)

                # Define paths (adjust these as necessary)
                structure_file = os.path.join(run_final_sim_dir, 'md.gro')
                trajectory_file = os.path.join(run_final_sim_dir, 'md_center.xtc')
                topology_file = os.path.join(run_final_sim_dir, "md.tpr")
                index_file = "NA_P_index.ndx"  # We'll create this file
                selection_group = "NA"
                # We'll produce distances like distances_NaP1.xvg -> distances_NaP6.xvg
                output_prefix = "distances_NaP"
                log_file = "distances_NaP_gmx.log"

                # List of phosphorus groups
                p_groups = ["P", "P1", "P2", "P3", "P4", "P5"]

                # Setup logging
                setup_logging(log_file)
                logging.info("Logging is set up.")

                # Check if input files exist
                input_files = [structure_file, trajectory_file, topology_file]
                for file in input_files:
                    if not os.path.isfile(file):
                        logging.error(f"Required input file '{file}' not found. Exiting.")
                        sys.exit(1)
                logging.info("All required input files are present.")

                # Step 1: Create the index file (includes P, P1, P2, P3, P4, P5)
                create_index_file(structure_file, index_file)

                # Step 2: Run gmx pairdist for each phosphorus group
                run_pairdist_commands(
                    trajectory=trajectory_file,
                    topology=topology_file,
                    index_file=index_file,
                    groups=p_groups,
                    selection_group=selection_group,
                    output_prefix=output_prefix
                )

                # Step 3: Plot the distance data
                plot_all_distances_subplots(
                    output_prefix=output_prefix,
                    num_groups=len(p_groups),
                    plot_filename='all_distances_subplots.pdf'
                )

                # NEW STEPS: Count how many ions are close to each phosphate group
                # Adjust cutoff as desired, e.g., 0.5 nm
                cutoff_distance = 0.5
                ion_count_prefix = "ion_count_"

                run_ion_count_commands(
                    trajectory=trajectory_file,
                    topology=topology_file,
                    index_file=index_file,
                    groups=p_groups,
                    selection_group=selection_group,
                    cutoff=cutoff_distance,
                    output_prefix=ion_count_prefix
                )

                # Plot the ion count data for each group
                plot_ion_counts_subplots(
                    output_prefix=ion_count_prefix,
                    groups=p_groups,
                    plot_filename='ion_counts_subplots.pdf'
                )

                # OPTIONAL: If you have an index group for the entire molecule
                # (e.g., "MOL") in your index file, you could do:
                run_ion_count_whole_molecule(
                    trajectory=trajectory_file,
                    topology=topology_file,
                    index_file=index_file,
                    whole_group="Other",  # must exist in your index
                    selection_group="NA",
                    cutoff=cutoff_distance,
                    output_xvg="ion_count_MOL.xvg"
                )
            
                plot_whole_molecule_ion_count(
                    xvg_file="ion_count_MOL.xvg",
                    plot_filename="ion_count_whole_mol.pdf"
                )

                logging.info("All tasks (distance + ion count) completed successfully.\n")

if __name__ == "__main__":
    main()