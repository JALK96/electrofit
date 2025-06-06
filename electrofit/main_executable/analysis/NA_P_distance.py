import os
import logging
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
# Set the style for seaborn
sns.set_context("talk")

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

from electrofit.commands.run_commands import run_command
from electrofit.helper.set_logging import setup_logging


def create_index_file(structure_file, index_file):
    """
    Creates an index file with groups P, P1, P2, P3, P4, P5 using gmx make_ndx.

    Parameters:
        structure_file (str): Path to the structure file (.gro).
        index_file (str): Path to the output index file (.ndx).
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
    # Using printf to handle the newline characters correctly across different shells
    make_ndx_command = f'printf "{make_ndx_commands}" | gmx make_ndx -f {structure_file} -o {index_file}'

    # Execute the make_ndx command
    run_command(make_ndx_command)
    logging.info("Index file created successfully.")

def run_pairdist_commands(trajectory, topology, index_file, groups, selection_group, output_prefix):
    """
    Runs gmx pairdist for each group in groups against the selection_group.

    Parameters:
        trajectory (str): Path to the trajectory file (.xtc).
        topology (str): Path to the topology file (.tpr).
        index_file (str): Path to the index file (.ndx).
        groups (list): List of reference groups (e.g., ['P', 'P1', ..., 'P5']).
        selection_group (str): The selection group to calculate distances against (e.g., 'NA').
        output_prefix (str): Prefix for the output files.
    """
    logging.info("Starting pair distance calculations with gmx pairdist...")

    for i, group in enumerate(groups, start=1):
        output_file = f"{output_prefix}{i}.xvg"
        pairdist_command = (
            f'gmx pairdist -f {trajectory} -s {topology} -n {index_file} '
            f'-ref "{group}" -sel "{selection_group}" -pbc no -rmpbc yes -o {output_file}'
        )
        logging.info(f"Running gmx pairdist for group '{group}' and outputting to '{output_file}'...")
        run_command(pairdist_command)
        logging.info(f"Distance calculation for group '{group}' completed.")

    logging.info("All gmx pairdist commands executed successfully.")

def plot_all_distances_subplots(output_prefix, num_groups, plot_filename='all_distances_subplots.png'):
    """
    Plots all distance files generated by gmx pairdist as individual subplots.

    Parameters:
        output_prefix (str): Prefix of the distance files (e.g., 'distances_P').
        num_groups (int): Number of groups to plot.
        plot_filename (str): Filename to save the plot.
    """
    logging.info("Starting to plot all distance data as subplots...")
    # Suppress logging
    logging.disable(logging.CRITICAL)

    # Define the number of rows and columns for subplots
    # For 6 plots, a 2x3 grid is suitable
    nrows, ncols = 2, 3

    # Create a figure with specified size
    fig, axes = plt.subplots(nrows, ncols, figsize=(8, 6), sharex=True, sharey=True)
    axes = axes.flatten()  # Flatten to easily iterate over

    for i in range(1, num_groups + 1):
        filename = f"{output_prefix}{i}.xvg"
        ax = axes[i-1]  # Select the subplot
        try:
            # Load data, skipping comment lines
            data = np.loadtxt(filename, comments=('#', '@'))
            time = data[:, 0] / 1000
            distance = data[:, 1]
            # Determine group label
            group_label = f"P{i}" 
            # Plot data on the subplot
            ax.plot(time, distance, color='black', linestyle='-', linewidth=0.5, marker=None)
            # Calculate mean distance
            mean_distance = np.mean(distance)

            # Plot horizontal red line at mean distance
            ax.axhline(y=mean_distance, color='red', linestyle='--', linewidth=1.5, label='Mean Distance')

            # Add red-colored text label for mean distance
            # Position the text near the top of the subplot
            ax.text(0.95, 0.95, f"Mean: {mean_distance:.3f} nm",
                    horizontalalignment='right', verticalalignment='top',
                    transform=ax.transAxes, color='red', fontsize=12,
                    bbox=dict(facecolor='white', alpha=0, edgecolor='none'))

            ax.set_title(f"{group_label} - Na$^+$", fontsize=14)
            ax.set_xlabel('Time (ns)')
            if i % ncols == 1:  # Set y-label only on the first column
                ax.set_ylabel('Distance (nm)')
            ax.grid(False)
        except Exception as e:
            ax.text(0.5, 0.5, 'Error loading data', horizontalalignment='center', 
                    verticalalignment='center', transform=ax.transAxes, color='red')
            ax.set_title(f"Distance: {group_label} - NA (Error)", fontsize=14)

    # Adjust layout to prevent overlap
    plt.tight_layout()

    # Add a main title for all subplots
    fig.suptitle('Minimum Distances Between Phosphorus Groups and Na-Ions', fontsize=18, y=1.05)

    # Save the figure
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    logging.info(f"Subplot figure saved as '{plot_filename}'.")

    # Show the plot
    #plt.show()
    # Re-enable logging if needed (optional)
    logging.disable(logging.NOTSET)

    logging.info("Subplot plotting completed successfully.")

# ----------------------------
# Main Execution Flow
# ----------------------------

def main():
    # Define the base process directory
    process_dir = os.path.join(project_path, "process")

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

                # Define paths (adjust these as necessary)
                structure_file = os.path.join(run_final_sim_dir, 'md.gro')
                
                trajectory_file = os.path.join(run_final_sim_dir, 'md_center.xtc')
                index_file = "NA_P_index.ndx"
                topology_file = os.path.join(run_final_sim_dir, "md.tpr")
                selection_group = "NA"
                output_prefix = "distances_NaP"  # Generates distances_P1.xvg to distances_P6.xvg
                log_file = "distances_NaP_gmx.log"

                # List of phosphorus groups
                p_groups = ["P", "P1", "P2", "P3", "P4", "P5"]

                # Set up logging
                setup_logging(log_file)
                logging.info("Logging is set up.")

                # Check if input files exist
                input_files = [structure_file, trajectory_file, topology_file]
                for file in input_files:
                    if not os.path.isfile(file):
                        logging.error(f"Required input file '{file}' not found. Exiting.")
                        sys.exit(1)
                logging.info("All required input files are present.")

                # Step 1: Create the index file
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

                # Optional Step 3: Plot the distance data
                plot_all_distances_subplots(output_prefix=output_prefix, num_groups=len(p_groups))

                logging.info("All tasks completed successfully.")

if __name__ == "__main__":
    main()