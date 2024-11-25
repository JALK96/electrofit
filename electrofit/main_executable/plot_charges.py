import os
import sys
import shutil
import glob
import fnmatch
import json
import pandas as pd
import matplotlib.pyplot as plt

def plot_histograms(df, title, filename, adjusted_average_charges=None, color='darkred'):
    """
    Plot histograms of the DataFrame and save the plot.

    Parameters:
        df (pd.DataFrame): DataFrame containing the charges data.
        title (str): Title of the plot.
        filename (str): Filename to save the plot.
        adjusted_average_charges (dict, optional): Dictionary of adjusted average charges.
        color (str): Color of the histogram bars.
    """
    axes = df.hist(bins=20, figsize=(15, 10), color=color, alpha=0.9, grid=False)

    for ax, col in zip(axes.flatten(), df.columns):
        ax.set_title('')
        ax.set_xlabel(col)
        ax.set_ylabel('')

        if adjusted_average_charges is not None:
            # Get the adjusted average charge from the dictionary
            mean_value = adjusted_average_charges.get(col, None)
        else:
            # Calculate the mean of the data
            mean_value = df[col].mean()

        if mean_value is not None:
            # Plot a vertical dashed line at the mean
            ax.axvline(mean_value, color='black', linestyle='dashed', linewidth=1)

            # Add a text label with the mean value
            ax.text(
                0.95, 0.95,
                f'{mean_value:.2f}',
                color='black',
                fontsize=10,
                ha='right', va='top',
                transform=ax.transAxes
            )

    plt.suptitle(title, fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(filename)
    plt.close()

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

# Set up the project paths
script_dir = os.path.dirname(os.path.abspath(__file__))
project_path = find_project_root(current_dir=script_dir)
sys.path.append(project_path)
process_dir = os.path.join(project_path, "process")

# Import necessary modules from the project
from electrofit.helper.file_manipulation import find_file_with_extension, parse_charges_from_mol2, extract_charges_from_subdirectories, adjust_atom_names, load_symmetry_groups
from electrofit.helper.config_parser import ConfigParser
from electrofit.helper.plotting import plot_charges_by_atom, plot_charges_by_symmetry

# Iterate over each subdirectory in the process directory
for sub_dir in os.listdir(process_dir):
    sub_dir_path = os.path.join(process_dir, sub_dir)
    
    # Check if sub_dir is a directory
    if not os.path.isdir(sub_dir_path):
        print(f"Skipping '{sub_dir_path}' as it is not a directory.")
        continue  # Skip to the next item in the loop

    # Base directory containing the subdirectories with mol2 files
    base_dir = os.path.join(process_dir, sub_dir, "extracted_conforms")

    # Initial processing directory 
    pis_dir = os.path.join(process_dir, sub_dir, "run_gau_create_gmx_in")

    # Acpype directory
    abstract_ac_dir = os.path.join(pis_dir, "*.acpype")

    # Use glob to find all matching files or directories
    ac_files = glob.glob(abstract_ac_dir)

    # Check if any .acpype directories/files are found
    if not ac_files:
        print(f"No '.acpype' directories/files found in '{pis_dir}'. Skipping '{sub_dir_path}'.")
        continue  # Skip to the next subdirectory

    print(f"Acpype files found in {pis_dir}: {ac_files}")

    ac_dir = ac_files[0]

    results_dir = os.path.join(process_dir, sub_dir, "results")
    os.makedirs(results_dir, exist_ok=True)

    # Create an additional subdirectory for plots
    plots_dir = os.path.join(results_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    # Change to base directory containing configuration file (.ef)
    os.chdir(base_dir)
    config_file_path = os.path.join(base_dir, find_file_with_extension("ef"))
    # Optionally copy it into the results directory
    # shutil.copy2(config_file_path, results_dir)

    # Define molecule parameters
    config = ConfigParser(config_file_path)
    molecule_name = config.MoleculeName
    print(f"Processing: {molecule_name}")
    charge = config.Charge
    print(f"Charge set to: {charge}")
    atom_type = config.AtomType
    print(f"AtomType set to: {atom_type}")
    adjust_sym = config.AdjustSymmetry
    print("AdjustSymmetry set to:", adjust_sym)
    ignore_sym = config.IgnoreSymmetry
    print("IgnoreSymmetry set to:", ignore_sym)
    calc_group_average = config.CalculateGroupAverage  # either True or False (bool)
    print("CalculateGroupAverage set to:", calc_group_average)

    mol2_file_pattern = f"*{atom_type}.mol2"

    for file_name in os.listdir(ac_dir):
        if fnmatch.fnmatch(file_name, mol2_file_pattern):
            mol2_source_file_path = os.path.join(ac_dir, file_name)
            mol2_file_name = file_name
            break  # Stop after finding the first matching file

    # Extract the initial charges
    initial_charges_dict = adjust_atom_names(parse_charges_from_mol2(mol2_source_file_path))

    # Extract charges from all subdirectories
    atoms_dict = extract_charges_from_subdirectories(base_dir, results_dir)

    # Remove code that saves dictionaries
    # with open(os.path.join(results_dir, "charges_dict.json"), "w") as file:
    #     json.dump(atoms_dict, file, indent=4)
    # with open(os.path.join(results_dir, "initial_charges_dict.json"), "w") as file:
    #     json.dump(initial_charges_dict, file, indent=4)

    # Plotting charges
    if adjust_sym:
        os.chdir(pis_dir)
        equiv_group_path = os.path.join(pis_dir, find_file_with_extension("json"))
        equiv_group = load_symmetry_groups(equiv_group_path)

        # Plot the charges and average charges
        plot_charges_by_symmetry(atoms_dict, initial_charges_dict, plots_dir, equiv_group)
    else:
        # Plot the charges by atom
        plot_charges_by_atom(atoms_dict, initial_charges_dict, plots_dir)

    # Convert atoms_dict to DataFrame for plotting histograms
    charges_data = {}
    for key, value in atoms_dict.items():
        charges_data[key] = value['charges']
    df = pd.DataFrame(charges_data)

    # Plot histograms and save to plots_dir
    plot_histograms(
        df,
        title='Charge Distribution',
        filename=os.path.join(plots_dir, 'histogram.pdf'),
        color='darkred'
    )