import os
import sys
import glob
import fnmatch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def create_atom_color_mapping(atom_names, symmetry_groups):
    # List of colors to use for the groups
    group_colors_list = ['darkred', 'darkgreen', 'darkorange', 'purple', 'royalblue',
                         'lightcoral', 'deepskyblue', 'mediumvioletred', 'orange',
                         'olive', 'teal', 'dodgerblue', 'darkkhaki', 'salmon',
                         'firebrick', 'olivedrab', 'palevioletred']
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

def plot_histograms(df_original, df_adjusted, title, filename, adjusted_average_charges=None, symmetric_atoms=None, atom_to_color=None, symmetric_groups=None, combine_original_data=False, remove_outlier=False):
    """
    Plot histograms of the DataFrame and save the plot.

    Parameters:
        df_original (pd.DataFrame): Original charges data.
        df_adjusted (pd.DataFrame): Adjusted charges data.
        title (str): Title of the plot.
        filename (str): Filename to save the plot.
        adjusted_average_charges (dict, optional): Dictionary of adjusted average charges.
        symmetric_atoms (set, optional): Set of atom names that are symmetric.
        atom_to_color (dict, optional): Mapping of atom names to colors.
        symmetric_groups (dict, optional): Dictionary of symmetry groups.
    """
    num_atoms = len(df_original.columns)
    num_cols = min(4, num_atoms)
    num_rows = (num_atoms + num_cols - 1) // num_cols  # Calculate rows needed

    fig, axes = plt.subplots(num_rows, num_cols, figsize=(4*num_cols, 3*num_rows))
    axes = axes.flatten()

    for i, col in enumerate(df_original.columns):
        ax = axes[i]
        data_original = df_original[col].dropna()  # Original data
        data_adjusted = df_adjusted[col].dropna()  # Adjusted data

        # Determine if the atom is symmetric
        is_symmetric = symmetric_atoms is not None and col in symmetric_atoms
        # Get the color for the atom
        if atom_to_color is not None:
            color = atom_to_color.get(col, 'darkblue')
        else:
            color = 'darkblue'

        # Calculate common bins for both datasets
        if is_symmetric and symmetric_groups is not None:
            # For symmetric atoms, combine data from the symmetry group
            group_found = False
            for rep, group in symmetric_groups.items():
                group_atoms = [rep] + group
                if col in group_atoms:
                    # Combine charges for the group
                    combined_data_adjusted = pd.Series(dtype=float)
                    combined_data_original = pd.Series(dtype=float)

                    for atom in group_atoms:
                        combined_data_adjusted = pd.concat(
                            [combined_data_adjusted, df_adjusted[atom].dropna()],
                            ignore_index=True
                        )
                        combined_data_original = pd.concat(
                            [combined_data_original, df_original[atom].dropna()],
                            ignore_index=True
                        )

                    if combine_original_data == False:
                        combined_data_original = data_original
                         # Use combined data to calculate bins
                        data_combined = pd.concat([combined_data_original, combined_data_adjusted], ignore_index=True)
                        bins = np.histogram_bin_edges(data_combined, bins=20)
                        # Plot original histogram in background
                        ax.hist(combined_data_original, bins=bins, color=color, alpha=0.9, edgecolor='black', label='Original')
                        # Plot combined adjusted histogram
                        ax.hist(combined_data_adjusted, bins=bins, color=color, alpha=0.5, edgecolor='black', label='Combined')
                        group_found = True
                        break

                    else:
                        # Use combined data to calculate bins
                        data_combined = pd.concat([combined_data_original, combined_data_adjusted], ignore_index=True)
                        bins = np.histogram_bin_edges(data_combined, bins=20)
                        # Plot original histogram in background
                        ax.hist(combined_data_original, bins=bins, color=color, alpha=0.5, edgecolor='black', label='Comb. Original')
                        # Plot combined adjusted histogram
                        ax.hist(combined_data_adjusted, bins=bins, color=color, alpha=0.9, edgecolor='black', label='Comb. Clipped')
                        group_found = True
                        break
            if not group_found:
                # Atom is symmetric but not found in groups (shouldn't happen)
                # Use data from the individual atom
                data_combined = pd.concat([data_original, data_adjusted], ignore_index=True)
                bins = np.histogram_bin_edges(data_combined, bins=20)
                # Plot original and adjusted histograms
                ax.hist(data_original, bins=bins, color=color, alpha=0.5, edgecolor='black', label='Original')
                if remove_outlier:
                    ax.hist(data_adjusted, bins=bins, color=color, alpha=0.9, edgecolor='black', label='Clipped')
        else:
            # Non-symmetric atom
            # Use data from the individual atom
            data_combined = pd.concat([data_original, data_adjusted], ignore_index=True)
            bins = np.histogram_bin_edges(data_combined, bins=20)
            # Plot original and adjusted histograms
            ax.hist(data_original, bins=bins, color=color, alpha=0.5, edgecolor='black', label='Original')
            if remove_outlier:
                ax.hist(data_adjusted, bins=bins, color=color, alpha=0.9, edgecolor='black', label='Clipped')

        ax.set_title(col)
        ax.set_xlabel('Charge')
        ax.set_ylabel('Frequency')

        if adjusted_average_charges is not None:
            # For symmetric atoms, use the group average charge
            if is_symmetric:
                mean_value = adjusted_average_charges.get(col, None)
            else:
                # For non-symmetric atoms, use their individual mean
                mean_value = data_adjusted.mean()
        else:
            # Calculate the mean of the adjusted data
            mean_value = data_adjusted.mean()

        if mean_value is not None:
            # Set line color based on symmetry
            line_color = 'red' if is_symmetric else 'black'
            # Plot a vertical dashed line at the mean
            ax.axvline(mean_value, color=line_color, linestyle='dashed', linewidth=2, label=f'{mean_value:.2f}')
            # Add a text label with the mean value
            #ax.text(
            #    0.95, 0.95,
            #    f'{mean_value:.2f}',
            #    color=line_color,
            #    fontsize=10,
            #    ha='right', va='top',
            #    transform=ax.transAxes
            #)

        # Add legend
        ax.legend()

    # Remove any unused subplots
    for j in range(i+1, len(axes)):
        fig.delaxes(axes[j])

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

def calculate_symmetric_group_averages(charges_dict, equivalent_groups):
    """
    Calculate the mean of the average charges for symmetric atoms.

    Parameters:
    - charges_dict: Dictionary containing charges data for atoms.
    - equivalent_groups: Dictionary containing symmetric atom groups.

    Returns:
    - updated_charges_dict: Dictionary with updated average charges for symmetric groups.
    """
    updated_charges_dict = charges_dict.copy()  # Create a copy to modify

    # Iterate over each group in the equivalent groups
    for representative, group in equivalent_groups.items():
        # Add the representative atom to the group
        full_group = [representative] + group

        # Collect the charges for all atoms in the group
        group_charges = []
        for atom in full_group:
            if atom in charges_dict:
                group_charges.extend(charges_dict[atom]['charges'])

        # Calculate the mean of the charges
        if group_charges:
            group_average = sum(group_charges) / len(group_charges)

            # Update the average charge for all atoms in the group
            for atom in full_group:
                if atom in updated_charges_dict:
                    updated_charges_dict[atom]['average_charge'] = group_average

    return updated_charges_dict

def combine_and_calculate_symmetric_group_averages(charges_dict, equivalent_groups):
    """
    Calculate the mean of the average charges for symmetric atoms, update the average charges,
    and create a combined dataset with expanded charges lists for symmetric atoms.

    Parameters:
    - charges_dict (dict): Dictionary containing charges data for atoms.
      Expected format:
          {
              'C1': {'charges': [...], 'average_charge': ...},
              'C2': {'charges': [...], 'average_charge': ...},
              ...
          }

    - equivalent_groups (dict): Dictionary containing symmetric atom groups.
      Expected format:
          {
              'C1': ['C2', 'C3'],
              'O1': ['O2', 'O3'],
              ...
          }

    Returns:
    - updated_charges_dict (dict): Dictionary with updated average charges for symmetric groups.
    - combined_charges_dict (dict): Dictionary with updated average charges and expanded charges lists of combined datasets for symmetric atom entries.
    """
    import copy
    # Deep copy to avoid modifying the original dictionaries
    updated_charges_dict = copy.deepcopy(charges_dict)
    combined_charges_dict = copy.deepcopy(charges_dict)

    # Iterate over each group in equivalent_groups
    for representative, group in equivalent_groups.items():
        # Form the complete group including the representative
        full_group = [representative] + group

        # Collect charges from each atom in the group
        group_charges = []
        for atom in full_group:
            if atom in charges_dict:
                group_charges.extend(charges_dict[atom]['charges'])
            else:
                print(f"Warning: Atom '{atom}' not found in charges_dict.")

        # Calculate the mean charge for the group
        if group_charges:
            group_average = sum(group_charges) / len(group_charges)

            # Update the average_charge for each atom in updated_charges_dict
            for atom in full_group:
                if atom in updated_charges_dict:
                    updated_charges_dict[atom]['average_charge'] = group_average
                else:
                    print(f"Warning: Atom '{atom}' not found in updated_charges_dict.")

            # Update the average_charge and charges list for each atom in combined_charges_dict
            for atom in full_group:
                if atom in combined_charges_dict:
                    combined_charges_dict[atom]['average_charge'] = group_average
                    combined_charges_dict[atom]['charges'] = group_charges.copy()
                else:
                    print(f"Warning: Atom '{atom}' not found in combined_charges_dict.")
        else:
            print(f"Warning: No charges found for group '{representative}'.")

    # Atoms not in symmetric groups remain unchanged in both dictionaries
    return updated_charges_dict, combined_charges_dict

# Set up the project paths
script_dir = os.path.dirname(os.path.abspath(__file__))
project_path = find_project_root(current_dir=script_dir)
sys.path.append(project_path)
process_dir = os.path.join(project_path, "process")

# Import necessary modules from the project
from electrofit.helper.file_manipulation import find_file_with_extension, parse_charges_from_mol2, extract_charges_from_subdirectories, adjust_atom_names, load_symmetry_groups
from electrofit.helper.config_parser import ConfigParser
from electrofit.helper.plotting import plot_charges_by_atom, plot_charges_by_symmetry, plot_charges_by_atom_new_4

# ------ Interact --------------
# Remove outlier functionality
remove_outlier = False  # Set to True to perform outlier removal based on 1.5 IQR rule
# ------ Interact --------------

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

    # Find the mol2 source file path
    mol2_source_file_path = None
    for file_name in os.listdir(ac_dir):
        if fnmatch.fnmatch(file_name, mol2_file_pattern):
            mol2_source_file_path = os.path.join(ac_dir, file_name)
            break  # Stop after finding the first matching file

    if mol2_source_file_path is None:
        print(f"No mol2 file matching pattern '{mol2_file_pattern}' found in '{ac_dir}'. Skipping.")
        continue

    # Extract the initial charges
    initial_charges_dict = adjust_atom_names(parse_charges_from_mol2(mol2_source_file_path))

    # Extract charges from all subdirectories
    atoms_dict = extract_charges_from_subdirectories(base_dir, results_dir)

    # Plotting charges
    if adjust_sym:
        os.chdir(pis_dir)
        equiv_group_path = os.path.join(pis_dir, find_file_with_extension("json"))
        equiv_group = load_symmetry_groups(equiv_group_path)

        # Plot the charges and average charges
        plot_charges_by_symmetry(atoms_dict, initial_charges_dict, plots_dir, equiv_group)
        plot_charges_by_atom_new_4(atoms_dict, initial_charges_dict, plots_dir)
    else:
        # Plot the charges by atom
        plot_charges_by_atom(atoms_dict, initial_charges_dict, plots_dir)

    # Convert atoms_dict to DataFrame for plotting histograms
    charges_data = {}
    for key, value in atoms_dict.items():
        charges_data[key] = value['charges']
    df = pd.DataFrame(charges_data)

    # Create atom-to-color mapping if adjust_sym is True
    if adjust_sym:
        # Load symmetry groups
        symmetric_groups = equiv_group  # Dictionary of symmetry groups
        # Create a set of symmetric atoms
        symmetric_atoms = set()
        for group in symmetric_groups.values():
            symmetric_atoms.update(group)
        symmetric_atoms.update(symmetric_groups.keys())
        # Create atom-to-color mapping
        atom_to_color = create_atom_color_mapping(df.columns.tolist(), symmetric_groups)
    else:
        symmetric_atoms = set()
        atom_to_color = {atom: 'darkblue' for atom in df.columns}

    # Plot histograms before removing outliers and save to plots_dir
    plot_histograms(
        df_original=df,
        df_adjusted=df,
        title='Charge Distribution',
        filename=os.path.join(plots_dir, 'hist.pdf'),
        atom_to_color=atom_to_color,
        symmetric_atoms=symmetric_atoms
    )



    # ---------------------- Remove outlier --------------------

    if remove_outlier:
        print("Remove Outliers ...")
        # Function to get outlier mask based on IQR
        def get_outlier_mask(series):
            Q1 = series.quantile(0.25)
            Q3 = series.quantile(0.75)
            IQR = Q3 - Q1
            is_outlier = (series < (Q1 - 1.5 * IQR)) | (series > (Q3 + 1.5 * IQR))
            return is_outlier

        # Initialize a boolean mask for all rows (conformers), default False
        outlier_mask = pd.Series(False, index=df.index)

        # Identify outliers for each atom and update the mask
        for column in df.columns:
            series = df[column]
            is_outlier = get_outlier_mask(series)
            outlier_mask = outlier_mask | is_outlier  # Combine masks using logical OR

        # Remove rows (conformers) where any charge is an outlier
        df_no_outliers = df[~outlier_mask].reset_index(drop=True)

        # Calculate new average charges after outlier removal
        adjusted_average_charges = df_no_outliers.mean().to_dict()

        # ------------------ New Section for Group Average Charges ------------------

        if calc_group_average and adjust_sym:
            # Calculate group average charges
            cleaned_atoms_dict = {}
            for atom in df_no_outliers.columns:
                cleaned_atoms_dict[atom] = {
                    'charges': df_no_outliers[atom].tolist(),
                    'average_charge': adjusted_average_charges[atom]
                }

            updated_charges_dict = calculate_symmetric_group_averages(cleaned_atoms_dict, symmetric_groups)

            # Prepare the adjusted average charges for plotting
            adjusted_average_charges = {atom: data['average_charge'] for atom, data in updated_charges_dict.items()}

            # Plot histograms with group average charges and overlay original distributions
            plot_histograms(
                df_original=df,
                df_adjusted=df_no_outliers,
                title='Charge Distribution with Group Average Charges',
                filename=os.path.join(plots_dir, 'hist_group_average_clipped_charges.pdf'),
                adjusted_average_charges=adjusted_average_charges,
                symmetric_atoms=symmetric_atoms,
                atom_to_color=atom_to_color,
                symmetric_groups=symmetric_groups,  # Pass the symmetric groups
                combine_original_data=True,
                remove_outlier=True
            )

            # Plot charges by symmetry using the updated charges
            plot_charges_by_symmetry(updated_charges_dict, initial_charges_dict, plots_dir, equiv_group)

        elif adjust_sym:
            # Plot histograms with adjusted average charges (without group averaging) and overlay original distributions
            plot_histograms(
                df_original=df,
                df_adjusted=df_no_outliers,
                title='Clipped Charge Distribution with Symmetry and Clipped Average Charges',
                filename=os.path.join(plots_dir, 'hist_average_clipped_charges.pdf'),
                adjusted_average_charges=adjusted_average_charges,
                atom_to_color=atom_to_color,
                symmetric_atoms=symmetric_atoms,
                remove_outlier=True
            )


        else:
            # Plot histograms with adjusted average charges (without group averaging) and overlay original distributions
            plot_histograms(
                df_original=df,
                df_adjusted=df_no_outliers,
                title='Clipped Charge Distribution with no Symmetry and Clipped Average Charges',
                filename=os.path.join(plots_dir, 'hist_average_clipped_charges.pdf'),
                adjusted_average_charges=adjusted_average_charges,
                atom_to_color=atom_to_color,
                remove_outlier=True
            )

            # Plot the charges by atom using the cleaned data
            cleaned_atoms_dict = {}
            for atom in df_no_outliers.columns:
                cleaned_atoms_dict[atom] = {
                    'charges': df_no_outliers[atom].tolist(),
                    'average_charge': adjusted_average_charges[atom]
                }
            plot_charges_by_atom(cleaned_atoms_dict, initial_charges_dict, plots_dir)
    # ---------------------- Remove outlier --------------------

    else:
        if calc_group_average and adjust_sym:
            # Compute average of symmetric atoms
            updated_charges_dict, combind_charges_dict = combine_and_calculate_symmetric_group_averages(atoms_dict, symmetric_groups)

            # Prepare the adjusted average charges for plotting
            group_average_charges = {atom: data['average_charge'] for atom, data in updated_charges_dict.items()}

            # Plot histograms with group average charges and overlay original distributions
            plot_histograms(
                df_original=df,
                df_adjusted=df,
                title='Charge Distribution with Group Average Charges',
                filename=os.path.join(plots_dir, 'hist_group_average_charges.pdf'),
                adjusted_average_charges=group_average_charges,
                symmetric_atoms=symmetric_atoms,
                atom_to_color=atom_to_color,
                symmetric_groups=symmetric_groups,  # Pass the symmetric groups
                combine_original_data=False
            )

            # Plot charges by symmetry using the updated charges
            # Saves the average charges in updated_charges_dict to average_charges.chg
            plot_charges_by_symmetry(updated_charges_dict, initial_charges_dict, plots_dir, equiv_group)

            # Note: This plot only works for IP6 Configs 
            plot_charges_by_atom_new_4(atoms_dict, initial_charges_dict, plots_dir, combind_charges_dict, equivalent_groups=symmetric_groups)
        # ---------------------------------------------------------------------------

