import os
import sys
import shutil
import glob
import fnmatch
import json
import pandas as pd


def calculate_symmetric_group_averages(charges_dict_file, equivalent_groups_file):
    """
    Calculate the mean of the average charges for symmetric atoms.

    Parameters:
    - charges_dict_file: Path to the JSON file containing charges data for atoms.
        Example JSON:
        {
            "H1": {"average_charge": 0.12, "charges": [...]},
            "H3": {"average_charge": 0.11, "charges": [...]},
            "H6": {"average_charge": 0.13, "charges": [...]},
            ...
        }
    - equivalent_groups_file: Path to the JSON file containing symmetric atom groups.
        Example JSON:
        {
            "H1": ["H3", "H6"],
            "C1": ["C2"]
        }

    Returns:
    - updated_charges_dict: Dictionary with updated average charges for symmetric groups.
        Example:
        {
            "H1": {"average_charge": 0.12, "charges": [...]},
            "H3": {"average_charge": 0.12, "charges": [...]},
            "H6": {"average_charge": 0.12, "charges": [...]},
            ...
        }
    """
    # Load equivalent groups from the JSON file
    with open(equivalent_groups_file, "r") as f:
        equivalent_groups = json.load(f)

    # Load charges from the JSON file
    with open(charges_dict_file, "r") as f:
        charges_dict = json.load(f)

    updated_charges_dict = charges_dict.copy()  # Create a copy to modify

    # Iterate over each group in the equivalent groups
    for representative, group in equivalent_groups.items():
        # Add the representative atom to the group
        full_group = [representative] + group

        # Collect the average charges for all atoms in the group
        group_average_charges = [
            charges_dict[atom]["average_charge"] for atom in full_group if atom in charges_dict
        ]

        # Calculate the mean of the average charges
        if group_average_charges:
            group_average = sum(group_average_charges) / len(group_average_charges)

            # Update the average charge for all atoms in the group
            for atom in full_group:
                if atom in updated_charges_dict:
                    updated_charges_dict[atom]["average_charge"] = group_average

    return updated_charges_dict


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
process_dir = os.path.join(project_path, "process")

# ---------------------- Remove outlier --------------------
remove_outlier = False
# ---------------------- Remove outlier --------------------

from electrofit.commands.run_commands import run_command, run_acpype
from electrofit.helper.file_manipulation import find_file_with_extension, parse_charges_from_mol2, extract_charges_from_subdirectories, adjust_atom_names, load_symmetry_groups
from electrofit.helper.config_parser import ConfigParser
from electrofit.helper.set_logging import setup_logging, reset_logging
from electrofit.helper.plotting import plot_charges_by_atom, plot_charges_by_symmetry

# Iterate over each subdirectory in the process directory
for sub_dir in os.listdir(process_dir):
    sub_dir_path = os.path.join(process_dir, sub_dir)
    
    # **Check if sub_dir is a directory**
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

    # **Check if any .acpype directories/files are found**
    if not ac_files:
        print(f"No '.acpype' directories/files found in '{pis_dir}'. Skipping '{sub_dir_path}'.")
        continue  # Skip to the next subdirectory

    print(f"Acpype files found in {pis_dir}: {ac_files}")

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
    calc_group_average=config.CalculateGroupAverage # ether true or false (bool)
    print("CalculateGroupAverage set to:", calc_group_average)

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
    initial_charges_dict = adjust_atom_names(parse_charges_from_mol2(mol2_source_file_path))
  

    # Extract charges from all subdirectories
    atoms_dict = extract_charges_from_subdirectories(base_dir, results_dir)

    # Save dictionary to a JSON file
    with open(os.path.join(results_dir, "charges_dict.json"), "w") as file:
        json.dump(atoms_dict, file, indent=4)  # 'indent=4' makes the file human-readable
    with open(os.path.join(results_dir, "initial_charges_dict.json"), "w") as file:
        json.dump(initial_charges_dict, file, indent=4)  

    if adjust_sym:
        os.chdir(pis_dir)
        equiv_group_path = os.path.join(pis_dir, find_file_with_extension("json"))
        equiv_group = load_symmetry_groups(equiv_group_path)

        # Plot the charges and average charges
        plot_charges_by_symmetry(atoms_dict, initial_charges_dict, results_dir, equiv_group)

        if calc_group_average and not remove_outlier:
            charges_dict_file = os.path.join(results_dir, "charges_dict.json")
            # Calculate updated charges with symmetric group averages
            updated_charges_dict = calculate_symmetric_group_averages(charges_dict_file, equiv_group_path)

            # Go to the results directory and save the updated dictionary to a JSON file
            os.chdir(results_dir)
            with open(os.path.join(results_dir, "group_average_charges_dict.json"), "w") as file:
                json.dump(updated_charges_dict, file, indent=4)  # 'indent=4' makes the file human-readable

            # Write the average charges to the output file
            output_file = os.path.join(results_dir, "group_average_charges.txt")
            try:
                with open(output_file, 'w') as f:
                    f.write("#Atom_Name\tAverage_Charge\n")
                    for atom_name, atom_data in updated_charges_dict.items():
                        f.write(f"{atom_name}\t{atom_data['average_charge']:.4f}\n")
                print(f"Average charges successfully written to {output_file}")
            except Exception as e:
                print(f"An error occurred while writing to {output_file}: {e}")

            atom_names = list(updated_charges_dict.keys())
            updated_avg_charges = [updated_charges_dict[atom]['average_charge'] for atom in atom_names]
            updated_avg_charges_path = os.path.join(results_dir, "group_average_charges.chg")
            with open(updated_avg_charges_path, "w") as output:
                for i in updated_avg_charges:
                    output.write(str(round(i, 4)) + '\n')
            
            setup_logging(log_file_path)
            update_mol2 = os.path.join(project_path, "electrofit/execution/update_mol2.py")
            run_command(f"python {update_mol2} {mol2_source_file_path} {updated_avg_charges_path} {updated_mol2_file}")
            run_acpype(updated_mol2_file, charge, results_dir, atom_type, charges="user")
            reset_logging()

# ---------------------- Remove outlier --------------------
        # ---------------------- Remove outlier --------------------
        if remove_outlier and calc_group_average:
            data = atoms_dict
            # Convert the data into a DataFrame
            charges_data = {}

            for key, value in data.items():
                charges_data[key] = value['charges']

            df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in charges_data.items()]))

            # Function to remove outliers based on IQR
            def remove_outliers_iqr(series):
                Q1 = series.quantile(0.25)
                Q3 = series.quantile(0.75)
                IQR = Q3 - Q1
                is_outlier = (series < (Q1 - 1.5 * IQR)) | (series > (Q3 + 1.5 * IQR))
                return series.mask(is_outlier)

            # Remove outliers for each column
            df_no_outliers = df.apply(remove_outliers_iqr)

            # Compute the total charge
            total_net_charge = sum(
                atom_data['average_charge'] for atom_data in data.values()
            )

            # Prepare the cleaned data dictionary with the original structure
            cleaned_data = {}

            for atom in df_no_outliers.columns:
                # Get the cleaned charges list, drop NaNs
                cleaned_charges = df_no_outliers[atom].dropna().tolist()
                # Compute the new average charge
                if cleaned_charges:
                    new_average_charge = sum(cleaned_charges) / len(cleaned_charges)
                else:
                    new_average_charge = 0.0  # Set to 0.0 if no charges remain
                # Store in cleaned_data with the same structure
                cleaned_data[atom] = {
                    "charges": cleaned_charges,
                    "average_charge": new_average_charge
                }

            # Calculate the new total net charge
            new_total_net_charge = sum(
                atom_data['average_charge'] for atom_data in cleaned_data.values()
            )
            print("New Total Net Charge: ", new_total_net_charge)

            # Calculate the deviation from the required net charge
            required_net_charge = int(total_net_charge)
            print("Required Net Charge: ", required_net_charge)

            charge_deviation = required_net_charge - new_total_net_charge

            # Adjust the average charges proportionally
            total_charge = new_total_net_charge

            if total_charge == 0:
                # Distribute the charge deviation equally
                num_atoms = len(cleaned_data)
                charge_adjustment = charge_deviation / num_atoms
                for atom_data in cleaned_data.values():
                    atom_data['average_charge'] += charge_adjustment
            else:
                # Adjust charges proportionally
                for atom_data in cleaned_data.values():
                    proportion = atom_data['average_charge'] / total_charge if total_charge != 0 else 0
                    adjustment = proportion * charge_deviation
                    atom_data['average_charge'] += adjustment

            # Verify the adjusted total net charge
            adjusted_total_net_charge = sum(
                atom_data['average_charge'] for atom_data in cleaned_data.values()
            )
            print(f"Adjusted Total Net Charge: {adjusted_total_net_charge}")

            # Save the cleaned data to JSON in results_dir
            cleaned_adjusted_charges_file = os.path.join(results_dir, 'cleaned_adjusted_charges.json')
            with open(cleaned_adjusted_charges_file, 'w') as f:
                json.dump(cleaned_data, f, indent=4)

            # Update the path to the charges_dict_file
            charges_dict_file = cleaned_adjusted_charges_file

            # Calculate updated charges with symmetric group averages
            updated_charges_dict = calculate_symmetric_group_averages(charges_dict_file, equiv_group_path)

            # Save the updated dictionary to a JSON file in results_dir
            cleaned_adjusted_group_average_charges_file = os.path.join(results_dir, "cleaned_adjusted_group_average_charges_dict.json")
            with open(cleaned_adjusted_group_average_charges_file, "w") as file:
                json.dump(updated_charges_dict, file, indent=4)

            # Write the average charges to the output file
            output_file = os.path.join(results_dir, "cleaned_adjusted_group_average_charges.txt")
            try:
                with open(output_file, 'w') as f:
                    f.write("#Atom_Name\tAverage_Charge\n")
                    for atom_name, atom_data in updated_charges_dict.items():
                        f.write(f"{atom_name}\t{atom_data['average_charge']:.4f}\n")
                print(f"Average charges successfully written to {output_file}")
            except Exception as e:
                print(f"An error occurred while writing to {output_file}: {e}")

            atom_names = list(updated_charges_dict.keys())
            updated_avg_charges = [updated_charges_dict[atom]['average_charge'] for atom in atom_names]
            updated_avg_charges_path = os.path.join(results_dir, "cleaned_adjusted_group_average_charges.chg")
            with open(updated_avg_charges_path, "w") as output:
                for i in updated_avg_charges:
                    output.write(str(round(i, 4)) + '\n')
            
            # Define the updated mol2 file path
            update_mol2_file = os.path.join(results_dir, f"averaged_{molecule_name}_cleaned.mol2")

            setup_logging(log_file_path)
            update_mol2 = os.path.join(project_path, "electrofit/execution/update_mol2.py")
            run_command(f"python {update_mol2} {mol2_source_file_path} {updated_avg_charges_path} {update_mol2_file}")
            run_acpype(update_mol2_file, charge, results_dir, atom_type, charges="user")
            reset_logging()
        # ---------------------- Remove outlier --------------------
# ---------------------- Remove outlier --------------------


    plot_charges_by_atom(atoms_dict, initial_charges_dict, results_dir)

    if calc_group_average == False:
        setup_logging(log_file_path)
        average_charges = os.path.join(results_dir, "average_charges.chg")
        update_mol2 = os.path.join(project_path, "electrofit/execution/update_mol2.py")
        run_command(f"python {update_mol2} {mol2_source_file_path} {average_charges} {updated_mol2_file}")
        run_acpype(updated_mol2_file, charge, results_dir, atom_type, charges="user")
        reset_logging()




