import os
import sys
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

script_dir = os.path.dirname(os.path.abspath(__file__))
project_path = find_project_root(current_dir=script_dir)
sys.path.append(project_path)
process_dir = os.path.join(project_path, "process")

def calculate_symmetric_group_averages(charges_dict_file, equivalent_groups_file):
    """
    Calculate the mean of the average charges for symmetric atoms.

    Parameters:
    - charges_dict_file: Path to the JSON file containing charges data for atoms.
        Example:
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


# Iterate over each subdirectory in the process directory
for sub_dir in os.listdir(process_dir):
    sub_dir_path = os.path.join(process_dir, sub_dir)
    
    # **Check if sub_dir is a directory**
    if not os.path.isdir(sub_dir_path):
        print(f"Skipping '{sub_dir_path}' as it is not a directory.")
        continue  # Skip to the next item in the loop

    # Create path to the results directory
    results_dir = os.path.join(sub_dir_path, "results")

    # Specify the input to calculate the symmetric group average
    charges_dict_file = os.path.join(results_dir, "charges_dict.json")
    equivalent_groups_file = os.path.join(results_dir, "symmetry_groups.json")

    # Calculate updated charges with symmetric group averages
    updated_charges_dict = calculate_symmetric_group_averages(charges_dict_file, equivalent_groups_file)

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

