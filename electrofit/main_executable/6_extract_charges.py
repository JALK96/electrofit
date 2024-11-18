import os
import sys
import shutil
import glob
import fnmatch
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

from electrofit.commands.run_commands import run_command, run_acpype
from electrofit.helper.file_manipulation import find_file_with_extension, parse_charges_from_mol2, extract_charges_from_subdirectories, adjust_atom_names, load_symmetry_groups
from electrofit.helper.config_parser import ConfigParser
from electrofit.helper.set_logging import setup_logging, reset_logging
from electrofit.helper.plotting import plot_charges_by_atom, plot_charges_by_symmetry

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
    print("AdjustSymmetry set to:", adjust_sym)

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




