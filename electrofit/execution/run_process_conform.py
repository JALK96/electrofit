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

# trunk-ignore(ruff/E402)
from electrofit.main.process_conform import process_conform
from electrofit.helper.file_manipulation import find_file_with_extension, strip_extension
from electrofit.helper.config_parser import ConfigParser



def main_conform_processing():
    """
    Main function to process initial structure and create GROMACS input.
    """


    # Find current working directory
    home = os.getcwd()
    print("Working Directory:", home)

    # Get parent directory 
    extracted_conforms_dir = os.path.dirname(home)
    print("Extracted Conform Directory::", extracted_conforms_dir)

    # Go to the parent folder (the extracted_conforms directory) to open the configuration file (.ef) 
    os.chdir(extracted_conforms_dir)

    input = find_file_with_extension("ef")
    config = ConfigParser(input)

    # Go back to the home directory
    os.chdir(home)

    # Define file and molecule name
    pdb_file = find_file_with_extension("pdb")
    molecule_name = strip_extension(pdb_file)
    print("Processing conform:", molecule_name)

    # Define 
    base_scratch_dir = config.BaseScratchDir
    print("Scratch directory set to:", base_scratch_dir)
    residue_name = config.ResidueName
    print("Residue Name:", residue_name)
    net_charge = config.Charge
    print("Charge set to:", net_charge)
    adjust_sym = config.AdjustSymmetry
    print("AdjustSymmetry set to:", adjust_sym)
    protocol = config.Protocol
    print("Charge fit protocol set to:", protocol)
    ignore_sym = config.IgnoreSymmetry
    print("IgnoreSymmetry set to:", ignore_sym)



    # Process the initial structure
    process_conform(
        molecule_name=molecule_name,
        pdb_file=pdb_file,
        base_scratch_dir=base_scratch_dir,
        net_charge=net_charge,
        residue_name=residue_name,
        adjust_sym=adjust_sym,
        protocol=protocol,
        ignore_sym=ignore_sym
    )


if __name__ == "__main__":
    main_conform_processing()
