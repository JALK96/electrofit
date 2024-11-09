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
from electrofit.main.process_initial_structure import process_initial_structure
from electrofit.helper.file_manipulation import find_file_with_extension
from electrofit.helper.config_parser import ConfigParser

def main_processing(net_charge=None):
    """
    Main function to process initial structure and create GROMACS input.
    """

    # Define file paths and molecule name
    input = find_file_with_extension("ef")

    config=ConfigParser(input)

    base_scratch_dir = config.BaseScratchDir
    molecule_name = config.MoleculeName
    mol2_file = f"{molecule_name}.mol2"
    net_charge=config.Charge
    residue_name = config.ResidueName
    adjust_sym=config.AdjustSymmetry
    atom_type=config.AtomType
    additional_input = [] 

    # Delet later - just for debugging:
    additional_input.extend([find_file_with_extension("ef"), find_file_with_extension("gcrt.log"), find_file_with_extension("gcrt"), find_file_with_extension("gesp"), find_file_with_extension("chk")])
    
    if adjust_sym==True:
        additional_input.extend([find_file_with_extension("json")])
    

    # Process the initial structure
    process_initial_structure(
        molecule_name=molecule_name,
        mol2_file=mol2_file,
        base_scratch_dir=base_scratch_dir,
        additional_input=additional_input,
        residue_name=residue_name, 
        net_charge=net_charge,
        adjust_sym=adjust_sym,
        atom_type=atom_type,
        exit_screen=True # Set this to True if u execute via pis.sh, else set this to False or delete 
    )


if __name__ == "__main__":
    main_processing()
