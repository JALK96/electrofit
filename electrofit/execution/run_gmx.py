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
from electrofit.main.gmx_simulation import set_up_production
from electrofit.helper.file_manipulation import find_file_with_extension
from electrofit.helper.config_parser import ConfigParser




def main_production():
    """
    Main function to run molecular dynmaics simulation, i.e. performing setup (box definition, solvation, add ions and equilibration) and running the production run.
    """

    input = find_file_with_extension("ef")

    config = ConfigParser(input)

    base_scratch_dir = config.BaseScratchDir
    box_type = config.BoxType
    cation=config.Cation
    anion=config.Anion
    distance=config.BoxEdgeDistance
    conc=config.IonConcentration
    molecule_name=config.MoleculeName

    m_gro = find_file_with_extension("gro")
    MDP_dir = "MDP"


    # Check if the provided paths exist
    for path in [MDP_dir]:
        if not os.path.exists(path):
            print(f"The path {path} does not exist.")
            return

    set_up_production(m_gro=m_gro, MDP_dir=MDP_dir, base_scratch_dir=base_scratch_dir, molecule_name=molecule_name, box_type=box_type, cation=cation, anion=anion, d=distance, conc=conc)


if __name__ == "__main__":
    main_production()
