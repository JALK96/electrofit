import os
import sys
import logging

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

def main():
    # Define the base process directory
    process_dir = os.path.join(project_path, "process")

    # Loop through each subdirectory in the process directory
    for folder_name in os.listdir(process_dir):
        
        folder_path = os.path.join(process_dir, folder_name)

        if os.path.isdir(folder_path):
            # Define the 'run_final_gmx_simulation' directory within this folder
            run_final_sim_dir = os.path.join(folder_path, 'run_final_gmx_simulation')
            log_file = os.path.join(run_final_sim_dir, 'nojump_center.log')

            os.chdir(run_final_sim_dir)

            # Setup logging
            setup_logging(log_file)
            logging.info("Logging is set up.")

            run_command(
                'echo "0\n" | gmx trjconv -s md.tpr -f md.xtc -o md_nojump.xtc -pbc nojump', cwd=run_final_sim_dir
            )

            # Center molecule
            run_command(
                'echo "1\n0\n" | gmx trjconv -s md.tpr -f md_nojump.xtc -o md_center_new.xtc -center -pbc mol',
                cwd=run_final_sim_dir,
            )

if __name__ == "__main__":
    main()

