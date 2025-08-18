import logging
import os

from electrofit.cli.run_commands import run_command
from electrofit.infra.logging import setup_logging

PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
project_path = PROJECT_PATH


def main():
    # Define the base process directory
    process_dir = os.path.join(project_path, "process")

    # Loop through each subdirectory in the process directory
    for folder_name in os.listdir(process_dir):
        folder_path = os.path.join(process_dir, folder_name)

        if os.path.isdir(folder_path):
            # Define the 'run_final_gmx_simulation' directory within this folder
            run_final_sim_dir = os.path.join(folder_path, "run_final_gmx_simulation")
            log_file = os.path.join(run_final_sim_dir, "nojump_center.log")

            os.chdir(run_final_sim_dir)

            # Setup logging
            setup_logging(log_file)
            logging.info("Logging is set up.")

            run_command(
                'echo "0\n" | gmx trjconv -s md.tpr -f md.xtc -o md_nojump.xtc -pbc nojump',
                cwd=run_final_sim_dir,
            )

            # Center molecule
            run_command(
                'echo "1\n0\n" | gmx trjconv -s md.tpr -f md_nojump.xtc -o md_center_new.xtc -center -pbc mol',
                cwd=run_final_sim_dir,
            )


if __name__ == "__main__":
    main()
