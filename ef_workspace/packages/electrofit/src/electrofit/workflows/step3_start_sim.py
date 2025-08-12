import os

from electrofit.cli.run_commands import run_command

PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
project_path = PROJECT_PATH

dst = os.path.join(project_path, "process")


# Walk through the dst directory and its subdirectories
for root, dirs, files in os.walk(dst):
    for item in files:
        if item == "gmx.sh" and os.path.basename(root) == "run_gmx_simulation":
            # Construct the full path to the script
            script_path = os.path.join(root, item)
            # Run the script using bash
            run_command(f"bash {script_path}")
