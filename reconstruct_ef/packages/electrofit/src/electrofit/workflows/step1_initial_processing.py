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

from electrofit.commands.run_commands import run_command

dst = os.path.join(project_path, "process")

# Walk through the dst directory and its subdirectories
for root, dirs, files in os.walk(dst):
    for item in files:
        if item == "pis.sh":
            # Construct the full path to the script
            script_path = os.path.join(root, item)
            # Run the script using bash
            run_command(f'bash {script_path}')