import os
import time 
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

# Step 1: Collect all pc.sh script paths
all_scripts = []
for root, dirs, files in os.walk(dst):
    for item in files:
        if item == "pc.sh":
            all_scripts.append(os.path.join(root, item))

# Step 2: Run scripts in batches of 20 every 15 minutes
batch_size = 20
executed_scripts = set()  # Track executed scripts

while len(executed_scripts) < len(all_scripts):
    # Select the next 20 scripts that haven't been executed yet
    scripts_to_run = [script for script in all_scripts if script not in executed_scripts][:batch_size]
    
    # Run each script in the batch
    for script_path in scripts_to_run:
        run_command(f'bash {script_path}')
        executed_scripts.add(script_path)  # Mark this script as executed
    
    # Check if there are more scripts left to run
    remaining_scripts = len(all_scripts) - len(executed_scripts)
    if remaining_scripts == 0:
        print("All scripts have been executed.")
        break

    print(f"Executed {len(executed_scripts)} scripts so far. Waiting 10 minutes before running the next batch of {batch_size}.")
    time.sleep(900)  # Wait 15 minutes before running the next batch