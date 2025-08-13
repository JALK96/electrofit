import os

from electrofit.io.files import copy_and_rename_folders

PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
project_path = PROJECT_PATH

# Define source and destination paths
src = os.path.join(project_path, "data/input")
dst = os.path.join(project_path, "process")
bash_src = os.path.join(project_path, "scripts/pis.sh")

print("Source path:", src)  # Debugging line
print("Destination path:", dst)  # Debugging line

# Check if src exists to prevent errors
if not os.path.exists(src):
    raise FileNotFoundError(f"Source directory '{src}' does not exist.")

# Execute the function
copy_and_rename_folders(source=src, destination=dst, bash_script_source=bash_src)
