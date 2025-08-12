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

print("Project root path:", project_path)  # Debugging line

# Add the project root to sys.path for module imports
sys.path.append(project_path)

from electrofit.helper.file_manipulation import copy_and_rename_folders

# Define source and destination paths
src = os.path.join(project_path, "data/input")
dst = os.path.join(project_path, "process")

print("Source path:", src)  # Debugging line
print("Destination path:", dst)  # Debugging line

# Check if src exists to prevent errors
if not os.path.exists(src):
    raise FileNotFoundError(f"Source directory '{src}' does not exist.")

# Execute the function
copy_and_rename_folders(src, dst)