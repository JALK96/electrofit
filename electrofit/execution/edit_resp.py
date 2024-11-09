import sys
import os

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

from electrofit.helper.file_manipulation import edit_resp_input

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python edit_resp.py <input_RESP_file> <equiv_groups_file> <output_RESP_file>")
        sys.exit(1)
    edit_resp_input(sys.argv[1], sys.argv[2], sys.argv[3])