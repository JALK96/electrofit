import sys
import os
import argparse

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

def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description='Edit RESP input file with equivalence groups.')
    
    # Define required positional arguments
    parser.add_argument('input_RESP_file', help='Path to the input RESP file')
    parser.add_argument('equiv_groups_file', help='Path to the equivalence groups JSON file')
    parser.add_argument('output_RESP_file', help='Path to save the output RESP file')
    
    # Define optional arguments
    parser.add_argument('--ignore_sym', action='store_true', 
                        help='Ignore symmetry groups when editing RESP file')
    
    # Parse the arguments
    args = parser.parse_args()
    
    # Call the edit_resp_input function with parsed arguments
    edit_resp_input(
        input_file=args.input_RESP_file,
        equiv_groups_file=args.equiv_groups_file,
        output_file=args.output_RESP_file,
        ignore_sym=args.ignore_sym
    )

if __name__ == "__main__":
    main()