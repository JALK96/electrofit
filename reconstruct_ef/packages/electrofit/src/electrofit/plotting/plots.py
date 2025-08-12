import sys
import os
import argparse
import json

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

from electrofit.helper.file_manipulation import load_symmetry_groups
from electrofit.helper.plotting import plot_charges_by_atom, plot_charges_by_symmetry


import sys
import os
import argparse
import json
import matplotlib.pyplot as plt
import re
import numpy as np

def main():
    # Set up the argument parser
    parser = argparse.ArgumentParser(description='Plot charges by atom or symmetry.')
    parser.add_argument('-ic', '--initial_charges', required=True, help='Path to initial_charges_dict.json')
    parser.add_argument('-c', '--charges', required=True, help='Path to charges_dict.json')
    parser.add_argument('-c2', '--charges2', help='Path to second charges_dict.json for comparison (optional)')
    parser.add_argument('-s', '--symmetry_groups', help='Path to symmetry_groups.json (optional)')
    parser.add_argument('-d', '--directory', default='.', help='Directory to save the plot (optional)')

    # Parse the arguments
    args = parser.parse_args()

    # Load the initial charges dictionary
    with open(args.initial_charges, 'r') as f:
        initial_charges_dict = json.load(f)

    # Load the first charges dictionary
    with open(args.charges, 'r') as f:
        atoms_dict1 = json.load(f)

    # Load the second charges dictionary if provided
    if args.charges2:
        with open(args.charges2, 'r') as f:
            atoms_dict2 = json.load(f)
    else:
        atoms_dict2 = None

    # Ensure the output directory exists
    base_dir = args.directory
    os.makedirs(base_dir, exist_ok=True)

    # Decide which plotting function to use
    if args.symmetry_groups:
        # Load symmetry groups if provided
        symmetry_groups = load_symmetry_groups(args.symmetry_groups)
        # Plot charges by symmetry
        plot_charges_by_symmetry(atoms_dict1, initial_charges_dict, base_dir, symmetry_groups, atoms_dict2)
    else:
        # Plot charges by atom
        plot_charges_by_atom(atoms_dict1, initial_charges_dict, base_dir, atoms_dict2)

if __name__ == '__main__':
    main()