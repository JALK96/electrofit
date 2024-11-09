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

from electrofit.helper.file_manipulation import parse_mol2, build_connections, create_equivalence_groups, write_equivalence_groups


def main():
    if len(sys.argv) != 3:
        print("Usage: python create_equiv_groups.py <input_MOL2_file> <output_JSON_file>")
        print("Example: python create_equiv_groups.py molecule.mol2 equiv_groups.json")
        sys.exit(1)

    mol2_file = sys.argv[1]
    output_json = sys.argv[2]

    atoms, bonds = parse_mol2(mol2_file)
    build_connections(atoms, bonds)
    equiv_groups = create_equivalence_groups(atoms)

    if not equiv_groups:
        print("No equivalence groups found. Please check your MOL2 file's structure.")
    else:
        write_equivalence_groups(equiv_groups, output_json)

if __name__ == "__main__":
    main()