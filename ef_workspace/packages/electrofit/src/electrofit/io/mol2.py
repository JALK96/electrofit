import os
import sys

from openmol import tripos_mol2 as mol2


def read_charges(chg_file_path):
    """
    Reads charges from a .chg file and returns a list of floats.
    Assumes charges are separated by spaces and may span multiple lines.
    """
    charges = []
    try:
        with open(chg_file_path, "r") as file:
            for line in file:
                # Split the line into tokens based on whitespace and convert to float
                line_charges = [float(charge) for charge in line.strip().split()]
                charges.extend(line_charges)
    except FileNotFoundError:
        print(f"Error: Charges file '{chg_file_path}' not found.")
        sys.exit(1)
    except ValueError as ve:
        print(
            f"Error: Non-numeric value encountered in charges file '{chg_file_path}'."
        )
        print(f"Details: {ve}")
        sys.exit(1)
    return charges


def update_mol2_charges(mol2_input_path, chg_file_path, mol2_output_path):
    """
    Updates the charges in a MOL2 file using charges from a .chg file.

    Parameters:
    - mol2_input_path: Path to the input MOL2 file.
    - chg_file_path: Path to the .chg file containing new charges.
    - mol2_output_path: Path where the updated MOL2 file will be saved.
    """
    # Step 1: Check if input MOL2 file exists
    if not os.path.isfile(mol2_input_path):
        print(f"Error: MOL2 input file '{mol2_input_path}' does not exist.")
        sys.exit(1)

    # Step 2: Read the MOL2 file
    print(f"Reading MOL2 file from '{mol2_input_path}'...")
    try:
        p = mol2.read(mol2_input_path)
    except Exception as e:
        print(f"Error: Failed to read MOL2 file '{mol2_input_path}'.")
        print(f"Details: {e}")
        sys.exit(1)

    num_atoms = len(p.atom_name)
    print(f"Number of atoms in MOL2 file: {num_atoms}")

    # Step 3: Read the charges from the .chg file
    print(f"Reading charges from '{chg_file_path}'...")
    charges = read_charges(chg_file_path)
    num_charges = len(charges)
    print(f"Number of charges read: {num_charges}")

    # Step 4: Validate the number of charges matches the number of atoms
    if num_charges != num_atoms:
        print(
            f"Error: Number of charges ({num_charges}) does not match number of atoms ({num_atoms})."
        )
        sys.exit(1)
    else:
        print("Number of charges matches the number of atoms.")

    # Step 5: Update the charges in the MOL2 structure
    print("Updating atom charges...")
    for i in range(num_atoms):
        original_charge = p.atom_q[i]
        p.atom_q[i] = charges[i]
        atom_name = p.atom_name[i]
        print(
            f"Atom {i + 1}: {atom_name} - Charge updated from {original_charge} to {charges[i]}"
        )

    # Step 6: Rebuild the MOL2 structure
    print("Rebuilding MOL2 structure...")
    try:
        p = mol2.build(p)
    except Exception as e:
        print("Error: Failed to rebuild MOL2 structure.")
        print(f"Details: {e}")
        sys.exit(1)

    # Step 7: Write the updated MOL2 file
    print(f"Writing updated MOL2 file to '{mol2_output_path}'...")
    try:
        mol2.Writer(p, mol2_output_path).write()
    except Exception as e:
        print(f"Error: Failed to write MOL2 file to '{mol2_output_path}'.")
        print(f"Details: {e}")
        sys.exit(1)
    print("MOL2 file update complete.")


def main():
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 4:
        print(
            "Usage: python update_mol2.py <input_MOL2_file> <chg_file> <output_MOL2_file>"
        )
        print(
            "Example: python update_mol2.py mol2_input.mol2 charges.chg fixed_out.mol2"
        )
        sys.exit(1)

    # Assign command-line arguments to variables
    mol2_file = sys.argv[1]
    chg_file = sys.argv[2]
    mol2_output = sys.argv[3]

    # Validate input file paths
    if not os.path.isfile(mol2_file):
        print(f"Error: Input MOL2 file '{mol2_file}' does not exist.")
        sys.exit(1)
    if not os.path.isfile(chg_file):
        print(f"Error: Charges file '{chg_file}' does not exist.")
        sys.exit(1)

    # Ensure the output directory exists
    output_dir = os.path.dirname(os.path.abspath(mol2_output))
    if output_dir and not os.path.exists(output_dir):
        print(
            f"Error: The directory for the output file '{output_dir}' does not exist."
        )
        sys.exit(1)

    # Update the MOL2 charges
    update_mol2_charges(mol2_file, chg_file, mol2_output)


if __name__ == "__main__":
    main()
