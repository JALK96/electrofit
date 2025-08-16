import os
import sys
import logging

# Delay heavy/native import until actually needed (for crash isolation)
try:
    from openmol import tripos_mol2 as mol2  # type: ignore
except Exception as _imp_err:  # pragma: no cover - diagnostic path
    mol2 = None  # will trigger fallback/explicit error when used
    logging.warning("Deferred openmol import failed at module import: %s", _imp_err)


class Mol2ChargeError(Exception):
    """Domain-specific exception for MOL2 charge update failures."""


def read_charges(chg_file_path):
    """
    Reads charges from a .chg file and returns a list of floats.
    Assumes charges are separated by spaces and may span multiple lines.
    """
    charges = []
    try:
        with open(chg_file_path, "r") as file:
            for line in file:
                line = line.strip()
                if not line:
                    continue
                try:
                    line_charges = [float(charge) for charge in line.split()]
                except ValueError as ve:
                    raise Mol2ChargeError(
                        f"Non-numeric value in charges file '{chg_file_path}': {ve}"
                    ) from ve
                charges.extend(line_charges)
    except FileNotFoundError as fnf:
        raise Mol2ChargeError(f"Charges file '{chg_file_path}' not found") from fnf
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
    logging.info("[mol2-update] ENTER mol2_input=%s chg=%s out=%s", mol2_input_path, chg_file_path, mol2_output_path)
    if not os.path.isfile(mol2_input_path):
        raise Mol2ChargeError(f"MOL2 input file '{mol2_input_path}' does not exist")

    # Step 2: Read the MOL2 file
    logging.info("[mol2-update] Reading MOL2 file ...")
    if mol2 is None:
        raise Mol2ChargeError("openmol.tripos_mol2 not available (import failed earlier)")
    try:
        p = mol2.read(mol2_input_path)
    except Exception as e:
        raise Mol2ChargeError(f"Failed to read MOL2 file '{mol2_input_path}': {e}") from e

    num_atoms = len(p.atom_name)
    logging.info("[mol2-update] Atoms=%d", num_atoms)

    # Step 3: Read the charges from the .chg file
    logging.info("[mol2-update] Reading charges ...")
    charges = read_charges(chg_file_path)
    num_charges = len(charges)
    logging.info("[mol2-update] Charges read=%d", num_charges)

    # Step 4: Validate the number of charges matches the number of atoms
    if num_charges != num_atoms:
        raise Mol2ChargeError(
            f"Charge/atom count mismatch: charges={num_charges} atoms={num_atoms}"
        )
    else:
        logging.info("[mol2-update] Count match OK")

    # Step 5: Update the charges in the MOL2 structure
    logging.info("[mol2-update] Updating atom charges ...")
    for i in range(num_atoms):
        original_charge = p.atom_q[i]
        p.atom_q[i] = charges[i]
        atom_name = p.atom_name[i]
        logging.debug(
            "[mol2-update] atom=%d name=%s old=%s new=%s",
            i + 1,
            atom_name,
            original_charge,
            charges[i],
        )

    # Step 6: Rebuild the MOL2 structure
    logging.info("[mol2-update] Rebuilding MOL2 structure ...")
    try:
        p = mol2.build(p)
    except Exception as e:
        raise Mol2ChargeError(f"Failed to rebuild MOL2 structure: {e}") from e

    # Step 7: Write the updated MOL2 file
    logging.info("[mol2-update] Writing updated MOL2 file ...")
    try:
        mol2.Writer(p, mol2_output_path).write()
    except Exception as e:
        raise Mol2ChargeError(f"Failed to write MOL2 file '{mol2_output_path}': {e}") from e
    logging.info("[mol2-update] DONE")


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
        raise Mol2ChargeError(f"Input MOL2 file '{mol2_file}' does not exist")
    if not os.path.isfile(chg_file):
        raise Mol2ChargeError(f"Charges file '{chg_file}' does not exist")

    # Ensure the output directory exists
    output_dir = os.path.dirname(os.path.abspath(mol2_output))
    if output_dir and not os.path.exists(output_dir):
        raise Mol2ChargeError(f"Output directory '{output_dir}' does not exist")

    # Update the MOL2 charges
    try:
        update_mol2_charges(mol2_file, chg_file, mol2_output)
    except Mol2ChargeError as e:
        logging.error("mol2 charge update failed: %s", e)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
