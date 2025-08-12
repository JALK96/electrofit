import logging
import os

from electrofit.config_parser import ConfigParser

# trunk-ignore(ruff/E402)
from electrofit.core.process_conform import process_conform
from electrofit.io.files import find_file_with_extension, strip_extension
from electrofit.logging import setup_logging


def main_conform_processing():
    """
    Main function to process initial structure and create GROMACS input.
    """

    # Find current working directory
    home = os.getcwd()

    # Initialize logging with log file in home
    log_file_path = os.path.join(home, "process.log")
    setup_logging(log_file_path)

    logging.info(f"Working Directory: {home}")

    # Get parent directory
    extracted_conforms_dir = os.path.dirname(home)
    logging.info(f"Extracted Conform Directory: {extracted_conforms_dir}")

    # Go to the parent folder (the extracted_conforms directory) to open the configuration file (.ef)
    os.chdir(extracted_conforms_dir)

    input = find_file_with_extension("ef")
    logging.info(
        f"=== Reading parameter from configuration file: {input} in {extracted_conforms_dir} ==="
    )
    config = ConfigParser(input)

    # Go back to the home directory
    os.chdir(home)

    # Define file and molecule name
    pdb_file = find_file_with_extension("pdb")
    molecule_name = strip_extension(pdb_file)
    logging.info(f"Processing conform: {molecule_name}")
    logging.info("------------------------------------------")

    # Define
    base_scratch_dir = config.BaseScratchDir
    logging.info(f"Scratch directory set to: {base_scratch_dir}")
    residue_name = config.ResidueName
    logging.info(f"Residue Name: {residue_name}")
    net_charge = config.Charge
    logging.info(f"Charge set to: {net_charge}")
    adjust_sym = config.AdjustSymmetry
    logging.info(f"AdjustSymmetry set to: {adjust_sym}")
    protocol = config.Protocol
    logging.info(f"Charge fit protocol set to: {protocol}")
    ignore_sym = config.IgnoreSymmetry
    logging.info(f"IgnoreSymmetry set to: {ignore_sym}")

    logging.info("=== Executing script 'process_conform'! ===")

    # Process the initial structure
    process_conform(
        molecule_name=molecule_name,
        pdb_file=pdb_file,
        base_scratch_dir=base_scratch_dir,
        net_charge=net_charge,
        residue_name=residue_name,
        adjust_sym=adjust_sym,
        protocol=protocol,
        ignore_sym=ignore_sym,
    )


if __name__ == "__main__":
    main_conform_processing()
