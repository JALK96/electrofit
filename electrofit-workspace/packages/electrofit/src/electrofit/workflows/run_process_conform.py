import logging
import os
from pathlib import Path

from electrofit.config.loader import load_config, dump_config
# Legacy path; actual logic resides in electrofit.domain.charges.process_conformer
from electrofit.core.process_conform import process_conform
from electrofit.io.files import find_file_with_extension, strip_extension
from electrofit.infra.logging import setup_logging


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

    # Determine project root (assumes run directory structure: .../process/<mol>/extracted_conforms/<confX>)
    run_dir = Path(home)
    # ascend until we hit 'process'
    p = run_dir
    project_root = None
    molecule_name = None
    while p.parent != p:
        if p.name == "process":
            project_root = p.parent
            break
        p = p.parent
    if project_root is None:
        project_root = run_dir  # fallback
    # molecule directory is immediate child under process
    try:
        molecule_name = run_dir.relative_to(project_root / "process").parts[0]
    except Exception:
        molecule_name = None

    cfg = load_config(project_root, context_dir=run_dir, molecule_name=molecule_name)
    dump_config(cfg, log_fn=logging.info)
    proj_cfg = cfg.project

    # Define file and molecule name
    pdb_file = find_file_with_extension("pdb")
    molecule_name = proj_cfg.molecule_name or strip_extension(pdb_file)
    logging.info(f"Processing conform: {molecule_name}")
    logging.info("------------------------------------------")

    # Define
    base_scratch_dir = cfg.paths.base_scratch_dir or "/tmp/electrofit_scratch"
    logging.info(f"Scratch directory set to: {base_scratch_dir}")
    residue_name = proj_cfg.residue_name or "LIG"
    logging.info(f"Residue Name: {residue_name}")
    net_charge = proj_cfg.charge or 0
    logging.info(f"Charge set to: {net_charge}")
    adjust_sym = getattr(proj_cfg, "adjust_symmetry", False)
    logging.info(f"AdjustSymmetry set to: {adjust_sym}")
    protocol = proj_cfg.protocol or "bcc"
    logging.info(f"Charge fit protocol set to: {protocol}")
    ignore_sym = getattr(proj_cfg, "ignore_symmetry", False)
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
