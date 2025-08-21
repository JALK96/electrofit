"""Deprecated shim for initial structure processing.

Real implementation moved to `electrofit.domain.prep.process_initial`.
This module keeps the original function signature for backwards compatibility.
"""
from __future__ import annotations

import logging
import warnings
import os
from electrofit.infra.scratch_manager import setup_scratch_directory
from electrofit.cli.safe_run import ensure_finalized
from electrofit.domain.prep.process_initial import (
    InitialPrepConfig,
    process_initial,
)
from electrofit.cli.run_commands import run_acpype  # re-export for legacy test monkeypatching

__all__ = ["process_initial_structure", "run_acpype"]
from electrofit.io.files import find_file_with_extension

PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
project_path = PROJECT_PATH


_DEPRECATION_WARNED = False


def process_initial_structure(molecule_name, mol2_file, base_scratch_dir, additional_input, net_charge, residue_name,
                              adjust_sym=False, ignore_sym=False, atom_type="gaff2", protocol="bcc"):
    global _DEPRECATION_WARNED
    if not _DEPRECATION_WARNED:  # only warn once per process
        warnings.warn(
            "electrofit.core.process_initial_structure.process_initial_structure is deprecated; use electrofit.domain.prep.process_initial.process_initial instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        _DEPRECATION_WARNED = True
    # Prepare scratch inputs (legacy behaviour retains additional_input inclusion).
    # Residue normalization now handled inside domain.prep.process_initial; keep removed here to avoid duplicate work.
    input_files = [mol2_file] + list(additional_input)
    scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch_dir)
    cfg = InitialPrepConfig(
        molecule_name=molecule_name,
        mol2_file=mol2_file,
        net_charge=net_charge,
        residue_name=residue_name,
        atom_type=atom_type,
        adjust_sym=adjust_sym,
        ignore_sym=ignore_sym,
        protocol=protocol,
    )
    try:
        from electrofit.cli.safe_run import ensure_finalized
        with ensure_finalized(original_dir=original_dir, scratch_dir=scratch_dir, input_files=input_files):
            # Forward the (possibly monkeypatched) run_acpype symbol for tests
            process_initial(cfg, scratch_dir, original_dir, input_files, run_acpype=run_acpype)
    except Exception as e:  # pragma: no cover
        logging.error("Error processing initial structure: %s", e)
        raise
