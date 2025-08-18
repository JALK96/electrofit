"""Deprecated shim for conformer charge processing.

Real implementation moved to :mod:`electrofit.domain.charges.process_conformer`.
This module keeps the original function signature for backwards compatibility
and emits a one-time :class:`DeprecationWarning` mirroring the strategy used in
``core.process_initial_structure``.

Notes:
- Single warning per process (``_DEPRECATION_WARNED`` flag) for noise control.
- No additional injection points yet; domain layer already factors sub‑steps.
- If future tests require monkeypatching of adapters, re-export symbols here.
"""
from __future__ import annotations

import logging
import warnings
import os
import subprocess
from electrofit.cli.safe_run import ensure_finalized
from electrofit.infra.scratch_manager import setup_scratch_directory
from electrofit.domain.charges.process_conformer import (
    ConformerConfig,
    process_conformer,
)

__all__ = ["process_conform"]

PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
project_path = PROJECT_PATH


_DEPRECATION_WARNED = False


def process_conform(molecule_name: str, pdb_file: str, base_scratch_dir: str, net_charge: int, residue_name: str,
                    adjust_sym: bool = False, ignore_sym: bool = False, exit_screen: bool = True, protocol: str = "bcc"):
    """Legacy entrypoint delegating to domain layer (deprecated).

    Parameters mirror historical signature; new code should construct a
    :class:`ConformerConfig` and call :func:`electrofit.domain.charges.process_conformer`.
    """
    global _DEPRECATION_WARNED
    if not _DEPRECATION_WARNED:
        warnings.warn(
            "electrofit.core.process_conform.process_conform is deprecated; use electrofit.domain.charges.process_conformer.process_conformer instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        _DEPRECATION_WARNED = True
    from electrofit.io.files import find_file_with_extension  # local import to avoid heavy module cost for non-users

    if protocol == "opt":  # keep previous behaviour for selecting input files
        respin1_base = "ANTECHAMBER_RESP1_MOD.IN" if adjust_sym else "ANTECHAMBER_RESP1.IN"
        input_files: list[str] = [pdb_file, respin1_base, "ANTECHAMBER_RESP2.IN"]
    else:
        if adjust_sym:
            json_file = find_file_with_extension("json")
            input_files = [pdb_file, json_file] if json_file else [pdb_file]
        else:
            input_files = [pdb_file]

    scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch_dir)
    logging.info("Following %s protocol.", protocol)
    defer_finalize = os.environ.get("ELECTROFIT_DEBUG_DEFER_FINALIZE") == "1"
    context_mgr = None if defer_finalize else ensure_finalized(
        original_dir=original_dir,
        scratch_dir=scratch_dir,
        input_files=input_files,
        remove_parent_if_empty=False,
    )
    cfg = ConformerConfig(
        molecule_name=molecule_name,
        pdb_file=pdb_file,
        net_charge=net_charge,
        residue_name=residue_name,
        adjust_sym=adjust_sym,
        ignore_sym=ignore_sym,
        protocol=protocol,
    )
    try:
        if context_mgr is not None:
            with context_mgr:
                process_conformer(cfg, scratch_dir, original_dir, input_files, defer_finalize=False)
        else:
            process_conformer(cfg, scratch_dir, original_dir, input_files, defer_finalize=True)
        if exit_screen:
            try:
                os.chdir(original_dir)
            except Exception:  # pragma: no cover - defensive
                pass
            sty = os.environ.get("STY")
            if sty:
                try:
                    subprocess.run(["screen", "-S", sty, "-X", "quit"], check=True)
                except subprocess.CalledProcessError:  # pragma: no cover - environment dependent
                    logging.warning("Failed to quit screen session %s", sty)
    except Exception as e:  # keep broad except for legacy parity
        logging.error("Error processing conform: %s", e)
        raise
    logging.debug("[process-conform-end] %s", molecule_name)
