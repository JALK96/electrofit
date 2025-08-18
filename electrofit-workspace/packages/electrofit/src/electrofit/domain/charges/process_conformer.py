"""Conformer charge computation orchestration (extracted from legacy core.process_conform).

Initial extraction: collapses duplicated logic branches (normal vs defer finalize) into a reusable
`process_conformer` function plus supporting helpers. Further decomposition (Gaussian/RESP adapters)
can follow in later commits.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import logging
import os
import shutil
import time

from electrofit.cli.run_commands import (
    run_command,
    run_espgen,
    run_python,
)
from electrofit.io.files import (
    find_file_with_extension,
    pdb_to_mol2,
    replace_charge_in_ac_file,
)
from electrofit.io.resp import edit_resp_input
from electrofit.io.mol2 import update_mol2_charges
from electrofit.domain.symmetry.write_symmetry import write_symmetry
from electrofit.domain.symmetry.resp_constraints import apply_and_optionally_modify, RespSymmetryConfig
from electrofit.adapters import gaussian as gaussian_adapter
from electrofit.adapters import resp as resp_adapter
from electrofit.adapters.results import ConformerChargeResult, GaussianResult, RespStageResult
from electrofit.infra.step_logging import log_relevant_config
from electrofit.infra.decisions import build_conformer_decision


@dataclass(slots=True)
class ConformerConfig:
    molecule_name: str
    pdb_file: str
    net_charge: int
    residue_name: str
    adjust_sym: bool = False
    ignore_sym: bool = False
    protocol: str = "bcc"  # or 'opt'


def _maybe_use_gaussian_cache(molecule_name: str, scratch_dir: str) -> bool:
    """Delegate cache hydration to gaussian adapter (kept for backward name)."""
    return gaussian_adapter.maybe_use_cache(molecule_name, scratch_dir)


def _prepare_inputs(cfg: ConformerConfig, scratch_dir: str) -> tuple[str, str]:
    mol2_file = f"{cfg.molecule_name}.mol2"
    pdb_to_mol2(cfg.pdb_file, mol2_file, cfg.residue_name, cwd=scratch_dir)
    gaussian_input = gaussian_adapter.build_input(mol2_file, cfg.molecule_name, cfg.net_charge, scratch_dir)
    return mol2_file, gaussian_input


def _generate_esp(molecule_name: str, scratch_dir: str):
    gesp_file = f"{molecule_name}.gesp"; esp_file = f"{molecule_name}.esp"
    run_espgen(gesp_file, esp_file, scratch_dir)
    return esp_file


def _prepare_resp_inputs(cfg: ConformerConfig, scratch_dir: str) -> tuple[str, str]:
    # Use adapter to produce ac + resp inputs
    resp_adapter.prepare_ac(f"{cfg.molecule_name}.mol2", cfg.molecule_name, cfg.net_charge, scratch_dir)
    resp1, resp2 = resp_adapter.generate_inputs(cfg.molecule_name, scratch_dir)
    if cfg.adjust_sym:
        resp_adapter.apply_symmetry(scratch_dir, cfg.adjust_sym, cfg.ignore_sym)
    return resp1, resp2


def _resolve_resp_inputs(cfg: ConformerConfig, scratch_dir: str) -> tuple[str, str]:
    # Determine effective RESP1 (possibly modified) + RESP2
    resp1 = os.path.join(scratch_dir, "ANTECHAMBER_RESP1_MOD.IN" if cfg.adjust_sym else "ANTECHAMBER_RESP1.IN")
    resp2 = os.path.join(scratch_dir, "ANTECHAMBER_RESP2.IN")
    for f in (resp1, resp2):
        if not os.path.isfile(f):
            raise FileNotFoundError(f"Missing RESP input file: {f}")
    return resp1, resp2


def _write_symmetry_file(respin1: str, adjust: bool, scratch_dir: str):
    """Write symmetry documentation file (purely informational)."""
    run_python(
        write_symmetry,
        respin1,
        "symmetry_resp_MOD.txt" if adjust else "symmetry_resp.txt",
        cwd=scratch_dir,
    )


def process_conformer(cfg: ConformerConfig, scratch_dir: str, original_dir: str, input_files: list[str], defer_finalize: bool) -> ConformerChargeResult:
    """Execute conformer charge pipeline inside prepared scratch directory and return structured result.

    Backward compatibility: Existing callers that ignored the previous (None) return value can continue
    to do so; they now simply discard the returned ConformerChargeResult. No behavioural differences are
    intended vs. the legacy path except for collecting timing metadata.
    """
    # Log minimal relevant config subset & derived decisions once at start
    try:
        log_relevant_config(
            step="step5",
            cfg=cfg,
            fields=[
                "molecule_name",
                "net_charge",
                "protocol",
                "adjust_sym",
                "ignore_sym",
                "residue_name",
            ],
        )
    # Structured decision model (ensemble context)
        try:
            build_conformer_decision(cfg.protocol, cfg.adjust_sym, cfg.ignore_sym).log('step5')
        except Exception:
            logging.debug('[step5][decisions] logging failed', exc_info=True)
    except Exception:  # defensive: logging must not break execution
        logging.debug("[step5][log] selective config logging failed", exc_info=True)

    # 1. Prepare mol2 (Gaussian input will be (re)built inside adapter run_with_result).
    mol2_file, _gaussian_input_legacy = _prepare_inputs(cfg, scratch_dir)

    # 2. Gaussian stage with structured result (includes optional cache hydration)
    gaussian_result: GaussianResult = gaussian_adapter.run_with_result(
        mol2_file=mol2_file,
        name=cfg.molecule_name,
        net_charge=cfg.net_charge,
        cwd=scratch_dir,
        try_cache=True,
    )
    if not gaussian_result.gesp_file or not os.path.isfile(gaussian_result.gesp_file):
        raise FileNotFoundError(f"Gaussian stage did not produce .gesp file for {cfg.molecule_name}")

    # 3. ESP generation from Gaussian .gesp
    esp_file = _generate_esp(cfg.molecule_name, scratch_dir)

    # 4. RESP input preparation (BCC protocol only)
    symmetry_applied = False
    if cfg.protocol == "bcc":
        _prepare_resp_inputs(cfg, scratch_dir)
        symmetry_applied = bool(cfg.adjust_sym)
    respin1, respin2 = _resolve_resp_inputs(cfg, scratch_dir)

    # 5. Symmetry documentation file (before RESP run for better debugging if RESP fails)
    _write_symmetry_file(respin1, cfg.adjust_sym, scratch_dir)
    symmetry_file = os.path.join(
        scratch_dir,
        "symmetry_resp_MOD.txt" if cfg.adjust_sym else "symmetry_resp.txt",
    )

    # 6. RESP two-stage fitting with structured result
    resp_result: RespStageResult = resp_adapter.run_two_stage_with_result(
        name=cfg.molecule_name,
        esp_file=esp_file,
        resp1=respin1,
        resp2=respin2,
        cwd=scratch_dir,
        symmetry_applied=symmetry_applied,
    )

    # 7. Optional deferred finalize
    if defer_finalize:
        from electrofit.infra.scratch_manager import finalize_scratch_directory
        finalize_scratch_directory(
            original_dir,
            scratch_dir,
            input_files,
            output_files=None,
            overwrite=True,
            remove_parent_if_empty=False,
            reason="defer-finalize",
        )
    logging.info("Conformer processing completed for %s", cfg.molecule_name)
    return ConformerChargeResult(
        gaussian=gaussian_result,
        resp=resp_result,
        protocol=cfg.protocol,
        residue_name=cfg.residue_name,
        scratch_dir=scratch_dir,
    symmetry_file=symmetry_file if os.path.isfile(symmetry_file) else None,
    )
