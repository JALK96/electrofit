# packages/electrofit/src/electrofit/workflows/run_process_initial_structure.py
import argparse
import os
import logging

from electrofit.config.loader import load_config, dump_config
# Legacy import path; logic now in electrofit.domain.prep.process_initial
from electrofit.core.process_initial_structure import process_initial_structure
from electrofit.io.files import find_file_with_extension

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    args = ap.parse_args()
    project_root = os.path.abspath(args.project)

    run_dir = os.getcwd()

    # Config resolution: prefer run-local snapshot; external overrides are injected earlier via builder in step1.
    run_local_cfg = os.path.join(run_dir, "electrofit.toml")
    cfg_path = run_local_cfg if os.path.isfile(run_local_cfg) else None

    cfg = load_config(project_root, config_path=cfg_path)
    dump_config(cfg, log_fn=logging.info)

    # ---- decide scratch dir with fallbacks ----
    base_scratch_dir = (
        getattr(cfg.paths, "base_scratch_dir", None)
        or os.environ.get("ELECTROFIT_SCRATCH_DIR")
        or "/tmp/electrofit_scratch"
    )

    # Figure out mol2 & molecule name
    mol2_from_name = None
    if cfg.project.molecule_name:
        cand = os.path.join(run_dir, f"{cfg.project.molecule_name}.mol2")
        if os.path.isfile(cand):
            mol2_from_name = cand

    mol2_file = mol2_from_name or find_file_with_extension("mol2")
    if not mol2_file:
        print("[step1][skip] no mol2 file detected in run directory; nothing to process")
        return
    molecule_name = os.path.splitext(os.path.basename(mol2_file))[0]

    additional_input = []
    if cfg.project.adjust_symmetry:
        try:
            json_file = find_file_with_extension("json")
            if json_file:
                additional_input.append(json_file)
        except Exception:
            pass

    process_initial_structure(
        molecule_name=molecule_name,
        mol2_file=os.path.basename(mol2_file),
        base_scratch_dir=base_scratch_dir,
        additional_input=additional_input,
        residue_name=cfg.project.residue_name or molecule_name[:3].upper(),
        net_charge=cfg.project.charge if cfg.project.charge is not None else 0,
        adjust_sym=cfg.project.adjust_symmetry,
        ignore_sym=cfg.project.ignore_symmetry,
        atom_type=cfg.project.atom_type or "gaff2",
        exit_screen=False,
        protocol=cfg.project.protocol or "bcc",
    )

if __name__ == "__main__":
    main()