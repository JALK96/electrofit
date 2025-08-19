"""(Deprecated) Step2 legacy entrypoint.

Moved to :mod:`electrofit.pipeline.steps.step2`. Use the new pipeline module.

Original docstring:

Step2: Create run_gmx_simulation directories per molecule.

Copies GROMACS input files from *.acpype subdirectories of run_gau_create_gmx_in/, propagates electrofit.toml snapshot,
copies MDP templates and builds a run.json manifest for Step3.

Now uses unified logging & snapshot layering (results > process/<mol> > data/input/<mol> > project).

CLI: electrofit step2 --project <path> [--config override.toml]
"""
import argparse
import fnmatch
import os
import shutil
import json
import logging
from pathlib import Path

import warnings
warnings.warn(
    "Importing 'electrofit.workflows.step2_setup_sim_dir' is deprecated; use 'electrofit.pipeline.steps.step2' instead.",
    DeprecationWarning,
    stacklevel=2,
)
from electrofit.infra.config_snapshot import compose_snapshot, CONFIG_ARG_HELP
from electrofit.infra.logging import setup_logging, log_run_header, reset_logging


def _write_manifest(dest_dir: str, files: dict[str, str], mdp_subdir: str = "MDP") -> None:
    """Write a small run.json manifest so step3 can run deterministically."""
    gro = files.get("gro")
    molecule = None
    if gro:
        base = os.path.splitext(os.path.basename(gro))[0]
        molecule = base[:-4] if base.endswith("_GMX") else base
    def _basename(val):
        return os.path.basename(val) if isinstance(val, str) and val else ""
    manifest = {
        "molecule": molecule or "unknown",
        "gro": _basename(files.get("gro")),
        "top": _basename(files.get("top")),
        "itp": _basename(files.get("itp")),
        "posres": _basename(files.get("posres")),
        "mdp_dir": mdp_subdir,
    }
    with open(os.path.join(dest_dir, "run.json"), "w") as f:
        json.dump(manifest, f, indent=2)
    logging.info(f"[step2] Wrote manifest: {os.path.join(dest_dir, 'run.json')}")


def main():  # pragma: no cover (CLI wrapper)
    ap = argparse.ArgumentParser(description="Step2: Prepare run_gmx_simulation directories")
    ap.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--log-console", action="store_true", help="Also echo logs to console")
    args = ap.parse_args()

    project_path = Path(args.project).resolve()
    process_dir = project_path / "process"
    mdp_source_dir = project_path / "data" / "MDP"
    if not process_dir.is_dir():
        print("[step2] No process directory.")
        return
    mol_dirs = [p for p in process_dir.iterdir() if p.is_dir()]
    multi_mol = len(mol_dirs) > 1

    # Global config (for defaults & echo at top-level)
    top_log = project_path / "step.log"
    setup_logging(str(top_log), also_console=args.log_console)
    log_run_header("step2")

    file_patterns = ["*GMX.gro", "*GMX.itp", "*GMX.top", "posre_*.itp"]

    prepared = 0
    for folder_path in mol_dirs:
        run_gau_dir = folder_path / "run_gau_create_gmx_in"
        if not run_gau_dir.is_dir():
            logging.info(f"[step2][skip] {folder_path.name}: no run_gau_create_gmx_in dir")
            continue
        dest_dir = folder_path / "run_gmx_simulation"
        dest_dir.mkdir(exist_ok=True)
        logging.info(f"[step2] Preparing {folder_path.name} -> {dest_dir}")
        # Per-molecule logging
        reset_logging()
        setup_logging(str(dest_dir / "process.log"), also_console=False)
        log_run_header("step2")
        upstream_snap = run_gau_dir / "electrofit.toml"
        molecule_input = project_path / "data" / "input" / folder_path.name / "electrofit.toml"
        process_cfg = folder_path / "electrofit.toml"
        project_defaults = project_path / "electrofit.toml"
        extra_override = Path(args.config) if args.config else None
        snap_dest = compose_snapshot(
            dest_dir,
            project_path,
            folder_path.name,
            multi_molecule=multi_mol,
            log_fn=logging.info,
            upstream=upstream_snap,
            process_cfg=process_cfg,
            molecule_input=molecule_input,
            project_defaults=project_defaults,
            extra_override=extra_override,
        )
        # Locate acpype
        acpype_dir = None
        for sub in run_gau_dir.iterdir():
            if sub.is_dir() and sub.name.endswith(".acpype"):
                acpype_dir = sub; break
        if not acpype_dir:
            logging.warning(f"[step2][skip] {folder_path.name}: no .acpype dir")
            continue
        # Copy matched topology/input files
        for pattern in file_patterns:
            for file_name in os.listdir(acpype_dir):
                if fnmatch.fnmatch(file_name, pattern):
                    shutil.copy(acpype_dir / file_name, dest_dir)
                    logging.info(f"[step2] Copied {file_name} -> {dest_dir}")
        # Copy MDP dir
        md_dest_dir = dest_dir / "MDP"
        if mdp_source_dir.is_dir():
            if md_dest_dir.exists():
                shutil.rmtree(md_dest_dir)
            shutil.copytree(mdp_source_dir, md_dest_dir)
            logging.info(f"[step2] Copied MDP -> {md_dest_dir}")
        else:
            logging.warning(f"[step2][warn] no MDP source {mdp_source_dir}")
        # Build manifest
        selected = {"gro": None, "itp": None, "top": None, "posres": None}
        for name in os.listdir(dest_dir):
            if name.endswith("GMX.gro"):
                selected["gro"] = os.path.join(dest_dir, name)
            elif name.endswith("GMX.itp") and not name.startswith("posre_"):
                selected["itp"] = os.path.join(dest_dir, name)
            elif name.endswith(".top"):
                selected["top"] = os.path.join(dest_dir, name)
            elif name.startswith("posre_") and name.endswith(".itp"):
                selected["posres"] = os.path.join(dest_dir, name)
        # Fallback scans
        if not selected["gro"]:
            for name in os.listdir(dest_dir):
                if name.endswith(".gro"):
                    selected["gro"] = os.path.join(dest_dir, name); break
        if not selected["itp"]:
            for name in os.listdir(dest_dir):
                if name.endswith(".itp") and not name.startswith("posre_"):
                    selected["itp"] = os.path.join(dest_dir, name); break
        if not selected["top"]:
            for name in os.listdir(dest_dir):
                if name.endswith(".top"):
                    selected["top"] = os.path.join(dest_dir, name); break
        if not selected["posres"]:
            for name in os.listdir(dest_dir):
                if name.startswith("posre_") and name.endswith(".itp"):
                    selected["posres"] = os.path.join(dest_dir, name); break
        _write_manifest(dest_dir, selected, mdp_subdir="MDP")
        prepared += 1
        reset_logging()

    print(f"[step2] Prepared simulation dirs for {prepared}/{len(mol_dirs)} molecules.")


if __name__ == "__main__":  # pragma: no cover
    main()