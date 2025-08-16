# packages/electrofit/src/electrofit/workflows/step1_initial_processing.py
import argparse
import os
import logging
from pathlib import Path
from contextlib import contextmanager

from electrofit.workflows.run_process_initial_structure import main as run_pis_main
from electrofit.workflows.snapshot import build_snapshot_with_layers, CONFIG_ARG_HELP
from electrofit.logging import setup_logging, log_run_header, reset_logging


TARGET_RUN_DIRS = {
    "run_gau_create_gmx_in",        # legacy deterministic name
    "init_system_create_gmx_in",    # future/clearer name (alias)
}


@contextmanager
def _pushd(path: str | None):
    """Temporarily change working directory."""
    prev = os.getcwd()
    if path:
        os.chdir(path)
    try:
        yield
    finally:
        if path:
            os.chdir(prev)


def _iter_run_dirs(project_root: str):
    """
    Yield run directories under '<project>/process/**/<TARGET_RUN_DIRS>'.
    """
    process_root = os.path.join(project_root, "process")
    if not os.path.isdir(process_root):
        return
    for root, dirs, _files in os.walk(process_root):
        for d in dirs:
            if d in TARGET_RUN_DIRS:
                yield os.path.join(root, d)


def _run_one_dir(run_dir: str, project_root: str, override_cfg: str | None, multi_mol: bool):
    """Execute initial structure processing in one run directory with snapshot merge."""
    os.environ["ELECTROFIT_PROJECT_PATH"] = project_root
    mol = Path(run_dir).parent.name
    # logging per run directory
    reset_logging()
    setup_logging(str(Path(run_dir) / "process.log"), also_console=False)
    log_run_header("step1")
    molecule_input = Path(project_root) / "data" / "input" / mol / "electrofit.toml"
    project_defaults = Path(project_root) / "electrofit.toml"
    snap = build_snapshot_with_layers(
        Path(run_dir),
        Path(project_root),
        mol,
        multi_molecule=multi_mol,
        log_fn=logging.info,
        molecule_input=molecule_input,
        project_defaults=project_defaults,
        extra_override=Path(override_cfg) if override_cfg else None,
    )

    with _pushd(run_dir):
        run_pis_main()


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Step1: Initial molecule processing (prepare Gaussian/Antechamber inputs). "
            "Discovers run_gau_create_gmx_in directories under process/ and creates/merges electrofit.toml snapshots."
        )
    )
    ap.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help=CONFIG_ARG_HELP)
    ap.add_argument("--log-console", action="store_true", help="Also echo logging to console")
    args = ap.parse_args()

    project_root = os.path.abspath(args.project)
    process_root = Path(project_root) / "process"
    multi_mol = False
    if process_root.is_dir():
        mol_dirs = [p for p in process_root.iterdir() if p.is_dir()]
        multi_mol = len(mol_dirs) > 1

    # 1) Deterministic discovery under project/process/**/<TARGET_RUN_DIRS>
    discovered = list(_iter_run_dirs(project_root))

    if discovered:
        for run_dir in discovered:
            _run_one_dir(run_dir, project_root, args.config, multi_mol)
        return

    # 2) Fallback: treat the current directory as a single run (keeps tests simple)
    run_dir = os.getcwd()
    _run_one_dir(run_dir, project_root, args.config, multi_mol)


if __name__ == "__main__":
    main()