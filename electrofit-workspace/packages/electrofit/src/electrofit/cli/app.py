# packages/electrofit/src/electrofit/cli/app.py
import argparse
import os
import runpy
import sys
import logging
from pathlib import Path

STEP_MODULES = {
    "step0": "electrofit.workflows.step0_setup",
    "step1": "electrofit.workflows.step1_initial_processing",
    "step2": "electrofit.workflows.step2_setup_sim_dir",
    "step3": "electrofit.workflows.step3_start_sim",
    "pis":   "electrofit.workflows.run_process_initial_structure",  # legacy alias
}

def _run_module(module: str, project: Path, config: Path | None, rest: list[str]):
    # Export for legacy modules
    os.environ["ELECTROFIT_PROJECT_PATH"] = str(project)
    if config:
        os.environ["ELECTROFIT_CONFIG_PATH"] = str(config)

    prev_argv = sys.argv
    try:
        # Only pass --project; DO NOT pass --config to legacy modules
        sys.argv = [module, "--project", str(project)]
        # Forward any extra user args (if any)
        sys.argv += rest
        runpy.run_module(module, run_name="__main__")
    finally:
        sys.argv = prev_argv

def _ensure_project_scaffold(project: str | os.PathLike) -> tuple[Path, Path]:
    root = Path(project).resolve()
    process = root / "process"
    if not process.exists():
        logging.info(f"'process' folder not found at {process}; creating it for you.")
        process.mkdir(parents=True, exist_ok=True)
    return root, process

def main():
    # Be talkative by default unless caller configured logging already
    if not logging.getLogger().handlers:
        logging.basicConfig(level=logging.INFO, format="%(message)s")

    p = argparse.ArgumentParser("electrofit")
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_cmd(name: str):
        sp = sub.add_parser(name)
        sp.add_argument("--project", required=True, help="Path to project/case folder")
        sp.add_argument("--config", help="Optional path to electrofit.toml")
        return sp

    for cmd in STEP_MODULES:
        add_cmd(cmd)

    args, rest = p.parse_known_args()

    # Ensure project scaffold (create <project>/process if missing)
    project_root, _process_dir = _ensure_project_scaffold(args.project)

    # Resolve config: use provided file, otherwise fallback to CWD/electrofit.toml (where the user is running)
    if args.config:
        config_path = Path(args.config).resolve()
    else:
        default_cfg = Path.cwd() / "electrofit.toml"
        config_path = default_cfg if default_cfg.exists() else None
        if config_path:
            logging.info(f"No --config given; using {config_path}")

    logging.info(f"Dispatching {args.cmd} | project={project_root}"
                 + (f" | config={config_path}" if config_path else " | config=<none>"))

    _run_module(STEP_MODULES[args.cmd], project_root, config_path, rest)

if __name__ == "__main__":
    main()