# packages/electrofit/src/electrofit/cli/app.py
import argparse
import os
import runpy
from pathlib import Path

STEP_MODULES = {
    "step0": "electrofit.workflows.step0_setup",
    "step1": "electrofit.workflows.step1_initial_processing",
    "step2": "electrofit.workflows.step2_setup_sim_dir",
    "step3": "electrofit.workflows.step3_start_sim",
    "pis":   "electrofit.workflows.run_process_initial_structure",  # keep legacy entry if you like
}

def _run_module(module: str, project: Path, config: Path | None):
    os.environ["ELECTROFIT_PROJECT_PATH"] = str(project)
    if config:
        os.environ["ELECTROFIT_CONFIG_PATH"] = str(config)
    runpy.run_module(module, run_name="__main__")

def main():
    p = argparse.ArgumentParser("electrofit")
    sub = p.add_subparsers(dest="cmd", required=True)

    def add_cmd(name: str):
        sp = sub.add_parser(name)
        sp.add_argument("--project", required=True, help="Path to project/case folder")
        sp.add_argument("--config", help="Optional path to electrofit.toml")
        return sp

    for cmd in STEP_MODULES:
        add_cmd(cmd)

    args = p.parse_args()
    project = Path(args.project).resolve()
    config = Path(args.config).resolve() if args.config else None
    _run_module(STEP_MODULES[args.cmd], project, config)

if __name__ == "__main__":
    main()