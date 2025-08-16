# packages/electrofit/src/electrofit/cli/app.py
import argparse
import os
import runpy
import importlib
import sys
import logging
from pathlib import Path

STEP_MODULES = {
    "step0": "electrofit.workflows.step0_setup",
    "step1": "electrofit.workflows.step1_initial_processing",
    "step2": "electrofit.workflows.step2_setup_sim_dir",
    "step3": "electrofit.workflows.step3_start_sim",
    "step4": "electrofit.workflows.step4_extract_conforms",
    "step5": "electrofit.workflows.step5_process_conforms",
    "step6": "electrofit.workflows.step6_extract_average_charges",
    "step7": "electrofit.workflows.step7_setup_final_sim",
    "step8": "electrofit.workflows.step8_start_final_sim",
    "pis":   "electrofit.workflows.run_process_initial_structure",  # legacy alias
}

def _run_module(module: str, project: Path, config: Path | None, rest: list[str]):
    # Export project path only (old ELECTROFIT_CONFIG_PATH removed)
    os.environ["ELECTROFIT_PROJECT_PATH"] = str(project)

    prev_argv = sys.argv
    try:
        # Only pass --project; DO NOT pass --config to legacy modules
        sys.argv = [module, "--project", str(project)] + rest
        if module == "electrofit.workflows.step4_extract_conforms":
            # Import normally so multiprocessing can pickle top-level functions
            mod = importlib.import_module(module)
            if hasattr(mod, "main"):
                mod.main()
            else:
                raise RuntimeError(f"Module {module} has no main()")
        else:
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

    DESCRIPTIONS = {
        'step0': 'Initial project setup: copy input molecule folders from data/input -> process (one directory per molecule).',
        'step1': 'Initial structure processing (Gaussian/Antechamber preparation) for all run_gau_create_gmx_in directories.',
        'step2': 'Create run_gmx_simulation directories, copy GMX files & MDP templates, write run.json manifest.',
        'step3': 'Launch GROMACS production setup (box, ions, topology) per molecule based on run.json.',
        'step4': 'Extract conformers from simulation trajectories (sampling: linear|random|maxmin) into extracted_conforms/.',
        'step5': 'Process extracted conformers (Gaussian/RESP pipeline) with batching & parallelisation options.',
    }
    def add_cmd(name: str):
        sp = sub.add_parser(name, help=DESCRIPTIONS.get(name,''), description=DESCRIPTIONS.get(name,''))
        sp.add_argument("--project", required=True, help="Path to project/case folder")
        sp.add_argument("--config", help="Optional path to electrofit.toml")
        return sp

    for cmd in STEP_MODULES:
        add_cmd(cmd)

    args, rest = p.parse_known_args()
    # Flag-only validation: if user mistakenly supplies a value after --isolate-conformer
    if args.cmd == 'step5' and '--isolate-conformer' in sys.argv:
        # finde Position
        for i, tok in enumerate(sys.argv):
            if tok == '--isolate-conformer':
                # if next token does not start with '-' and is not a subcommand -> warning
                if i+1 < len(sys.argv) and not sys.argv[i+1].startswith('-'):
                    print("[step5][warn] '--isolate-conformer' is a flag-only option; ignored value: '"+sys.argv[i+1]+"'", file=sys.stderr)
                break

    # Ensure project scaffold (create <project>/process if missing)
    project_root, _process_dir = _ensure_project_scaffold(args.project)

    # Resolve config path: explicit only (no legacy env fallback anymore).
    if args.config:
        config_path = Path(args.config).resolve()
    else:
        config_path = None

    logging.info(f"Dispatching {args.cmd} | project={project_root}"
                 + (f" | config={config_path}" if config_path else " | config=<none>"))

    _run_module(STEP_MODULES[args.cmd], project_root, config_path, rest)

if __name__ == "__main__":
    main()