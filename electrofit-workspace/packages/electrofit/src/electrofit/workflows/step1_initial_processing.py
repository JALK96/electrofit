# packages/electrofit/src/electrofit/workflows/step1_initial_processing.py
import argparse
import os
from contextlib import contextmanager

from electrofit.workflows.run_process_initial_structure import main as run_pis_main


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


def _pick_config(run_dir: str, project_root: str, cli_cfg: str | None) -> str | None:
    """
    Priority: run-local TOML > CLI --config > project-level TOML > None
    """
    run_local = os.path.join(run_dir, "electrofit.toml")
    if os.path.isfile(run_local):
        return run_local
    if cli_cfg:
        return cli_cfg
    project_toml = os.path.join(project_root, "electrofit.toml")
    if os.path.isfile(project_toml):
        return project_toml
    return None


def _run_one_dir(run_dir: str, project_root: str, cfg_path: str | None):
    """
    Execute the Python pipeline in a single run directory.
    """
    os.environ["ELECTROFIT_PROJECT_PATH"] = project_root
    if cfg_path:
        os.environ["ELECTROFIT_CONFIG_PATH"] = cfg_path
    else:
        os.environ.pop("ELECTROFIT_CONFIG_PATH", None)

    # run_process_initial_structure expects cwd == run_dir
    with _pushd(run_dir):
        run_pis_main()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--config", help="Optional electrofit.toml to use")
    args = ap.parse_args()

    project_root = os.path.abspath(args.project)

    # 1) Deterministic discovery under project/process/**/<TARGET_RUN_DIRS>
    discovered = list(_iter_run_dirs(project_root))

    if discovered:
        for run_dir in discovered:
            cfg_path = _pick_config(run_dir, project_root, args.config)
            _run_one_dir(run_dir, project_root, cfg_path)
        return

    # 2) Fallback: treat the current directory as a single run (keeps tests simple)
    run_dir = os.getcwd()
    cfg_path = _pick_config(run_dir, project_root, args.config)
    _run_one_dir(run_dir, project_root, cfg_path)


if __name__ == "__main__":
    main()