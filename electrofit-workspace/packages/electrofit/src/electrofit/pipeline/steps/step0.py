"""Pipeline Step 0: Initialize process workspace.

Copies input molecule subdirectories from ``data/input/*`` into ``process/*`` as
working copies. Pure orchestration: delegates filesystem copying to
``electrofit.io.files.copy_and_rename_folders`` and logging to infra layer.

CLI equivalent (legacy): ``electrofit step0 --project <path>``
"""
from __future__ import annotations

from pathlib import Path
import os
import argparse
import logging

from electrofit.io.files import copy_and_rename_folders
from electrofit.infra.logging import setup_logging, log_run_header

__all__ = ["main"]

def main():  # pragma: no cover (CLI orchestration)
    parser = argparse.ArgumentParser(description="Step0: initialize process/ workspace by copying data/input")
    parser.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    parser.add_argument("--log-console", action="store_true", help="Also echo logs to console")
    args = parser.parse_args()

    project_path = Path(args.project).resolve()
    src = project_path / "data" / "input"
    dst = project_path / "process"
    dst.mkdir(exist_ok=True)
    setup_logging(str(project_path / "step.log"), also_console=args.log_console)
    log_run_header("step0")
    logging.info("[step0] Source: %s", src)
    logging.info("[step0] Destination: %s", dst)
    if not src.exists():
        raise FileNotFoundError(f"Source directory '{src}' does not exist.")
    copy_and_rename_folders(source=str(src), destination=str(dst))
    logging.info("[step0] Copy complete.")

if __name__ == "__main__":  # pragma: no cover
    main()
