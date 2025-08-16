"""Step0: Copy input molecule folders (data/input/*) into process/* as working copies.

Usage: electrofit step0 --project <path> [--log-console]
"""
from pathlib import Path
import os
import argparse
import logging

from electrofit.io.files import copy_and_rename_folders
from electrofit.logging import setup_logging, log_run_header


def main():  # pragma: no cover (CLI wrapper)
    ap = argparse.ArgumentParser(description="Step0: initialize process/ workspace by copying data/input")
    ap.add_argument("--project", default=os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd()))
    ap.add_argument("--log-console", action="store_true", help="Also echo logs to console")
    args = ap.parse_args()
    project_path = Path(args.project).resolve()
    src = project_path / "data" / "input"
    dst = project_path / "process"
    dst.mkdir(exist_ok=True)
    setup_logging(str(project_path / "step0.log"), also_console=args.log_console)
    log_run_header("step0")
    logging.info(f"[step0] Source: {src}")
    logging.info(f"[step0] Destination: {dst}")
    if not src.exists():
        raise FileNotFoundError(f"Source directory '{src}' does not exist.")
    copy_and_rename_folders(source=str(src), destination=str(dst))
    logging.info("[step0] Copy complete.")


if __name__ == "__main__":
    main()
