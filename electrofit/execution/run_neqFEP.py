#!/usr/bin/env python3
import logging
import os
import subprocess
import shutil

import matplotlib.pyplot as plt
import pandas as pd

import os
import sys

def find_project_root(current_dir, project_name="electrofit"):
    root = None
    while True:
        parent_dir = os.path.dirname(current_dir)
        if os.path.basename(current_dir) == project_name:
            root = current_dir  # Set root to the current project_name directory
        if parent_dir == current_dir:
            # We've reached the filesystem root
            if root is None:
                raise FileNotFoundError(f"Project root directory '{project_name}' not found.")
            return root  # Return the outermost match found
        current_dir = parent_dir


script_dir = os.path.dirname(os.path.abspath(__file__))
project_path = find_project_root(current_dir=script_dir)

sys.path.append(project_path)

# Import helper functions from your electrofit modules
from electrofit.commands.run_commands import run_command
from electrofit.helper.file_manipulation import (
    include_ff,
    include_ions,
    include_tip3p,
    remove_defaults_section_lines,
    replace_posres_in_file,
    strip_extension,
)
from electrofit.helper.set_logging import setup_logging
from electrofit.helper.setup_finalize_scratch import (
    setup_scratch_directory,
    finalize_scratch_directory,
)

base_scratch_dir = "/scratch/johannal96/tmp/"

# Set up logging
fullpath = os.getcwd()
log_file_path = os.path.join(fullpath, "neqfep_process.log")
setup_logging(log_file_path)
logging.info(f"Logging initialized. Log file: {log_file_path}")
input_dir = "input"

# Set up scratch directory and copy input files
input_files = [input_dir]
scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch_dir)
os.chdir(scratch_dir)
logging.info(f"Changed working directory to scratch directory: {scratch_dir}")

try:
    run_command(
        'python /home/johannal96/MA/pmx/arthur/scripts/run_neqFEP.py',
        cwd=scratch_dir,
    )

except Exception as err:
    logging.error(f"Error occurred during non-equilibrium FEP: {err}")
    logging.info("Finalizing scratch directory due to error.")
    finalize_scratch_directory(original_dir, scratch_dir, input_files)
    raise  # Optionally re-raise the error for further handling

# Finalize the scratch directory (if we reach this point without errors)
finalize_scratch_directory(original_dir, scratch_dir, input_files)
