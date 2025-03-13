#!/usr/bin/env python3
"""
run_eqFEP.py

Driver script for running the eqFEP workflow for multi-lambda free energy calculations in water.
Optionally, it sets up a scratch directory for the run if the --scratch flag is provided.

Usage examples:
  # Cascade mode for two states, with scratch management:
  python run_eqFEP.py --cascade -s stateA stateB -e 011111 101111 -c 40 -n 20 -r 1 -m input/mdp --scratch -bsd /scratch/user/tmp/ --id input

  # Parallel mode for stateA only (no scratch):
  python run_eqFEP.py -s stateA -e 011111 101111 -c 40 -n 10 -r 1
"""

import argparse
import logging
import sys
import os
import numpy as np
import glob

# Helper: find project root (if needed)
def find_project_root(current_dir, project_name="electrofit"):
    root = None
    while True:
        parent_dir = os.path.dirname(current_dir)
        if os.path.basename(current_dir) == project_name:
            root = current_dir
        if parent_dir == current_dir:
            if root is None:
                raise FileNotFoundError(f"Project root directory '{project_name}' not found.")
            return root
        current_dir = parent_dir

# Set up project root and include it in PYTHONPATH
script_dir = os.path.dirname(os.path.abspath(__file__))
project_path = find_project_root(script_dir)
sys.path.append(project_path)

from electrofit.helper.eqFEP import eqFEP
from electrofit.helper.set_logging import setup_logging
from electrofit.helper.setup_finalize_scratch import setup_scratch_directory, finalize_scratch_directory

def main():
    parser = argparse.ArgumentParser(
        description="Run eqFEP in cascade or parallel mode, processing one or more states and edges."
    )
    parser.add_argument('--cascade', action='store_true', default=False,
                        help='Run transitions in cascade mode (each window uses previous confout).')
    parser.add_argument('-s', '--state', nargs='+', type=str, default=['stateA'],
                        help='One or more states to run. E.g.: "-s stateA" or "-s stateA stateB".')
    parser.add_argument('-w', '--workPath', type=str, default='workpath',
                        help='Main working directory for eqFEP output.')
    parser.add_argument('-m', '--mdpPath', type=str, default='mdp',
                        help='Path to MDP files.')
    parser.add_argument('-l', '--ligandPath', type=str, default='input/molecules',
                        help='Path to your ligand directories (e.g. IP_011111).')
    parser.add_argument('-e', '--edge', action='append', nargs=2, metavar=('LIG_A', 'LIG_B'),
                        help='Specify a ligand pair edge as two ligand IDs. Can be used multiple times.')
    parser.add_argument('-b', '--windowsBatch', type=int, default=4,
                        help='Number of lambda windows to run concurrently (if not cascade).')
    parser.add_argument('-r', '--replicas', type=int, default=1,
                        help='Number of replicate simulations.')
    parser.add_argument('-c', '--total_cpus', type=int, default=40,
                        help='Total CPU cores available (used to compute ntomp).')
    parser.add_argument('-n', '--n_lambdas', type=int, default=10,
                        help='Number of equally spaced lambda windows (default is 10).')
    # Optional scratch management arguments:
    parser.add_argument('--scratch', action='store_true', default=False,
                        help='If set, set up a scratch directory and copy input files there.')
    parser.add_argument('-bsd', '--base_scratch_dir', type=str, default='/scratch/johannal96/tmp/',
                        help='Base scratch directory to use if --scratch is set.')
    parser.add_argument('-id', '--input_dir', type=str, default='input',
                        help='Input directory to copy to scratch (e.g. containing the "molecules" folder).')
    args = parser.parse_args()

    # Set up logging; create a log file in the current directory
    fullpath = os.getcwd()
    log_file_path = os.path.join(fullpath, "eqfep_process.log")
    setup_logging(log_file_path)
    logging.info(f"Logging initialized. Log file: {log_file_path}")

    # Parse edges; if none provided, use a default.
    if args.edge is None:
        edges = [('011111', '101111')]
    else:
        edges = args.edge

    # Use the provided states list.
    states_list = args.state

    # Create a lambda schedule from 0 to 1 with n_lambdas values.
    lambda_list = list(np.linspace(0, 1, args.n_lambdas))

    # Optional: Set up scratch directory if requested.
    if args.scratch:
        input_files = [args.input_dir]
        scratch_dir, original_dir = setup_scratch_directory(input_files, args.base_scratch_dir)
        os.chdir(scratch_dir)
        logging.info(f"Changed working directory to scratch directory: {scratch_dir}")
    else:
        scratch_dir = None
        original_dir = None

    # Instantiate eqFEP with the parsed parameters.
    fe = eqFEP(
        workPath=args.workPath,
        mdpPath=args.mdpPath,
        ligandPath=args.ligandPath,
        edges=edges,
        replicas=args.replicas,
        windowsBatch=args.windowsBatch,
        total_cpus=args.total_cpus,
        lambdaStates=lambda_list,
        cascade=args.cascade,
        states=states_list
    )

    logging.info("Starting eqFEP pipeline...")

    # Run the workflow.
    fe.prepareFreeEnergyDir()
    fe.atom_mapping(bVerbose=True)
    fe.hybrid_structure_topology(bVerbose=True)
    fe.assemble_systems()
    fe.boxWaterIons()
    fe.prepare_simulation(simType='em')
    fe.run_simulation_locally(simType='em')
    fe.prepare_simulation(simType='nvt')
    fe.run_simulation_locally(simType='nvt')
    fe.prepare_simulation(simType='npt')
    fe.run_simulation_locally(simType='npt')
    fe.prepare_transitions(bGenTpr=True)
    fe.run_simulation_locally(simType='transitions', bVerbose=True)
    fe.run_analysis(bVerbose=True, start_time=100)

    logging.info("eqFEP run complete.")

    # If scratch was used, finalize the scratch directory.
    if args.scratch:
        finalize_scratch_directory(original_dir, scratch_dir, input_files)
        logging.info("Scratch directory finalized.")

if __name__ == '__main__':
    main()