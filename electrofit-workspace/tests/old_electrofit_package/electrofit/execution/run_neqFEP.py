#!/usr/bin/env python3
"""
run_neqFEP.py

Driver script for running the neqFEP workflow for non-equilibrium alchemical free energy calculations of ligand A to B in water..

Usage examples:
  python run_neqFEP.py -s stateA stateB -e 011111 101111 -r 3 -m input/mdp -l input/molecules -w workpath --scratch --base_scratch_dir /scratch/johannal96/tmp/ --input_dir input

"""

import argparse
import logging
import os
import sys

class StreamToLogger(object):
    """
    Redirects writes to a logger instance.
    """
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''
    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())
    def flush(self):
        pass

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

# Import helper functions and class 
from electrofit.helper.neqFEP import neqFEP
from electrofit.helper.set_logging import setup_logging
from electrofit.helper.setup_finalize_scratch import setup_scratch_directory, finalize_scratch_directory

# Set up logging as before.
fullpath = os.getcwd()
log_file_path = os.path.join(fullpath, "neqfep_process.log")
setup_logging(log_file_path)
logging.info(f"Logging initialized. Log file: {log_file_path}")

# Redirect stdout and stderr:
sys.stdout = StreamToLogger(logging.getLogger('STDOUT'), logging.INFO)
sys.stderr = StreamToLogger(logging.getLogger('STDERR'), logging.ERROR)

def main():
    parser = argparse.ArgumentParser(
        description="Run neqFEP workflow for non-equilibrium alchemical free energy calculations of ligand A to B in water."
    )
    # State arguments
    parser.add_argument('-s', '--state', nargs='+', type=str, default=['stateA'],
                        help='One or more states to run, e.g. "-s stateA" or "-s stateA stateB".')
    # Working directories and input
    parser.add_argument('-w', '--workPath', type=str, default='workpath',
                        help='Main working directory for eqFEP output.')
    parser.add_argument('-m', '--mdpPath', type=str, default='input/mdp',
                        help='Path to MDP files.')
    parser.add_argument('-l', '--ligandPath', type=str, default='input/molecules',
                        help='Path to ligand directories (e.g. directories named like IP_011111).')
    parser.add_argument('-e', '--edge', action='append', nargs=2, metavar=('LIG_A', 'LIG_B'),
                        help='Specify a ligand pair edge as two ligand IDs. Can be used multiple times.')
    parser.add_argument('-r', '--replicas', type=int, default=3,
                        help='Number of replicate simulations.')
    # Force field and simulation parameters
    parser.add_argument('--ff', type=str, default='amber14sb.ff',
                        help='Force field directory.')
    parser.add_argument('--boxshape', type=str, default='dodecahedron',
                        help='Box shape.')
    parser.add_argument('--boxd', type=float, default=1.2,
                        help='Box distance (in nm).')
    parser.add_argument('--water', type=str, default='tip3p',
                        help='Water model.')
    parser.add_argument('--conc', type=float, default=0.15,
                        help='Ion concentration.')
    parser.add_argument('--pname', type=str, default='NA',
                        help='Positive ion name.')
    parser.add_argument('--nname', type=str, default='CL',
                        help='Negative ion name.')
    parser.add_argument('--start_time', type=str, default=2250,
                        help='Negative ion name.')
    # Scratch management
    parser.add_argument('--scratch', action='store_true', default=False,
                        help='If set, set up a scratch directory and copy input files there.')
    parser.add_argument('--base_scratch_dir', type=str, default='/scratch/johannal96/tmp/',
                        help='Base scratch directory to use if --scratch is set.')
    parser.add_argument('--input_dir', type=str, default='input',
                        help='Input directory to copy to scratch.')
    args = parser.parse_args()

    # Process edges; if none provided, use a default.
    if args.edge is None:
        edges = [('molA', 'molB')]
    else:
        edges = [tuple(e) for e in args.edge]
    
    # Optional: Set up scratch directory if requested.
    if args.scratch:
        input_files = [args.input_dir]
        base_dir = args.base_scratch_dir
        scratch_dir, original_dir = setup_scratch_directory(input_files, base_dir)
        os.chdir(scratch_dir)
        logging.info(f"Changed working directory to scratch directory: {scratch_dir}")
    else:
        scratch_dir = None
        original_dir = None

    # Instantiate eqFEP object
    fe = neqFEP(
        workPath=args.workPath,
        mdpPath=args.mdpPath,
        ligandPath=args.ligandPath,
        edges=edges,
        replicas=args.replicas,
        ff=args.ff,
        boxshape=args.boxshape,
        boxd=args.boxd,
        water=args.water,
        conc=args.conc,
        pname=args.pname,
        nname=args.nname,
        states=args.state
    )


    # Run the workflow steps
    fe.prepareFreeEnergyDir()
    fe.atom_mapping(bVerbose=True)
    fe.hybrid_structure_topology(bVerbose=True)
    fe.assemble_systems()
    fe.boxWaterIons()
    fe.prepare_simulation(simType='em')
    fe.run_simulation_locally(simType='em', bVerbose=False)
    # Run NVT equilibration (input from em)
    fe.prepare_simulation(simType='nvt')
    fe.run_simulation_locally(simType='nvt', bVerbose=False)
    # Run NPT equilibration (input from nvt) - run for longer time untill fully equilibrated
    fe.prepare_simulation(simType='npt')
    fe.run_simulation_locally(simType='npt', bVerbose=False)
    if args.start_time:
        fe.prepare_transitions(startTime=args.start_time)
    else:
        fe.prepare_transitions()
    fe.run_simulation_locally(simType='transitions', bVerbose=False)
    fe.run_analysis(bVerbose=True)
    fe.analysis_summary()

    logging.info("neqFEP run complete.")

    # If scratch was used, finalize the scratch directory.
    if args.scratch:
        finalize_scratch_directory(original_dir, scratch_dir, input_files)
        logging.info("Scratch directory finalized.")

if __name__ == '__main__':
    main()