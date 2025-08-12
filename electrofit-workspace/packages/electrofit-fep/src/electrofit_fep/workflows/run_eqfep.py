#!/usr/bin/env python3
"""
run_eqFEP.py

Driver script for running the eqFEP workflow for multi-lambda free energy calculations in water or vacuum.
Optionally, it sets up a scratch directory for the run if the --scratch flag is provided.
"""

import argparse
import logging
import os

import numpy as np
from electrofit.logging import setup_logging
from electrofit.scratch.manager import (
    finalize_scratch_directory,
    setup_scratch_directory,
)
from electrofit_fep.core.eqfep import eqFEP, generate_cos_lambda_scheme


def main():
    parser = argparse.ArgumentParser(
        description="Run eqFEP in cascade or parallel mode, processing one or more states and edges."
    )
    parser.add_argument(
        "--cascade",
        action="store_true",
        default=False,
        help="Run transitions in cascade mode (each window uses previous confout).",
    )
    parser.add_argument(
        "-s",
        "--state",
        nargs="+",
        type=str,
        default=["stateA"],
        help='One or more states to run, e.g. "-s stateA" or "-s stateA stateB".',
    )
    parser.add_argument(
        "-w",
        "--workPath",
        type=str,
        default="workpath",
        help="Main working directory for eqFEP output.",
    )
    parser.add_argument(
        "-m", "--mdpPath", type=str, default="mdp", help="Path to MDP files."
    )
    parser.add_argument(
        "-l",
        "--ligandPath",
        type=str,
        default="input/molecules",
        help="Path to your ligand directories (e.g., IP_011111).",
    )
    parser.add_argument(
        "-e",
        "--edge",
        action="append",
        nargs=2,
        metavar=("LIG_A", "LIG_B"),
        help="Specify a ligand pair edge as two ligand IDs. Can be used multiple times.",
    )
    parser.add_argument(
        "-b",
        "--windowsBatch",
        type=int,
        default=4,
        help="Number of lambda windows to run concurrently (if not cascade).",
    )
    parser.add_argument(
        "-r", "--replicas", type=int, default=1, help="Number of replicate simulations."
    )
    parser.add_argument(
        "-c",
        "--total_cpus",
        type=int,
        default=40,
        help="Total CPU cores available (used to compute ntomp).",
    )
    parser.add_argument(
        "-n",
        "--n_lambdas",
        type=int,
        default=10,
        help="Number of equally spaced lambda windows (default = 10).",
    )
    # New argument: boxd (box distance from solute)
    parser.add_argument(
        "--boxd",
        type=float,
        default=1.2,
        help="Box distance from solute in nm (default is 1.2 nm).",
    )

    # Scratch management arguments:
    parser.add_argument(
        "--scratch",
        action="store_true",
        default=False,
        help="If set, set up a scratch directory and copy input files there.",
    )
    parser.add_argument(
        "-bsd",
        "--base_scratch_dir",
        type=str,
        default="/scratch/johannal96/tmp/",
        help="Base scratch directory to use if --scratch is set.",
    )
    parser.add_argument(
        "-id",
        "--input_dir",
        type=str,
        default="input",
        help='Input directory to copy to scratch (e.g. containing the "molecules" folder).',
    )

    # Lambda spacing arguments:
    parser.add_argument(
        "--cos_lambda",
        action="store_true",
        default=False,
        help="Use a non-uniform lambda spacing based on a cosine scheme.",
    )
    parser.add_argument(
        "-exp",
        "--exponent",
        type=float,
        default=1.5,
        help="Exponent for the cosine lambda scheme (only used if --cos_lambda is set).",
    )

    # Vacuum mode argument:
    parser.add_argument(
        "--vacuum",
        action="store_true",
        default=False,
        help="Run in vacuum mode (no solvation or ions).",
    )

    args = parser.parse_args()

    # Set up logging
    fullpath = os.getcwd()
    log_file_path = os.path.join(fullpath, "eqfep_process.log")
    setup_logging(log_file_path)
    logging.info(f"Logging initialized. Log file: {log_file_path}")

    # Parse edges
    if args.edge is None:
        edges = [("011111", "101111")]
    else:
        edges = args.edge

    # States list
    states_list = args.state

    # Configure lambda schedule
    if args.cos_lambda:
        lambda_list = list(
            generate_cos_lambda_scheme(n=args.n_lambdas, exponent=args.exponent)
        )
        logging.info(
            f"Using cosine lambda scheme with exponent {args.exponent}: {lambda_list}"
        )
    else:
        lambda_list = list(np.linspace(0, 1, args.n_lambdas))
        logging.info(f"Using equidistant lambda scheme: {lambda_list}")

    # Set up scratch directory if requested.
    if args.scratch:
        input_files = [args.input_dir]
        scratch_dir, original_dir = setup_scratch_directory(
            input_files, args.base_scratch_dir
        )
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
        states=states_list,
        boxd=args.boxd,
        boxshape="cubic",
    )

    logging.info("Starting eqFEP pipeline...")

    # Main pipeline steps
    fe.prepareFreeEnergyDir()
    fe.atom_mapping(bVerbose=True)
    fe.hybrid_structure_topology(bVerbose=True)
    fe.assemble_systems()

    # Conditionally skip solvation/ion steps if vacuum.
    if args.vacuum:
        # For vacuum simulation, disable water and ion steps.
        fe.boxWaterIons(bIon=False, bSolvate=False)
    else:
        fe.boxWaterIons(bIon=True, bSolvate=True)

    # Prepare and run the energy minimization (EM), NVT, and NPT steps.
    fe.prepare_simulation(simType="em")
    fe.run_simulation_locally(simType="em")
    fe.prepare_simulation(simType="nvt")
    fe.run_simulation_locally(simType="nvt")
    fe.prepare_simulation(simType="npt")
    fe.run_simulation_locally(simType="npt")

    # Prepare and run the alchemical transitions
    fe.prepare_transitions(bGenTpr=True)
    fe.run_simulation_locally(simType="transitions", bVerbose=True)

    # Run BAR analysis on the results.
    fe.run_analysis(bVerbose=True, start_time=500)
    logging.info("eqFEP run complete.")

    # Finalize scratch directory if used.
    if args.scratch:
        finalize_scratch_directory(original_dir, scratch_dir, [args.input_dir])
        logging.info("Scratch directory finalized.")


if __name__ == "__main__":
    main()
