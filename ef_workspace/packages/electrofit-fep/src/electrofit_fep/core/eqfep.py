#!/usr/bin/env python3
import concurrent.futures
import glob
import logging
import os

import numpy as np
import pandas as pd
from electrofit.cli.run_commands import run_command
from pmx import gmx, ligand_alchemy
from pmx.utils import create_folder


def generate_cos_lambda_scheme(n, exponent=1.5):
    """
    Generate a lambda scheme with non-uniform spacing, denser at both lambda=0 and lambda=1.

    Parameters
    ----------
    n : int
        Total number of lambda windows.
    exponent : float, optional
        Exponent for biasing the distribution.
        - If exponent = 1.0, the scheme is the standard cosine scheme.
        - If exponent > 1.0, points are denser at both ends.
        Default is 2.0.

    Returns
    -------
    lambdas : np.ndarray
        Array of lambda values in [0, 1] with enhanced density at the endpoints.
    """
    # Step 1: Generate standard cosine spacing:
    x = np.linspace(0, np.pi, n)
    u = 0.5 * (1 - np.cos(x))

    # Step 2: Apply the double-end enrichment transform.
    # This transformation maps u in [0,1] to lambda such that:
    #   For u near 0, lambda ~ u^exponent, and for u near 1, lambda ~ 1 - (1-u)^exponent.
    lambdas = u**exponent / (u**exponent + (1 - u) ** exponent)
    return lambdas


class eqFEP:
    """
    eqFEP: A class for running multi-lambda free energy calculations (TI-FEP) in water using GROMACS, PMX,
    and a specified force field.

    Overview
    --------
    This workflow automates the setup, execution, and analysis of alchemical free-energy calculations
    for ligand transformations in an aqueous environment. It supports two main approaches:
    1. Parallel (batch) mode (default): All lambda windows start from the same equilibrated NPT snapshot.
    2. Cascade mode (cascade=True): Each lambda window depends on the final confout of the previous window.

    Key Stages:
    -----------
    1. Directory Preparation
    - Create consistent folder structures for each edge (ligand pair), state, replicate, and simulation type.

    2. Ligand Processing
    - Use PMX to map atoms between two ligands and create hybrid topologies.

    3. System Assembly in Water
    - Build a minimal top file, solvate the system, and add ions.

    4. Equilibration Steps
    - Energy minimization (EM), then NVT, then NPT simulations.

    5. Alchemical Transitions (TI-FEP)
    - Prepare separate .tpr files for each lambda window (either up front or on the fly).
    - Run them either in parallel batches (parallel approach) or sequentially (cascade approach).

    6. Analysis
    - Invoke GROMACS 'gmx bar' on the dhdl.xvg files to estimate free-energy differences.

    Main Features:
    --------------
    - Automatic directory structure creation for each replicate/state/edge.
    - Handling of GROMACS commands (grompp, mdrun) via direct function calls or run_command.
    - Flexibility in running transitions (cascade vs. parallel/batch).
    - Optional positional restraints for water and ligands.
    - A built-in basic BAR analysis routine.

    Usage Example:
    --------------
    .. code-block:: python

        from eqFEP import eqFEP

        fe = eqFEP(
            workPath="my_workpath",
            mdpPath="my_mdp_files",
            ligandPath="my_ligands",
            edges=[("011111", "101111")],  # list of ligand ID pairs
            replicas=1,
            windowsBatch=4,
            total_cpus=40,
            lambdaStates=[0.00, 0.05, ..., 1.00],
            cascade=False
        )

        # 1) Prepare directory structures
        fe.prepareFreeEnergyDir()

        # 2) Atom mapping & hybrid topology
        fe.atom_mapping()
        fe.hybrid_structure_topology()

        # 3) System setup in water
        fe.assemble_systems()
        fe.boxWaterIons()

        # 4) EM, NVT, NPT
        fe.prepare_simulation(simType='em')
        fe.run_simulation_locally(simType='em')
        fe.prepare_simulation(simType='nvt')
        fe.run_simulation_locally(simType='nvt')
        fe.prepare_simulation(simType='npt')
        fe.run_simulation_locally(simType='npt')

        # 5) Transitions (TI-FEP windows)
        fe.prepare_transitions(bGenTpr=True)
        fe.run_simulation_locally(simType='transitions')

        # 6) BAR analysis
        fe.run_analysis(bVerbose=True, start_time=100)


    Class Attributes:
    ----------------
    - workPath (str): Main working directory.
    - mdpPath (str): Directory with MDP files/templates (e.g., em_l0.mdp, NVT_l0.mdp, etc.).
    - ligandPath (str): Directory containing subfolders like 'IP_011111' for each ligand.
    - edges (list or dict): Ligand pairs to transform. e.g., [("011111", "101111")].
    - replicas (int): Number of replicate simulations per edge and state.
    - simTypes (list): Default ["em", "nvt", "npt", "transitions"].
    - states (list): Default ["stateA", "stateB"].
    - boxshape (str): GROMACS box type (e.g., "dodecahedron").
    - boxd (float): Box distance from solute (nm).
    - water (str): Water model (e.g., "tip3p").
    - conc (float): Ion concentration (e.g., 0.15 M).
    - pname (str): Positive ion name (e.g., "NA").
    - nname (str): Negative ion name (e.g., "CL").
    - windowsBatch (int): Number of lambda windows to run concurrently in parallel mode.
    - ntmpi (int): MPI ranks for each job (usually 1).
    - total_cpus (int): Total CPU cores. Used to set ntomp in parallel mode.
    - lambdaStates (list): The lambda schedule for alchemical transformations.
    - cascade (bool): If True, run transitions in cascade mode (each window depends on the previous).

    Public Methods:
    --------------
    - prepareFreeEnergyDir(): Builds directory structures and logs initial info.
    - atom_mapping(), hybrid_structure_topology(): PMX-based ligand processing.
    - assemble_systems(), boxWaterIons(): Prepare a water box, solvate, add ions.
    - prepare_simulation(simType): Create TPR files for EM/NVT/NPT.
    - run_simulation_locally(simType): Run GROMACS locally for EM/NVT/NPT or transitions (cascade or parallel).
    - prepare_transitions(bGenTpr=True): Set up TPR files for each lambda window (parallel approach).
    - run_analysis(bVerbose=False, start_time=200): Basic BAR analysis with GROMACS 'gmx bar'.

    Notes:
    ------
    - For advanced users, you can override states, edges, or pass a custom lambda schedule.
    - The script logs to Python's logging, so set your log level or handlers as needed.
    - Adjust your HPC submission approach if you want to run states/edges concurrently.

    Dependencies:
    -------------
    - Python 3.7+
    - GROMACS 2022+ (gmx in PATH)
    - pmx (for ligand mapping, hybrid generation)
    - numpy, pandas, concurrent.futures, etc.

    Author:
    -------
    JALK
    """

    def __init__(self, **kwargs):
        # Initialize PMX gmx path
        gmx.set_gmxlib()

        # DataFrame placeholders for results
        self.resultsAll = pd.DataFrame()
        self.resultsSummary = pd.DataFrame()

        # Default directories / parameters
        self.workPath = "./"
        self.mdpPath = "./mdp"
        self.ligandPath = None
        self.edges = {}

        self.replicas = 3
        self.simTypes = ["em", "nvt", "npt", "transitions"]
        self.states = ["stateA", "stateB"]
        self.thermCycleBranches = ["water"]

        # Force field and box settings
        self.ff = "amber14sb.ff"
        self.boxshape = "dodecahedron"
        self.boxd = 1.2
        self.water = "tip3p"
        self.conc = 0.15
        self.pname = "NA"
        self.nname = "CL"

        # Parallel transitions settings (by default)
        self.windowsBatch = 4
        self.ntmpi = 1
        self.total_cpus = 20
        # Example default lambda schedule
        self.lambdaStates = [
            0.00,
            0.05,
            0.10,
            0.20,
            0.30,
            0.40,
            0.50,
            0.60,
            0.70,
            0.80,
            0.90,
            0.95,
            1.00,
        ]
        # Cascade mode flag
        self.cascade = False

        # Allow user overrides
        for key, val in kwargs.items():
            setattr(self, key, val)

        # Recalculate ntomp (OpenMP threads) from total_cpus & windowsBatch
        if not self.cascade:
            self.ntomp = self.total_cpus // self.windowsBatch
        else:
            self.ntomp = self.total_cpus

    # -------------------------------------------------------
    # 1) Directory Setup
    # -------------------------------------------------------
    def prepareFreeEnergyDir(self):
        self.workPath = os.path.abspath(self.workPath)
        self.mdpPath = os.path.abspath(self.mdpPath)
        if self.ligandPath is None:
            raise ValueError("ligandPath must be defined.")
        self.ligandPath = os.path.abspath(self.ligandPath)

        create_folder(self.workPath)
        self._read_ligands()
        self._read_edges()
        self._create_folder_structure()

        logging.info("Summary of the setup:")
        logging.info(f"   workPath:       {self.workPath}")
        logging.info(f"   mdpPath:        {self.mdpPath}")
        logging.info(f"   ligandPath:     {self.ligandPath}")
        logging.info(f"   replicas:       {self.replicas}")
        logging.info(f"   windowsBatch:   {self.windowsBatch}")
        logging.info(f"   cascade mode?:  {self.cascade}")
        logging.info("   edges:")
        for e in self.edges:
            logging.info(f"        {e}")

        self._print_folder_structure()
        logging.info("Finished preparing free energy directory.\n")

    def _read_ligands(self):
        self.ligands = {}
        lig_folders = glob.glob(f"{self.ligandPath}/*")
        for folder in lig_folders:
            basename = os.path.basename(folder)
            if basename.startswith("IP_"):
                lname = basename[3:]
            else:
                lname = basename
            self.ligands[lname] = os.path.abspath(folder)

    def _read_edges(self):
        # If edges is a list of tuples, convert to dict. If edges is a file path, read them, etc.
        if isinstance(self.edges, list):
            edge_dict = {}
            for e in self.edges:
                key = f"edge_{e[0]}_{e[1]}"
                edge_dict[key] = e
            self.edges = edge_dict
        elif isinstance(self.edges, str) and os.path.isfile(self.edges):
            # parse file
            self.edges = "Edges read from file"

    def _create_folder_structure(self):
        for edge in self.edges:
            edgepath = os.path.join(self.workPath, edge)
            create_folder(edgepath)
            create_folder(os.path.join(edgepath, "hybridStrTop"))
            for wp in self.thermCycleBranches:
                wppath = os.path.join(edgepath, wp)
                create_folder(wppath)
                for state in self.states:
                    statepath = os.path.join(wppath, state)
                    create_folder(statepath)
                    for r in range(1, self.replicas + 1):
                        runpath = os.path.join(statepath, f"run{r}")
                        create_folder(runpath)
                        for sim in ["em", "nvt", "npt"]:
                            create_folder(os.path.join(runpath, sim))
                        transBase = os.path.join(runpath, "transitions")
                        create_folder(transBase)
                        for i, _ in enumerate(self.lambdaStates):
                            create_folder(os.path.join(transBase, f"lambda_{i}"))

    def _print_folder_structure(self):
        logging.info("Partial folder structure:")
        logging.info(f"{self.workPath}/")
        logging.info(" |-- edge_X_Y")
        logging.info(" |    |-- hybridStrTop")
        logging.info(" |    |-- water")
        logging.info(" |         |-- stateA")
        logging.info(" |         |    |-- run1/2/3")
        logging.info(" |         |         |-- em/nvt/npt/transitions/")
        logging.info(" |         |              |-- lambda_0..N")
        logging.info(" |         |-- stateB")
        logging.info(" |              |-- run1/2/3")
        logging.info(" |                   |-- em/nvt/npt/transitions/")
        logging.info(" |                        |-- lambda_0..N")

    # -------------------------------------------------------
    # 2) Atom Mapping & Hybrid Structure/Topology
    # -------------------------------------------------------
    def atom_mapping(self, edges=None, bVerbose=False):
        logging.info("Performing atom mapping")
        if edges is None:
            edges = self.edges
        for edge in edges:
            logging.info(f"Atom mapping for edge: {edge}")
            ligA, ligB = self.edges[edge]
            ligApath = os.path.join(self.ligandPath, f"IP_{ligA}")
            ligBpath = os.path.join(self.ligandPath, f"IP_{ligB}")
            ligAname = f"IP_{ligA}"
            ligBname = f"IP_{ligB}"
            outpath = os.path.join(self.workPath, edge, "hybridStrTop")
            i1 = os.path.join(ligApath, f"{ligAname}.pdb")
            i2 = os.path.join(ligBpath, f"{ligBname}.pdb")
            o1 = os.path.join(outpath, "pairs1.dat")
            o2 = os.path.join(outpath, "pairs2.dat")
            opdb1 = os.path.join(outpath, "out_pdb1.pdb")
            opdb2 = os.path.join(outpath, "out_pdb2.pdb")
            opdbm1 = os.path.join(outpath, "out_pdbm1.pdb")
            opdbm2 = os.path.join(outpath, "out_pdbm2.pdb")
            score = os.path.join(outpath, "score.dat")
            logf = os.path.join(outpath, "mapping.log")
            cmd = (
                "pmx atomMapping "
                f"-i1 {i1} -i2 {i2} "
                f"-o1 {o1} -o2 {o2} "
                f"-opdb1 {opdb1} -opdb2 {opdb2} "
                f"-opdbm1 {opdbm1} -opdbm2 {opdbm2} "
                f"-score {score} -log {logf} "
                "--H2Hpolar --no-mcs"
            )
            out = run_command(cmd)
            if bVerbose and out:
                logging.info(out)
        logging.info("Finished atom mapping.")

    def hybrid_structure_topology(self, edges=None, bVerbose=False):
        logging.info("Creating hybrid structure/topology")
        if edges is None:
            edges = self.edges
        for edge in edges:
            logging.info(f"Hybrid structure for edge: {edge}")
            ligA, ligB = self.edges[edge]
            ligApath = os.path.join(self.ligandPath, f"IP_{ligA}")
            ligBpath = os.path.join(self.ligandPath, f"IP_{ligB}")
            ligAname = f"IP_{ligA}"
            ligBname = f"IP_{ligB}"
            outpath = os.path.join(self.workPath, edge, "hybridStrTop")
            i1 = os.path.join(ligApath, f"{ligAname}.pdb")
            i2 = os.path.join(ligBpath, f"{ligBname}.pdb")
            itp1 = os.path.join(ligApath, f"{ligAname}_GMX.itp")
            itp2 = os.path.join(ligBpath, f"{ligBname}_GMX.itp")
            pairs = os.path.join(outpath, "pairs1.dat")
            oA = os.path.join(outpath, "mergedA.pdb")
            oB = os.path.join(outpath, "mergedB.pdb")
            oitp = os.path.join(outpath, "merged.itp")
            offitp = os.path.join(outpath, "ffmerged.itp")
            logf = os.path.join(outpath, "hybrid.log")
            cmd = (
                "pmx ligandHybrid "
                f"-i1 {i1} -i2 {i2} "
                f"-itp1 {itp1} -itp2 {itp2} "
                f"-pairs {pairs} "
                f"-oA {oA} -oB {oB} "
                f"-oitp {oitp} -offitp {offitp} "
                f"-log {logf}"
            )
            out = run_command(cmd)
            if bVerbose and out:
                logging.info(out)
        logging.info("Finished creating hybrid structures.")

    # -------------------------------------------------------
    # 3) System Assembly (Water Only)
    # -------------------------------------------------------
    def assemble_systems(self, edges=None):
        logging.info("Assembling the systems (water only).")
        if edges is None:
            edges = self.edges
        for edge in edges:
            logging.info(f"Assembling system for edge: {edge}")
            ligA, ligB = self.edges[edge]
            hybridStrTopPath = os.path.join(self.workPath, edge, "hybridStrTop")
            outLigPath = os.path.join(self.workPath, edge, "water")

            # 1) Create init.pdb from mergedA
            self._make_clean_pdb(
                os.path.join(hybridStrTopPath, "mergedA.pdb"),
                os.path.join(outLigPath, "init.pdb"),
            )

            # 2) Merge force field files
            ffitpOut = os.path.join(hybridStrTopPath, "ffmerged.itp")
            ffitpIn1 = os.path.join(self.ligandPath, f"IP_{ligA}", f"ff_{ligA}.itp")
            ffitpIn2 = os.path.join(self.ligandPath, f"IP_{ligB}", f"ff_{ligB}.itp")
            ffitpIn3 = os.path.join(hybridStrTopPath, "ffmerged.itp")
            ligand_alchemy._merge_FF_files(
                ffitpOut, ffsIn=[ffitpIn1, ffitpIn2, ffitpIn3]
            )

            # 3) Create topol.top
            ligTopFname = os.path.join(outLigPath, "topol.top")
            ligFFitp = os.path.join(hybridStrTopPath, "ffmerged.itp")
            ligItp = os.path.join(hybridStrTopPath, "merged.itp")
            ligposresItp = os.path.join(self.ligandPath, "posres", "posres.itp")
            itps = {"FFitp": ligFFitp, "ligItp": ligItp, "posresItp": ligposresItp}
            self._create_top(ligTopFname, itps, f"IP_{ligA}")
        logging.info("System assembly complete.")

    def _make_clean_pdb(self, fnameIn, fnameOut):
        with open(fnameIn, "r") as f:
            lines = f.readlines()
        keep = [line for line in lines if line.startswith("ATOM") or line.startswith("HETATM")]
        with open(fnameOut, "w") as w:
            w.writelines(keep)

    def _create_top(self, fname, itp, ligAname, systemName="ligand in water"):
        with open(fname, "w") as fp:
            fp.write(f'#include "{self.ff}/forcefield.itp"\n')
            for key in itp:
                if key != "posresItp":
                    fp.write(f'#include "{itp[key]}"\n')
            fp.write("#ifdef POSRES\n")
            fp.write(f'#include "{itp["posresItp"]}"\n')
            fp.write("#endif\n")
            fp.write(f'#include "{self.ff}/{self.water}.itp"\n')
            fp.write("#ifdef POSRES_WATER\n")
            fp.write("; Position restraint for each water oxygen\n")
            fp.write("[ position_restraints ]\n")
            fp.write(";  i   funct    fcx    fcy    fcz\n")
            fp.write("   1    1      1000   1000   1000\n")
            fp.write("#endif\n\n")
            fp.write(f'#include "{self.ff}/ions.itp"\n\n')
            fp.write("[ system ]\n")
            fp.write(f"{systemName}\n\n")
            fp.write("[ molecules ]\n")
            fp.write(f"{ligAname} 1\n")

    # -------------------------------------------------------
    # 4) Box, Solvate, Ions
    # -------------------------------------------------------
    def boxWaterIons(self, edges=None, bBox=True, bSolvate=True, bIon=True):
        logging.info("Performing box creation, solvation, and ion addition.")
        if edges is None:
            edges = self.edges
        for edge in edges:
            logging.info(f"Edge: {edge}")
            outLigPath = os.path.join(self.workPath, edge, "water")
            if bBox:
                inStr = os.path.join(outLigPath, "init.pdb")
                outStr = os.path.join(outLigPath, "box.pdb")
                gmx.editconf(inStr, o=outStr, bt=self.boxshape, d=self.boxd)
            if bSolvate:
                inBox = os.path.join(outLigPath, "box.pdb")
                outWater = os.path.join(outLigPath, "water.pdb")
                top = os.path.join(outLigPath, "topol.top")
                gmx.solvate(inBox, cs="spc216.gro", p=top, o=outWater)
            if bIon:
                inWater = os.path.join(outLigPath, "water.pdb")
                outIons = os.path.join(outLigPath, "ions.pdb")
                mdp = os.path.join(self.mdpPath, "em_l0.mdp")
                tpr = os.path.join(outLigPath, "tpr.tpr")
                top = os.path.join(outLigPath, "topol.top")
                mdout = os.path.join(outLigPath, "mdout.mdp")
                gmx.grompp(
                    f=mdp,
                    c=inWater,
                    p=top,
                    o=tpr,
                    maxwarn=4,
                    other_flags=f" -po {mdout}",
                )
                gmx.genion(
                    s=tpr,
                    p=top,
                    o=outIons,
                    conc=self.conc,
                    neutral=True,
                    other_flags=f" -pname {self.pname} -nname {self.nname}",
                )
        logging.info("Done with box/solvate/ions.")

    # -------------------------------------------------------
    # 5) Prepare & Run Simulations (EM/NVT/NPT)
    # -------------------------------------------------------
    def prepare_simulation(self, edges=None, simType="em"):
        logging.info(f"Preparing simulation: {simType}")
        if edges is None:
            edges = self.edges
        predecessor = None
        if simType == "nvt":
            predecessor = "em"
        elif simType == "npt":
            predecessor = "nvt"
        elif simType == "em":
            predecessor = None
        else:
            raise ValueError(f"Unknown simType requested: {simType}")

        for edge in edges:
            logging.info(f"   Edge: {edge}")
            outLigPath = os.path.join(self.workPath, edge, "water")
            for state in self.states:
                for r in range(1, self.replicas + 1):
                    simpath = os.path.join(outLigPath, state, f"run{r}", simType)
                    if predecessor is None:
                        # Choose a fallback depending on what you actually have (for vacuum)
                        if os.path.exists(os.path.join(outLigPath, "ions.pdb")):
                            inStruct = os.path.join(outLigPath, "ions.pdb")
                        elif os.path.exists(os.path.join(outLigPath, "water.pdb")):
                            inStruct = os.path.join(outLigPath, "water.pdb")
                        elif os.path.exists(os.path.join(outLigPath, "box.pdb")):
                            inStruct = os.path.join(outLigPath, "box.pdb")
                        else:
                            inStruct = os.path.join(outLigPath, "init.pdb")
                    else:
                        inStruct = os.path.join(
                            outLigPath, state, f"run{r}", predecessor, "confout.gro"
                        )
                    if simType == "em":
                        mdpFile = (
                            os.path.join(self.mdpPath, "em_l0.mdp")
                            if state == "stateA"
                            else os.path.join(self.mdpPath, "em_l1.mdp")
                        )
                    elif simType == "nvt":
                        mdpFile = (
                            os.path.join(self.mdpPath, "NVT_l0.mdp")
                            if state == "stateA"
                            else os.path.join(self.mdpPath, "NVT_l1.mdp")
                        )
                    elif simType == "npt":
                        mdpFile = (
                            os.path.join(self.mdpPath, "NPT_l0.mdp")
                            if state == "stateA"
                            else os.path.join(self.mdpPath, "NPT_l1.mdp")
                        )
                    outTPR = os.path.join(simpath, "tpr.tpr")
                    mdout = os.path.join(simpath, "mdout.mdp")
                    topFile = os.path.join(outLigPath, "topol.top")

                    if not os.path.isfile(inStruct):
                        raise FileNotFoundError(
                            f"Cannot find input structure for {simType}: {inStruct}\n"
                            "Did you run the previous step or check naming?"
                        )

                    gmx.grompp(
                        f=mdpFile,
                        c=inStruct,
                        p=topFile,
                        o=outTPR,
                        maxwarn=4,
                        other_flags=f" -po {mdout}",
                    )
                    self._clean_backup_files(simpath)

        logging.info(f"Done preparing {simType} simulation(s).")

    def run_simulation_locally(self, edges=None, simType="em", bVerbose=False):
        """
        Runs EM/NVT/NPT or transitions. If 'cascade' is True, transitions are run in cascade mode.
        Otherwise, they are run in parallel/batch mode by default.
        """
        if edges is None:
            edges = self.edges
        if simType not in ["em", "nvt", "npt", "transitions"]:
            raise ValueError(f"Unknown simType for local run: {simType}")

        if simType in ["em", "nvt", "npt"]:
            logging.info(f"Running simulation locally: {simType}")
            self._run_regular_sim(edges, simType, bVerbose)
        else:
            # simType == 'transitions'
            if hasattr(self, "cascade") and self.cascade:
                logging.info("Running transitions in cascade mode.")
                self._run_multi_lambda_cascade(edges, bVerbose)
            else:
                logging.info("Running transitions in parallel batches.")
                self._run_multi_lambda_in_batches(edges, bVerbose)

    def _run_regular_sim(self, edges, simType, bVerbose):
        for edge in edges:
            for state in self.states:
                for r in range(1, self.replicas + 1):
                    base_simpath = os.path.join(
                        self.workPath, edge, "water", state, f"run{r}", simType
                    )
                    tpr = os.path.join(base_simpath, "tpr.tpr")
                    ener = os.path.join(base_simpath, "ener.edr")
                    confout = os.path.join(base_simpath, "confout.gro")
                    mdlog = os.path.join(base_simpath, "md.log")
                    trr = os.path.join(base_simpath, "traj.trr")
                    xtc = os.path.join(base_simpath, "traj.xtc")

                    if not os.path.isfile(tpr):
                        raise FileNotFoundError(
                            f"Cannot find {simType} TPR file at: {tpr}\n"
                            f"Did you forget to run prepare_simulation(simType='{simType}')?"
                        )

                    self._run_mdrun(tpr, ener, confout, mdlog, trr, bVerbose, xtc=xtc)
                    self._clean_backup_files(base_simpath)

    # -------------------------------------------------------
    # 6) Prepare & Run Multi-Lambda Transitions
    # -------------------------------------------------------
    def prepare_transitions(self, edges=None, bGenTpr=True):
        """
        Creates ti.tpr for each lambda window in transitions using the same final NPT confout.
        If running in cascade mode, you might skip this step and build TPRs on the fly.
        Otherwise, we generate all TPRs up front (parallel approach).
        """
        logging.info(
            f"Preparing multi-lambda transitions (bGenTpr={bGenTpr}). Cascade? {self.cascade}"
        )
        if edges is None:
            edges = self.edges
        if self.cascade:
            logging.info(
                "In cascade mode, you can skip prepare_transitions() if desired."
            )

        for edge in edges:
            # ligA, ligB = edges[edge]
            for state in self.states:
                if state == "stateA":
                    template_mdp = os.path.join(self.mdpPath, "template_ti_l0.mdp")
                    lam_schedule = " ".join([f"{x:.6f}" for x in self.lambdaStates])
                else:
                    template_mdp = os.path.join(self.mdpPath, "template_ti_l1.mdp")
                    lam_schedule = " ".join(
                        [f"{x:.6f}" for x in reversed(self.lambdaStates)]
                    )

                for r in range(1, self.replicas + 1):
                    transBase = os.path.join(
                        self.workPath, edge, "water", state, f"run{r}", "transitions"
                    )
                    npt_out = os.path.join(
                        self.workPath,
                        edge,
                        "water",
                        state,
                        f"run{r}",
                        "npt",
                        "confout.gro",
                    )
                    if not os.path.isfile(npt_out):
                        raise FileNotFoundError(f"NPT confout missing: {npt_out}")
                    for i in range(len(self.lambdaStates)):
                        lamFolder = os.path.join(transBase, f"lambda_{i}")
                        mdpFile = os.path.join(lamFolder, f"ti_lambda{i}.mdp")
                        with open(template_mdp, "r") as f:
                            mdp_txt = f.read()
                        # Replace placeholders
                        # mdp_txt = mdp_txt.replace('MOL_TYPE', f'IP_{ligA}')
                        mdp_txt = mdp_txt.replace(
                            "init_lambda_state        = REPLACE_LAMBDA",
                            f"init_lambda_state        = {i}",
                        )
                        mdp_txt = mdp_txt.replace("REPLACE_FEP_LAMBDAS", lam_schedule)
                        with open(mdpFile, "w") as f:
                            f.write(mdp_txt)

                        if bGenTpr and not self.cascade:
                            self._prepare_single_lambda_tpr(
                                lamFolder,
                                mdpFile,
                                npt_out,
                                toppath=os.path.join(self.workPath, edge, "water"),
                            )
        logging.info("Finished preparing transitions.")

    def _prepare_single_lambda_tpr(self, lamFolder, mdpFile, inStruct, toppath):
        logging.info(
            f"Preparing singele lamda tpr for folder: {lamFolder} with input structure: {inStruct} using topology: {toppath}"
        )
        if not os.path.isfile(inStruct):
            raise FileNotFoundError(
                f"Cannot find the input structure for transitions: {inStruct}\n"
                "Check that NPT was run and confout.gro is present."
            )
        topFile = os.path.join(toppath, "topol.top")
        if not os.path.isfile(topFile):
            raise FileNotFoundError(
                f"Cannot find topol.top for transitions: {topFile}\n"
                "Check that assemble_systems() and boxWaterIons() were run correctly."
            )
        outTPR = os.path.join(lamFolder, "ti.tpr")
        mdout = os.path.join(lamFolder, "mdout.mdp")
        gmx.grompp(
            f=mdpFile,
            c=inStruct,
            p=topFile,
            o=outTPR,
            maxwarn=4,
            other_flags=f" -po {mdout}",
        )
        logging.info("Finshed preparing singele lamda tpr!")

        self._clean_backup_files(lamFolder)

    def _run_multi_lambda_in_batches(self, edges, bVerbose=False):
        """
        Parallel/batch approach: all windows use the same NPT confout.
        We run them in chunks of 'windowsBatch' threads at once.
        """
        if edges is None:
            edges = self.edges
        for edge in edges:
            for state in self.states:
                for r in range(1, self.replicas + 1):
                    transBase = os.path.join(
                        self.workPath, edge, "water", state, f"run{r}", "transitions"
                    )
                    npt_out = os.path.join(
                        self.workPath,
                        edge,
                        "water",
                        state,
                        f"run{r}",
                        "npt",
                        "confout.gro",
                    )
                    if not os.path.isfile(npt_out):
                        raise FileNotFoundError(
                            f"Missing NPT confout for {state}/run{r}. Expected: {npt_out}"
                        )
                    numLams = len(self.lambdaStates)
                    iStart = 0
                    while iStart < numLams:
                        iEnd = min(iStart + self.windowsBatch, numLams)
                        commands = []
                        for i in range(iStart, iEnd):
                            lamFolder = os.path.join(transBase, f"lambda_{i}")
                            tpr = os.path.join(lamFolder, "ti.tpr")
                            confout = os.path.join(lamFolder, "confout.gro")
                            mdlog = os.path.join(lamFolder, "md.log")
                            trr = os.path.join(lamFolder, "traj.trr")
                            xtc = os.path.join(lamFolder, "traj.xtc")
                            dhdl = os.path.join(lamFolder, "dhdl.xvg")
                            ener = os.path.join(lamFolder, "ener.edr")
                            if not os.path.isfile(tpr):
                                raise FileNotFoundError(
                                    f"ti.tpr not found in {lamFolder}. Did you run prepare_transitions()?"
                                )
                            pinoffset = (i - iStart) * self.ntomp
                            cmd_list = [
                                "gmx",
                                "mdrun",
                                "-s",
                                tpr,
                                "-e",
                                ener,
                                "-c",
                                confout,
                                "-o",
                                trr,
                                "-g",
                                mdlog,
                                "-ntmpi",
                                str(self.ntmpi),
                                "-ntomp",
                                str(self.ntomp),
                                "-pin",
                                "on",
                                "-pinoffset",
                                str(pinoffset),
                                "-x",
                                xtc,
                                "-dhdl",
                                dhdl,
                            ]
                            cmd_str = " ".join(cmd_list)
                            commands.append((cmd_str, lamFolder))

                        # Run these in parallel with a ThreadPoolExecutor
                        with concurrent.futures.ThreadPoolExecutor(
                            max_workers=self.windowsBatch
                        ) as executor:
                            future_map = {}
                            for cmd_str, lamFolder in commands:
                                future = executor.submit(
                                    run_command, cmd_str, lamFolder
                                )
                                future_map[future] = lamFolder

                            for future in concurrent.futures.as_completed(future_map):
                                lamFolder = future_map[future]
                                try:
                                    out = future.result()
                                    if bVerbose and out:
                                        logging.info(out)
                                except Exception as e:
                                    logging.error(
                                        f"Error in lambda window {lamFolder}: {e}"
                                    )
                                self._clean_backup_files(lamFolder)

                        iStart += self.windowsBatch

    # -------------------------------------------------------
    # Cascade version for transitions: sequential windows per replicate/state
    # -------------------------------------------------------
    def _run_multi_lambda_cascade(self, edges, bVerbose=False):
        """
        Cascade approach: each window uses the final structure from the previous window.
        We can optionally rebuild TPR each time or rely on pre-built TPR if we called prepare_transitions.
        """
        logging.info("Running multi-lambda transitions in cascade mode...")
        if edges is None:
            edges = self.edges
        for edge in edges:
            for state in self.states:
                for r in range(1, self.replicas + 1):
                    transBase = os.path.join(
                        self.workPath, edge, "water", state, f"run{r}", "transitions"
                    )
                    input_structure = os.path.join(
                        self.workPath,
                        edge,
                        "water",
                        state,
                        f"run{r}",
                        "npt",
                        "confout.gro",
                    )
                    if not os.path.isfile(input_structure):
                        raise FileNotFoundError(
                            f"Missing NPT confout for {state}/run{r}. Expected at {input_structure}"
                        )
                    for i in range(len(self.lambdaStates)):
                        lamFolder = os.path.join(transBase, f"lambda_{i}")
                        # If you want to rebuild TPR on the fly for each step, do so here:
                        mdpFile = os.path.join(lamFolder, f"ti_lambda{i}.mdp")
                        if not os.path.isfile(mdpFile):
                            logging.warning(
                                f"No MDP file found at {mdpFile}, creating one on the fly."
                            )
                            # you'd fill in template logic here if skipping prepare_transitions
                        # Re-grompp using the current input structure
                        self._prepare_single_lambda_tpr(
                            lamFolder,
                            mdpFile,
                            input_structure,
                            toppath=os.path.join(self.workPath, edge, "water"),
                        )
                        tpr = os.path.join(lamFolder, "ti.tpr")
                        ener = os.path.join(lamFolder, "ener.edr")
                        mdlog = os.path.join(lamFolder, "md.log")
                        trr = os.path.join(lamFolder, "traj.trr")
                        xtc = os.path.join(lamFolder, "traj.xtc")
                        dhdl = os.path.join(lamFolder, "dhdl.xvg")
                        cmd_list = [
                            "gmx",
                            "mdrun",
                            "-s",
                            tpr,
                            "-e",
                            ener,
                            "-c",
                            os.path.join(lamFolder, "confout.gro"),
                            "-o",
                            trr,
                            "-g",
                            mdlog,
                            "-ntmpi",
                            str(self.ntmpi),
                            "-ntomp",
                            str(self.ntomp),
                            "-pin",
                            "on",
                            "-pinoffset",
                            "0",  # or could do something else
                        ]
                        if xtc:
                            cmd_list += ["-x", xtc]
                        if dhdl:
                            cmd_list += ["-dhdl", dhdl]
                        cmd_str = " ".join(cmd_list)
                        logging.info(
                            f"[CASCADE] Running window {i} with command: {cmd_str}"
                        )
                        try:
                            out = run_command(cmd_str, cwd=lamFolder)
                            if bVerbose and out:
                                logging.info(out)
                        except Exception as e:
                            logging.error(f"Error in cascade window {lamFolder}: {e}")
                            break
                        # Update input structure for the next window
                        new_input = os.path.join(lamFolder, "confout.gro")
                        if os.path.isfile(new_input):
                            input_structure = new_input
                        else:
                            logging.error(
                                f"Missing confout.gro from cascade step {i} in {lamFolder}!"
                            )
                            break
                        self._clean_backup_files(lamFolder)
        logging.info("DONE running multi-lambda transitions in cascade mode.")

    # -------------------------------------------------------
    # 7) Low-level run & cleanup
    # -------------------------------------------------------
    def _run_mdrun(
        self, tpr, ener, confout, mdlog, trr, bVerbose=False, xtc=None, dhdl=None
    ):
        cmd = [
            "gmx",
            "mdrun",
            "-s",
            tpr,
            "-e",
            ener,
            "-c",
            confout,
            "-o",
            trr,
            "-g",
            mdlog,
            "-ntmpi",
            "1",
        ]
        if xtc:
            cmd += ["-x", xtc]
        if dhdl:
            cmd += ["-dhdl", dhdl]
        cmd_str = " ".join(cmd)
        logging.debug(f"Running: {cmd_str}")
        out = run_command(cmd_str)
        if bVerbose and out:
            logging.info(out)

    def _clean_backup_files(self, path):
        backups = glob.glob(os.path.join(path, "*#"))
        for b in backups:
            os.remove(b)

    # -------------------------------------------------------
    # 8) Basic Analysis Example
    # -------------------------------------------------------
    def run_analysis(self, edges=None, bVerbose=False, start_time=200):
        """
        Performs a basic BAR analysis on the non-equilibrium free energy data by invoking
        GROMACS' gmx bar command on the dhdl.xvg files from each lambda window.
        """
        if edges is None:
            edges = self.edges
        for edge in edges:
            for r in range(1, self.replicas + 1):
                for state in self.states:
                    logging.info(f"Analyzing edge={edge}, replicate={r}, state={state}")
                    transitions_dir = os.path.join(
                        self.workPath, edge, "water", state, f"run{r}", "transitions"
                    )
                    dhdl_files = []
                    for i in range(len(self.lambdaStates)):
                        lam_dir = os.path.join(transitions_dir, f"lambda_{i}")
                        dhdl_file = os.path.join(lam_dir, "dhdl.xvg")
                        if os.path.isfile(dhdl_file):
                            dhdl_files.append(dhdl_file)
                        else:
                            logging.warning(f"Missing dhdl.xvg in {lam_dir}")
                    if not dhdl_files:
                        logging.warning(
                            f"No dhdl.xvg files for replicate={r}, state={state}, skipping."
                        )
                        continue
                    file_list_str = " ".join(dhdl_files)
                    cmd = f"gmx bar -b {start_time} -f {file_list_str}"
                    logging.info(f"Running BAR analysis: {cmd}")
                    try:
                        out = run_command(cmd)
                    except Exception as e:
                        logging.error(
                            f"Error during BAR analysis for edge {edge}, replicate={r}, state={state}: {e}"
                        )
                        continue
                    analysis_dir = os.path.join(
                        self.workPath, edge, "water", f"analyse_r{r}", state
                    )
                    create_folder(analysis_dir)
                    out_file = os.path.join(analysis_dir, "bar_results.txt")
                    with open(out_file, "w") as f:
                        f.write(out)
                    if bVerbose and out:
                        lines = out.splitlines()
                        n_lines_to_keep = 20 + 2 * (len(self.lambdaStates) - 1)
                        if len(lines) > n_lines_to_keep:
                            lines = lines[-n_lines_to_keep:]
                        logging.info("\n".join(lines))
        logging.info("Analysis done.")
