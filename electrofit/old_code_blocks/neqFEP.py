import pmx
from pmx.utils import create_folder
from pmx import gmx, ligand_alchemy, jobscript
import os
import subprocess
import glob
import pandas as pd
import numpy as np

class neqFEP:
    """Class for setting up non-equilibrium free energy calculations (neFEP) in water.

    This workflow runs:
      1) Atom mapping / hybrid structure generation
      2) Box creation & solvation
      3) Energy minimization, followed by two‐step equilibration:
           a) NVT equilibration (NVT_l0/NVT_l1)
           b) NPT equilibration (NPT_l0/NPT_l1)
      4) Production equilibration run (prod) for sampling equilibrium conditions
      5) Transitions: extraction of snapshots from the production run for free energy sampling
      6) Analysis of free energy differences

    Attributes
    ----------
    workPath : str
        Main working directory.
    mdpPath : str
        Directory holding .mdp files for different simulation phases:
          - em_l0.mdp, em_l1.mdp  (minimization)
          - NVT_l0.mdp, NVT_l1.mdp (first equilibration step)
          - NPT_l0.mdp, NPT_l1.mdp (second equilibration step)
          - PROD_l0.mdp, PROD_l1.mdp (production equilibration run)
          - ti_l0.mdp, ti_l1.mdp   (transitions)
    ligandPath : str
        Directory with ligand structures and topologies.
    edges : dict or list
        Pairs of ligands to be transformed, e.g. [ ["lig1", "lig2"] ].
    """

    def __init__(self, **kwargs):
        # Initialize PMX gmx path
        gmx.set_gmxlib()

        # DataFrame placeholders for results
        self.resultsAll = pd.DataFrame()
        self.resultsSummary = pd.DataFrame()

        # Default paths
        self.workPath = './'
        self.mdpPath = './mdp'
        self.ligandPath = None  # Must be provided by user
        self.edges = {}         # Edges: e.g., { "edge_ligA_ligB": ["ligA","ligB"] }

        # This workflow is water-only (no protein)
        self.replicas = 3
        # Define simulation types.
        self.simTypes = ['em', 'nvt', 'npt', 'prod', 'transitions']
        self.states = ['stateA', 'stateB']
        self.thermCycleBranches = ['water']

        # Simulation parameters
        self.ff = 'amber14sb.ff'
        self.boxshape = 'dodecahedron'
        self.boxd = 1.2
        self.water = 'tip3p'
        self.conc = 0.15
        self.pname = 'NA'
        self.nname = 'CL'

        # Job submission defaults (if needed)
        self.JOBqueue = 'SLURM'
        self.JOBsimtime = 24   # hours
        self.JOBsimcpu = 8
        self.JOBbGPU = True
        self.JOBmodules = []
        self.JOBsource = []
        self.JOBexport = []
        self.JOBgmx = 'gmx mdrun'

        # Allow user overrides
        for key, val in kwargs.items():
            setattr(self, key, val)

    # -------------------------------------------------------
    # 1) Directory Setup
    # -------------------------------------------------------
    def prepareFreeEnergyDir(self):
        # Ensure absolute paths
        self.workPath = os.path.abspath(self.workPath)
        self.mdpPath = os.path.abspath(self.mdpPath)
        if self.ligandPath is None:
            raise ValueError("ligandPath must be defined for ligand transformations in water.")
        self.ligandPath = os.path.abspath(self.ligandPath)

        # Create main work directory
        create_folder(self.workPath)

        # Read ligand directories and edges
        self._read_ligands()
        self._read_edges()

        # Create folder structure for each edge, state, replicate, and simulation type.
        self._create_folder_structure()

        self._print_summary()
        self._print_folder_structure()
        print("DONE")
    
    def _get_specific_path(self, edge=None, bHybridStrTop=False, wp=None, state=None, r=None, sim=None):
        if edge is None:
            return self.workPath       
        edgepath = f"{self.workPath}/{edge}"
        if bHybridStrTop:
            return f"{edgepath}/hybridStrTop"
        if wp is None:
            return edgepath
        wppath = f"{edgepath}/{wp}"
        if state is None:
            return wppath
        statepath = f"{wppath}/{state}"
        if r is None:
            return statepath
        runpath = f"{statepath}/run{r}"
        if sim is None:
            return runpath
        return f"{runpath}/{sim}"
                
    def _read_path(self, path):
        return os.path.abspath(path)

    def _read_ligands(self):
        """Scan ligandPath for folders and store them."""
        self.ligands = {}
        lig_folders = glob.glob(f"{self.ligandPath}/*")
        for folder in lig_folders:
            basename = os.path.basename(folder)
            # Assume ligand folders start with 'IP_'
            lname = basename[3:] if basename.startswith('IP_') else basename
            self.ligands[lname] = os.path.abspath(folder)

    def _read_edges(self):
        """Convert edges list to dictionary if needed."""
        if isinstance(self.edges, list):
            foo = {}
            for e in self.edges:
                key = f"edge_{e[0]}_{e[1]}"
                foo[key] = e
            self.edges = foo
        elif isinstance(self.edges, str) and os.path.isfile(self.edges):
            self.edges = "Edges read from file"

    def _create_folder_structure(self):
        """Create directories for each edge, state, replicate, and simType."""
        for edge in self.edges:
            edgepath = f"{self.workPath}/{edge}"
            create_folder(edgepath)
            # Folder for hybrid ligand structures
            create_folder(f"{edgepath}/hybridStrTop")
            # We only use 'water'
            for wp in self.thermCycleBranches:
                wppath = f"{edgepath}/{wp}"
                create_folder(wppath)
                for state in self.states:
                    statepath = f"{wppath}/{state}"
                    create_folder(statepath)
                    for r in range(1, self.replicas + 1):
                        runpath = f"{statepath}/run{r}"
                        create_folder(runpath)
                        # Create subfolders for each simulation type: em, nvt, npt, prod, transitions
                        for sim in self.simTypes:
                            create_folder(f"{runpath}/{sim}")

    def _print_summary(self):
        print("\n---------------------\nSummary of the setup:\n---------------------\n")
        print(f"   workPath:   {self.workPath}")
        print(f"   mdpPath:    {self.mdpPath}")
        print(f"   ligandPath: {self.ligandPath}")
        print(f"   replicas:   {self.replicas}")
        print("   edges:")
        for e in self.edges:
            print(f"        {e}")

    def _print_folder_structure(self):
        print("\n---------------------\nDirectory structure:\n---------------------\n")
        print(f"{self.workPath}/")
        print("|-- edge_X_Y")
        print("|   |-- hybridStrTop")
        print("|   |-- water")
        print("|       |-- stateA")
        print("|       |   |-- run1/2/3")
        print("|       |       |-- em/nvt/npt/prod/transitions")
        print("|       |-- stateB")
        print("|           |-- run1/2/3")
        print("|               |-- em/nvt/prod/npt/transitions")
        print("|-- edge_...")
        print()

    # -------------------------------------------------------
    # 2) Atom Mapping & Hybrid Structure/Topology
    # -------------------------------------------------------
    def atom_mapping(self, edges=None, bVerbose=False):
        print("-----------------------")
        print("Performing atom mapping")
        print("-----------------------")
        if edges is None:
            edges = self.edges
        for edge in edges:
            print(edge)
            ligA, ligB = self.edges[edge][0], self.edges[edge][1]
            ligApath = f"{self.ligandPath}/IP_{ligA}"
            ligBpath = f"{self.ligandPath}/IP_{ligB}"
            ligAname = f"IP_{ligA}"
            ligBname = f"IP_{ligB}"
            outpath = f"{self.workPath}/{edge}/hybridStrTop"

            i1 = f"{ligApath}/{ligAname}.pdb"
            i2 = f"{ligBpath}/{ligBname}.pdb"
            o1 = f"{outpath}/pairs1.dat"
            o2 = f"{outpath}/pairs2.dat"
            opdb1 = f"{outpath}/out_pdb1.pdb"
            opdb2 = f"{outpath}/out_pdb2.pdb"
            opdbm1 = f"{outpath}/out_pdbm1.pdb"
            opdbm2 = f"{outpath}/out_pdbm2.pdb"
            score = f"{outpath}/score.dat"
            log = f"{outpath}/mapping.log"

            cmd = [
                "pmx", "atomMapping",
                "-i1", i1,
                "-i2", i2,
                "-o1", o1,
                "-o2", o2,
                "-opdb1", opdb1,
                "-opdb2", opdb2,
                "-opdbm1", opdbm1,
                "-opdbm2", opdbm2,
                "-score", score,
                "-log", log,
                "--H2Hpolar",
                "--no-mcs"
            ]
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self._be_verbose(process, bVerbose)
            process.wait()
        print("DONE")

    def hybrid_structure_topology(self, edges=None, bVerbose=False):
        print("----------------------------------")
        print("Creating hybrid structure/topology")
        print("----------------------------------")
        if edges is None:
            edges = self.edges
        for edge in edges:
            print(edge)
            ligA, ligB = self.edges[edge][0], self.edges[edge][1]
            ligApath = f"{self.ligandPath}/IP_{ligA}"
            ligBpath = f"{self.ligandPath}/IP_{ligB}"
            ligAname = f"IP_{ligA}"
            ligBname = f"IP_{ligB}"
            outpath = f"{self.workPath}/{edge}/hybridStrTop"

            i1 = f"{ligApath}/{ligAname}.pdb"
            i2 = f"{ligBpath}/{ligBname}.pdb"
            itp1 = f"{ligApath}/{ligAname}_GMX.itp"
            itp2 = f"{ligBpath}/{ligBname}_GMX.itp"
            pairs = f"{outpath}/pairs1.dat"
            oA = f"{outpath}/mergedA.pdb"
            oB = f"{outpath}/mergedB.pdb"
            oitp = f"{outpath}/merged.itp"
            offitp = f"{outpath}/ffmerged.itp"
            log = f"{outpath}/hybrid.log"

            cmd = [
                "pmx", "ligandHybrid",
                "-i1", i1,
                "-i2", i2,
                "-itp1", itp1,
                "-itp2", itp2,
                "-pairs", pairs,
                "-oA", oA,
                "-oB", oB,
                "-oitp", oitp,
                "-offitp", offitp,
                "-log", log
            ]
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self._be_verbose(process, bVerbose)
            process.wait()
        print("DONE")

    def _be_verbose(self, process, bVerbose):
        out = process.communicate()
        if bVerbose:
            for line in out[0].decode().splitlines():
                print(line)
        for err in out[1].decode().splitlines():
            print(err)

    # -------------------------------------------------------
    # 3) System Assembly (Water Only)
    # -------------------------------------------------------
    def assemble_systems(self, edges=None):
        print("Assembling the systems (water only).")
        if edges is None:
            edges = self.edges
        for edge in edges:
            print(f"Assembling system for edge: {edge}")
            ligA, ligB = self.edges[edge]
            hybridStrTopPath = os.path.join(self.workPath, edge, "hybridStrTop")
            outLigPath = os.path.join(self.workPath, edge, "water")

            # 1) Create init.pdb from mergedA
            self._make_clean_pdb(
                os.path.join(hybridStrTopPath, "mergedA.pdb"),
                os.path.join(outLigPath, "init.pdb")
            )

            # 2) Merge force field files
            ffitpOut = os.path.join(hybridStrTopPath, "ffmerged.itp")
            ffitpIn1 = os.path.join(self.ligandPath, f"IP_{ligA}", f"ff_{ligA}.itp")
            ffitpIn2 = os.path.join(self.ligandPath, f"IP_{ligB}", f"ff_{ligB}.itp")
            ffitpIn3 = os.path.join(hybridStrTopPath, "ffmerged.itp")
            ligand_alchemy._merge_FF_files(ffitpOut, ffsIn=[ffitpIn1, ffitpIn2, ffitpIn3])

            # 3) Create topol.top
            ligTopFname = os.path.join(outLigPath, "topol.top")
            ligFFitp = os.path.join(hybridStrTopPath, "ffmerged.itp")
            ligItp   = os.path.join(hybridStrTopPath, "merged.itp")
            ligposresItp = os.path.join(self.ligandPath, "posres", "posres.itp")
            itps = {"FFitp": ligFFitp, "ligItp": ligItp, "posresItp": ligposresItp}
            self._create_top(ligTopFname, itps, f"IP_{ligA}")
        print("System assembly complete.")

    def _make_clean_pdb(self, fnameIn, fnameOut):
        with open(fnameIn, 'r') as f:
            lines = f.readlines()
        keep = [l for l in lines if l.startswith('ATOM') or l.startswith('HETATM')]
        with open(fnameOut, 'w') as w:
            w.writelines(keep)

    def _create_top(self, fname, itp, ligAname, systemName="ligand in water"):
        with open(fname, 'w') as fp:
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
        print("----------------")
        print("Box, water, ions")
        print("----------------")
        if edges is None:
            edges = self.edges
        for edge in edges:
            print(edge)
            outLigPath = f"{self.workPath}/{edge}/water"
            # Define box
            if bBox:
                inStr = f"{outLigPath}/init.pdb"
                outStr = f"{outLigPath}/box.pdb"
                gmx.editconf(inStr, o=outStr, bt=self.boxshape, d=self.boxd)
            # Solvate
            if bSolvate:
                inBox = f"{outLigPath}/box.pdb"
                outWater = f"{outLigPath}/water.pdb"
                top = f"{outLigPath}/topol.top"
                gmx.solvate(inBox, cs='spc216.gro', p=top, o=outWater)
            # Add ions
            if bIon:
                inWater = f"{outLigPath}/water.pdb"
                outIons = f"{outLigPath}/ions.pdb"
                mdp = f"{self.mdpPath}/em_l0.mdp"
                tpr = f"{outLigPath}/tpr.tpr"
                top = f"{outLigPath}/topol.top"
                mdout = f"{outLigPath}/mdout.mdp"
                gmx.grompp(f=mdp, c=inWater, p=top, o=tpr, maxwarn=4,
                           other_flags=f" -po {mdout}")
                gmx.genion(s=tpr, p=top, o=outIons, conc=self.conc, neutral=True,
                           other_flags=f" -pname {self.pname} -nname {self.nname}")
        print("DONE")

    # -------------------------------------------------------
    # 5) Prepare & Run Simulations
    # -------------------------------------------------------
    def prepare_simulation(self, edges=None, simType='em'):
        """
        Create .tpr files for a given simulation type.
        For non-equilibrium FEP we now use:
         - 'em' for energy minimization,
         - 'nvt' for first equilibration (NVT),
         - 'npt' for second equilibration (NPT),
         - 'prod' for the production equilibration run,
         - 'transitions' for extracting snapshots from the production run.
        The input structure for each step comes from the previous step:
         - nvt uses output from em,
         - npt uses output from nvt,
         - prod uses output from npt,
         - transitions uses output from prod.
        """
        print("-----------------------------------------")
        print(f"Preparing simulation: {simType}")
        print("-----------------------------------------")
        if edges is None:
            edges = self.edges

        # Map simType to its predecessor
        predecessor = None
        if simType == "nvt":
            predecessor = "em"
        elif simType == "npt":
            predecessor = "nvt"
        elif simType == "prod":
            predecessor = "npt"
        elif simType == "transitions":
            predecessor = "prod"

        for edge in edges:
            print(edge)
            outLigPath = f"{self.workPath}/{edge}/water"
            for state in self.states:
                for r in range(1, self.replicas+1):
                    simpath = f"{outLigPath}/{state}/run{r}/{simType}"
                    emp = None if simType == "em" else f"{outLigPath}/{state}/run{r}/{predecessor}"
                    self._prepare_single_tpr(simpath=simpath,
                                             toppath=outLigPath,
                                             state=state,
                                             simType=simType,
                                             empath=emp)
        print("DONE")

    def _prepare_single_tpr(self, simpath, toppath, state, simType, empath=None, frameNum=0):
        """
        Creates a single .tpr file using grompp.
        The MDP file is chosen based on simType and state.
        """
        if simType == "em":
            mdpFile = f"{self.mdpPath}/em_l0.mdp" if state == "stateA" else f"{self.mdpPath}/em_l1.mdp"
            inStruct = f"{toppath}/ions.pdb"
            outTPR = f"{simpath}/tpr.tpr"
        elif simType == "nvt":
            mdpFile = f"{self.mdpPath}/NVT_l0.mdp" if state == "stateA" else f"{self.mdpPath}/NVT_l1.mdp"
            inStruct = f"{empath}/confout.gro"
            outTPR = f"{simpath}/tpr.tpr"
        elif simType == "npt":
            mdpFile = f"{self.mdpPath}/NPT_l0.mdp" if state == "stateA" else f"{self.mdpPath}/NPT_l1.mdp"
            inStruct = f"{empath}/confout.gro"
            outTPR = f"{simpath}/tpr.tpr"
        elif simType == "prod":
            mdpFile = f"{self.mdpPath}/PROD_l0.mdp" if state == "stateA" else f"{self.mdpPath}/PROD_l1.mdp"
            inStruct = f"{empath}/confout.gro"
            outTPR = f"{simpath}/tpr.tpr"
        elif simType == "transitions":
            mdpFile = f"{self.mdpPath}/ti_l0.mdp" if state == "stateA" else f"{self.mdpPath}/ti_l1.mdp"
            inStruct = f"{simpath}/frame{frameNum}.gro"
            outTPR = f"{simpath}/ti{frameNum}.tpr"
        else:
            raise ValueError("Unknown simType.")

        mdout = f"{simpath}/mdout.mdp"
        topFile = f"{toppath}/topol.top"

        gmx.grompp(f=mdpFile, c=inStruct, p=topFile, o=outTPR, maxwarn=4,
                   other_flags=f" -po {mdout}")
        self._clean_backup_files(simpath)

    def run_simulation_locally(self, edges=None, simType='em', bVerbose=False):
        """Run gmx mdrun locally for the specified simType."""
        print("-------------------------------------------")
        print(f"Run simulation locally: {simType}")
        print("-------------------------------------------")
        if edges is None:
            edges = self.edges
        for edge in edges:
            for state in self.states:
                for r in range(1, self.replicas+1):
                    base_simpath = f"{self.workPath}/{edge}/water/{state}/run{r}/{simType}"
                    if simType == "transitions":
                        for i in range(1, 81):
                            simpath = base_simpath
                            tpr     = f"{simpath}/ti{i}.tpr"
                            ener    = f"{simpath}/ener{i}.edr"
                            confout = f"{simpath}/frame{i}.gro"
                            mdlog   = f"{simpath}/md{i}.log"
                            trr     = f"{simpath}/traj{i}.trr"
                            xtc     = f"{simpath}/traj{i}.xtc"
                            dhdl    = f"{simpath}/dhdl{i}.xvg"
                            self._run_mdrun(tpr, ener, confout, mdlog, trr, bVerbose, xtc=xtc, dhdl=dhdl)
                            self._clean_backup_files(simpath)
                    else:
                        simpath = base_simpath
                        tpr     = f"{simpath}/tpr.tpr"
                        ener    = f"{simpath}/ener.edr"
                        confout = f"{simpath}/confout.gro"
                        mdlog   = f"{simpath}/md.log"
                        trr     = f"{simpath}/traj.trr"
                        xtc     = f"{simpath}/traj.xtc"
                        self._run_mdrun(tpr, ener, confout, mdlog, trr, bVerbose, xtc=xtc)
                        self._clean_backup_files(simpath)
        print("DONE")

    def _run_mdrun(self, tpr, ener, confout, mdlog, trr, bVerbose=False, cpo=None, xtc=None, dhdl=None):
        cmd = [
            "gmx", "mdrun",
            "-s", tpr,
            "-e", ener,
            "-c", confout,
            "-o", trr,
            "-g", mdlog,
            "-ntmpi", "1",
        ]
        if xtc:
            cmd.extend(["-x", xtc])
        if dhdl:
            cmd.extend(["-dhdl", dhdl])
        
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self._be_verbose(process, bVerbose)
        process.wait()

    def _clean_backup_files(self, path):
        backups = glob.glob(f"{path}/*#")
        for b in backups:
            os.remove(b)

    # -------------------------------------------------------
    # 6) Prepare Transitions
    # -------------------------------------------------------
    def prepare_transitions(self, edges=None, bGenTpr=True, startTime=2250):
        """
        For non-equilibrium FEP, production runs (or slow growth runs) are created after the production equilibration stage.
        This method extracts frames from the production run and prepares .tpr files for transitions.
        """
        print("---------------------")
        print("Preparing transitions")
        print("---------------------")
        if edges is None:
            edges = self.edges
        for edge in edges:
            outLigPath = f"{self.workPath}/{edge}/water"
            for state in self.states:
                for r in range(1, self.replicas+1):
                    # Use production output as input for transitions
                    eqpath = f"{outLigPath}/{state}/run{r}/prod"
                    tipath = f"{outLigPath}/{state}/run{r}/transitions"
                    self._extract_snapshots(eqpath, tipath, startTime)
                    if bGenTpr:
                        for i in range(1, 81):
                            self._prepare_single_tpr(simpath=tipath,
                                                     toppath=outLigPath,
                                                     state=state,
                                                     simType="transitions",
                                                     empath=None,
                                                     frameNum=i)
        print("DONE")

    def _extract_snapshots(self, eqpath, tipath, startTime):
        tpr = f"{eqpath}/tpr.tpr"
        trr = f"{eqpath}/traj.trr"
        frame = f"{tipath}/frame.gro"
        gmx.trjconv(s=tpr, f=trr, o=frame, sep=True, ur='compact', pbc='mol',
                    other_flags=f" -b {startTime}")
        # Optionally rename frame0.gro to frame80.gro if desired
        zeroFrame = f"{tipath}/frame0.gro"
        eightyFrame = f"{tipath}/frame80.gro"
        if os.path.exists(zeroFrame):
            os.rename(zeroFrame, eightyFrame)
        self._clean_backup_files(tipath)

    # -------------------------------------------------------
    # 7) Analysis
    # -------------------------------------------------------
    def run_analysis(self, edges=None, bVerbose=False):
        """
        Gathers production (.xvg) files from stateA and stateB transitions,
        and calls pmx analyse to compute free energy differences.
        """
        print("----------------")
        print("Running analysis")
        print("----------------")
        if edges is None:
            edges = self.edges
        for edge in edges:
            print(edge)
            for r in range(1, self.replicas+1):
                analysispath = f"{self.workPath}/{edge}/water/analyse{r}"
                create_folder(analysispath)
                stateApath = f"{self.workPath}/{edge}/water/stateA/run{r}/transitions"
                stateBpath = f"{self.workPath}/{edge}/water/stateB/run{r}/transitions"
                self._run_analysis_script(analysispath, stateApath, stateBpath, bVerbose)
        print("DONE")

    def _run_analysis_script(self, analysispath, stateApath, stateBpath, bVerbose):
        fA = ' '.join(glob.glob(f"{stateApath}/*.xvg"))
        fB = ' '.join(glob.glob(f"{stateBpath}/*.xvg"))
        if not fA or not fB:
            print("No .xvg files found for analysis. Skipping.")
            return
        oA = f"{analysispath}/integ0.dat"
        oB = f"{analysispath}/integ1.dat"
        wplot = f"{analysispath}/wplot.png"
        outRes = f"{analysispath}/results.txt"
        cmd = f"pmx analyse -fA {fA} -fB {fB} -o {outRes} -oA {oA} -oB {oB} -w {wplot} -t 310 -b 100"
        subprocess.call(cmd, shell=True)
        if bVerbose and os.path.exists(outRes):
            with open(outRes, 'r') as fp:
                print(fp.read())

    def analysis_summary(self, edges=None):
        """
        Aggregates results for each replicate and summarizes the free energy differences.
        """
        if edges is None:
            edges = self.edges
        for edge in edges:
            for r in range(1, self.replicas+1):
                analysispath = f"{self.workPath}/{edge}/water/analyse{r}"
                resultsfile = f"{analysispath}/results.txt"
                if os.path.exists(resultsfile):
                    res = self._read_neq_results(resultsfile)
                    self._fill_resultsAll(res, edge, r)
        self._summarize_results(edges)
        self.save_results_txt()

    def save_results_txt(self):
        """Save resultsAll and resultsSummary as text files in an 'analysis_all' folder."""
        out_folder = os.path.join(self.workPath, "analysis_all")
        create_folder(out_folder)
        results_all_file = os.path.join(out_folder, "resultsAll.txt")
        results_summary_file = os.path.join(out_folder, "resultsSummary.txt")
        self.resultsAll.to_csv(results_all_file, sep='\t', index=True)
        self.resultsSummary.to_csv(results_summary_file, sep='\t', index=True)
        print(f"Saved resultsAll to {results_all_file}")
        print(f"Saved resultsSummary to {results_summary_file}")

    def _read_neq_results(self, fname):
        with open(fname, 'r') as fp:
            lines = fp.readlines()
        outvals = [0, 0, 0, 0, 0]  # [frames0->1, frames1->0, dG, err_analytical, err_bootstrap]
        for line in lines:
            line = line.strip()
            if 'BAR: dG' in line:
                tokens = line.split()
                dG = float(tokens[3]) if tokens[2] == "=" else float(tokens[2])
                outvals[2] = dG
            elif 'BAR: Std Err (bootstrap)' in line:
                tokens = line.split()
                try:
                    value = float(tokens[-1])
                except ValueError:
                    for token in tokens:
                        try:
                            value = float(token)
                            break
                        except ValueError:
                            continue
                outvals[4] = value
            elif 'BAR: Std Err (analytical)' in line:
                tokens = line.split()
                try:
                    value = float(tokens[-1])
                except ValueError:
                    for token in tokens:
                        try:
                            value = float(token)
                            break
                        except ValueError:
                            continue
                outvals[3] = value
            elif '0->1' in line:
                tokens = line.split()
                outvals[0] = int(tokens[-1])
            elif '1->0' in line:
                tokens = line.split()
                outvals[1] = int(tokens[-1])
        return outvals

    def _fill_resultsAll(self, res, edge, r):
        rowName = f"{edge}_water_{r}"
        framesA, framesB, dG, errAna, errBoot = res
        self.resultsAll.loc[rowName, 'val'] = dG
        self.resultsAll.loc[rowName, 'err_analyt'] = errAna
        self.resultsAll.loc[rowName, 'err_boot'] = errBoot
        self.resultsAll.loc[rowName, 'framesA'] = framesA
        self.resultsAll.loc[rowName, 'framesB'] = framesB

    def _summarize_results(self, edges):
        for edge in edges:
            replicate_rows = [f"{edge}_water_{r}" for r in range(1, self.replicas+1)
                              if f"{edge}_water_{r}" in self.resultsAll.index]
            dGs = self.resultsAll.loc[replicate_rows, 'val'].values
            errAnalytical = self.resultsAll.loc[replicate_rows, 'err_analyt'].values
            errBoot = self.resultsAll.loc[replicate_rows, 'err_boot'].values
            if len(dGs) == 0:
                continue
            if len(dGs) == 1:
                avg_dG = dGs[0]
                avg_errAna = errAnalytical[0]
                avg_errBoot = errBoot[0]
            else:
                avg_dG = np.mean(dGs)
                avg_errAna  = np.sqrt(np.mean(errAnalytical**2))
                avg_errBoot = np.sqrt(np.mean(errBoot**2))
            rowName = f"{edge}_water"
            self.resultsAll.loc[rowName, 'val'] = avg_dG
            self.resultsAll.loc[rowName, 'err_analyt'] = avg_errAna
            self.resultsAll.loc[rowName, 'err_boot'] = avg_errBoot
            self.resultsSummary.loc[edge, 'dG_water'] = avg_dG
            self.resultsSummary.loc[edge, 'errAna_water'] = avg_errAna
            self.resultsSummary.loc[edge, 'errBoot_water'] = avg_errBoot