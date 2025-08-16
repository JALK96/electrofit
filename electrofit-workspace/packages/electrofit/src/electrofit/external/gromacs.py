import logging
import os
import subprocess

import matplotlib.pyplot as plt
import pandas as pd
from electrofit.cli.run_commands import run_command
from electrofit.io.files import (
    include_ff,
    include_ions,
    include_tip3p,
    remove_defaults_section_lines,
    replace_posres_in_file,
    strip_extension,
)
from electrofit.logging import setup_logging
from electrofit.cli.safe_run import ensure_finalized
from electrofit.scratch.manager import (
    setup_scratch_directory,
)


def plot_svg(svg):
    # Suppress logging
    logging.disable(logging.CRITICAL)
    name = strip_extension(svg)

    df = pd.read_csv(svg, sep=r"\s+", header=None, names=["time", name])

    plt.figure(figsize=(8, 6))
    plt.plot(
        df["time"], df[name], color="darkblue", linestyle="-", label=f"{name}"
    )
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel(f"{name}", fontsize=14)
    plt.axhline(0, color="black", linewidth=0.8, linestyle="--")
    plt.tight_layout()
    plt.savefig(f"{name}.pdf", format="pdf")
    plt.close()

    # Re-enable logging if needed (optional)
    logging.disable(logging.NOTSET)


def set_up_production(
    m_gro,
    MDP_dir,
    base_scratch_dir,
    molecule_name,
    box_type="dodecahedron",
    cation="NA",
    anion="CL",
    d="1.2",
    conc="0.15",
    exit_screen=True,
    ff="amber14sb.ff",
    *,
    threads: int | None = None,
    pin: bool | None = None,
):
    """
    Set up and execute a production molecular dynamics (MD) simulation using GROMACS following input generation by acpype.

    This function automates the preparation and execution of a molecular dynamics simulation by performing the following steps:

    1. **File Preparation**:
       - Strips the extension from the provided `.gro` file to derive base filenames for topology (`.itp` and `.top`) files.
       - Modifies the topology file by replacing position restraints (`POSRES_LIG` with `POSRES`).
       - Includes TIP3P water and ion parameters into the topology file.

    2. **Scratch Directory Setup**:
       - Creates a scratch directory for processing data and running simulations, copying necessary input files into it.

    3. **System Configuration**:
       - Defines a simulation box with the ligand at the center, ensuring sufficient distance (`d`) from the edges. The box type can be specified (default is `'dodecahedron'`).
       - Solvates the system with water molecules.
       - Adds ions to neutralize the system and set the specified ion concentration (`conc`).

    4. **Energy Minimization and Equilibration**:
       - Performs energy minimization to relieve any steric clashes or unfavorable interactions.
       - Conducts NVT (constant Number of particles, Volume, and Temperature) and NPT (constant Number of particles, Pressure, and Temperature) equilibration steps to stabilize the system.

    5. **Production Run**:
       - Prepares input files for the production MD run.
       - Executes the production MD simulation.

    6. **Post-processing**:
       - Centers the molecule in the simulation box.
       - Calculates the minimum distance between the molecule and its periodic images.
       - Computes the radius of gyration to assess the compactness of the molecule.

    7. **Cleanup**:
       - Finalizes the scratch directory by moving or copying output files back to the original directory and cleaning up temporary files.

    Args:
        m_gro (str):
            Path to the GROMACS `.gro` input file after RESP fitting.
            Example: `"input.gro"`

        MDP_dir (str):
            Path to the directory containing GROMACS MDP (Molecular Dynamics Parameter) files.
            Example: `"/path/to/MDP"`

        base_scratch_dir (str):
            Base directory for creating a scratch space where data will be processed and simulations will be performed.
            Example: `"/scratch/user/tmp"`

        box_type (str, optional):
            Type of the simulation box to use. Default is `"dodecahedron"`.
            Examples include `"cubic"`, `"triclinic"`, etc.

        cation (str, optional):
            The cation to add to the system to neutralize the charge. Default is `"NA"`.
            Ensure that the cation name matches the force field parameters.

        anion (str, optional):
            The anion to add to the system to neutralize the charge. Default is `"CL"`.
            Ensure that the anion name matches the force field parameters.

        d (str or float, optional):
            Distance (in nanometers) from the molecule to the edge of the simulation box.
            Default is `"1.2"`. This defines the minimum distance between any atom of the molecule and the box edge.

        conc (str or float, optional):
            Ion concentration (in mol/L) to add to the system.
            Default is `"0.15"`. This sets the molar concentration of ions in the simulation.
        
        ff (str, optional):
            Force field folder used for includes (e.g., `tip3p.itp`, `ions.itp`, `forcefield.itp`).
            Default is `"amber14sb.ff"`
        threads (int | None, optional): Number of threads for gmx mdrun. If None, omit -nt.
        pin (bool | None, optional): Whether to add `-pin on/off`. If None, omit -pin.

    Raises:
        FileNotFoundError:
            If any of the specified input files (`m_gro`, MDP files) or directories (`MDP_dir`, `base_scratch_dir`) do not exist.

        RuntimeError:
            If any of the GROMACS commands executed within the function fail. This includes errors during box definition, solvation, ion addition, energy minimization, equilibration, or production runs.

        Exception:
            For any other unforeseen errors that may occur during the execution of the function.

    Returns:
        None

    Example:
        ```python
        set_up_production(
            m_gro="IP_010101_resp_GMX.gro",
            MDP_dir="/path/to/MDP",
            base_scratch_dir="/scratch/user/tmp",
            box_type="dodecahedron",
            cation="NA",
            anion="CL",
            d="1.2",
            conc="0.15"
        )
        ```

    Notes:
        - Ensure that all GROMACS executables (`gmx`) are accessible in the system's PATH.
        - The `MDP_dir` must contain the necessary MDP files for each simulation step (`em_steep.mdp`, `NVT.mdp`, `NPT.mdp`, `Production.mdp`).
        - The parameters `box_type`, `cation`, `anion`, `d`, and `conc` allow customization of the simulation setup, including the box geometry, ion types, box size, and ion concentration.
        - The ion names (`cation`, `anion`) should match the names used in your force field files (`ions.itp`) and the GROMACS library.
        - The function generates several intermediate files (e.g., `ions.tpr`, `em_steep.tpr`, `nvt.tpr`, `npt.tpr`, `md.tpr`) which are essential for the simulation workflow.
        - Ensure that the input files and directories exist and have the correct permissions before running the function.
    ```
    """
    fullpath = os.getcwd()
    # Initialize logging with log file in cwd
    log_file_path = os.path.join(fullpath, "process.log")
    # Suppress duplicate 'initialized' line if a handler already exists
    suppress = any(h for h in logging.getLogger().handlers if isinstance(h, logging.FileHandler))
    setup_logging(log_file_path, suppress_initial_message=suppress)
    logging.info(f"Logging initialized. Log file: {log_file_path}")

    # Exptract name of your input
    m_gro_name = strip_extension(m_gro)
    m_itp = f"{m_gro_name}.itp"
    m_top = f"{m_gro_name}.top"
    m_posres = f"posre_{m_gro_name.replace('_GMX', '')}.itp"
    mdp_dirname = os.path.basename(os.path.normpath(MDP_dir))

    if not os.path.isfile(m_itp):
        raise FileNotFoundError(f"The file {m_itp} does not exist.")
    if not os.path.isfile(m_top):
        raise FileNotFoundError(f"The file {m_top} does not exist.")
    if not os.path.isfile(m_posres):
        raise FileNotFoundError(f"The file {m_posres} does not exist.")

    # Replace "POSRES_LIG" with "POSRES" in topology
    replace_posres_in_file(m_top)

    # Add parameters for tip3p and ions to the topology file (include atomtypes.itp/tip3p.itp/ions.itp)
    include_ff(m_top, ff)
    include_tip3p(m_top, ff)
    include_ions(m_top, ff)
    remove_defaults_section_lines(m_top)

    # Setup scratch directory
    input_files = [
        m_gro,
        m_itp,
        m_top,
        m_posres,
        MDP_dir,
    ]  
    
    # Only include existing input files/directories
    scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch_dir)

    # Work inside a managed finalize context to avoid duplicate cleanup
    with ensure_finalized(
        original_dir=original_dir,
        scratch_dir=scratch_dir,
        input_files=input_files,
    ):
        os.chdir(scratch_dir)

        # Define dodecahedron box with ligand at center, > d nm to edge and rotate the system to the principal axes (-princ)
        run_command(
            f"echo 0 | gmx editconf -f {m_gro} -o {m_gro_name}_box.gro -bt {box_type} -d {d} -c -princ -nobackup",
            cwd=scratch_dir,
        )

        # Solvate the system
        run_command(
            f"gmx solvate -cp {m_gro_name}_box.gro -cs spc216 -o {m_gro_name}_tip3p.gro -p {m_gro_name}.top -nobackup",
            cwd=scratch_dir,
        )

        # Visualize changes made to topology
        run_command(f"tail {m_gro_name}.top", cwd=scratch_dir)

        # Create empty ions.mdp file (neccessary in order to generate ions.tpr file)
        run_command("touch ions.mdp", cwd=scratch_dir)

        # Add Ions
        run_command(
            f"gmx grompp -f ions.mdp -c {m_gro_name}_tip3p.gro -p {m_gro_name}.top -o ions.tpr",
            cwd=scratch_dir,
        )

        # Replace Solvent with Ions (i.e. NA, CL)
        run_command(
            f'echo "SOL" | gmx genion -s ions.tpr -o {m_gro_name}_tip3p_ions.gro -conc {conc} -p {m_gro_name}.top -pname {cation} -nname {anion} -neutral',
            cwd=scratch_dir,
        )

        # Visualize changes made to topology
        run_command(f"tail {m_gro_name}.top", cwd=scratch_dir)

        # Energy minimization
        run_command(
            f"gmx grompp -f {mdp_dirname}/em_steep.mdp -c {m_gro_name}_tip3p_ions.gro -p {m_gro_name}.top -o em_steep.tpr",
            cwd=scratch_dir,
        )
        run_command("gmx mdrun -deffnm em_steep -nobackup", cwd=scratch_dir)

        # Visualize energy convergence
        run_command(
            'echo "Potential\n0\n" | gmx energy -f em_steep.edr -o potential.xvg -xvg none',
            cwd=scratch_dir,
        )
        plot_svg("potential.xvg")

        # Helper to assemble runtime flags
        def _mdrun_flags():
            flags = []
            if threads is not None:
                flags += ["-nt", str(threads)]
            if pin is not None:
                flags += ["-pin", "on" if pin else "off"]
            flags.append("-nobackup")
            return " ".join(flags)

        # Equilibrate the system
        # NVT Equilibration
        run_command(
            f"gmx grompp -f {mdp_dirname}/NVT.mdp -c em_steep.gro -r em_steep.gro -p {m_gro_name}.top -o nvt.tpr -nobackup",
            cwd=scratch_dir,
        )
        run_command(f"gmx mdrun -deffnm nvt {_mdrun_flags()}", cwd=scratch_dir)
        run_command(
            'echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg -xvg none -b 20',
            cwd=scratch_dir,
        )
        plot_svg("temperature.xvg")

        # NPT Equilibration
        run_command(
            f"gmx grompp -f {mdp_dirname}/NPT.mdp -c nvt.gro -r nvt.gro -p {m_gro_name}.top -o npt.tpr -nobackup",
            cwd=scratch_dir,
        )
        run_command(f"gmx mdrun -deffnm npt {_mdrun_flags()}", cwd=scratch_dir)
        run_command(
            'echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg -xvg none',
            cwd=scratch_dir,
        )
        plot_svg("pressure.xvg")

        # Production: create input for production run
        run_command(
            f"gmx grompp -f {mdp_dirname}/Production.mdp -c npt.gro -t npt.cpt -p {m_gro_name}.top -o md.tpr",
            cwd=scratch_dir,
        )

        # Production run
        run_command(
            f"gmx mdrun -deffnm md {_mdrun_flags()}", cwd=scratch_dir
        )  # -v optional

        # Remove periodic jumps
        run_command(
            'echo 0 | gmx trjconv -s md.tpr -f md.xtc -o md_nojump.xtc -pbc nojump',
            cwd=scratch_dir,
        )

        # Center molecule
        run_command(
            'echo "1\n0\n" | gmx trjconv -s md.tpr -f md_nojump.xtc -o md_center.xtc -center -pbc mol',
            cwd=scratch_dir,
        )

        # Calculate distance between molecule and its periodic image
        run_command(
            'echo "1\n" | gmx mindist -s md.tpr -f md_center.xtc -pi -od mindist.xvg',
            cwd=scratch_dir,
        )

        # Radius of gyration (compactness)
        run_command(
            'echo "1" | gmx gyrate -f md_center.xtc -s md.tpr -o gyrate.xvg -xvg none',
            cwd=scratch_dir,
        )
        plot_svg("gyrate.xvg")

    if exit_screen:
        os.chdir(original_dir)
        subprocess.run(["screen", "-S", molecule_name, "-X", "quit"])
