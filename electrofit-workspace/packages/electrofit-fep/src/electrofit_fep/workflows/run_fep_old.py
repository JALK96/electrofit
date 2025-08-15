#!/usr/bin/env python3
import logging
import os
import shutil
import subprocess

import matplotlib.pyplot as plt
import pandas as pd

# Import helper functions from your electrofit modules
from electrofit.cli.run_commands import run_command
from electrofit.io.files import (
    strip_extension,
)
from electrofit.logging import setup_logging
from electrofit.scratch.manager import (
    finalize_scratch_directory,
    setup_scratch_directory,
)
from electrofit.cli.safe_run import register_scratch


def plot_svg(svg):
    """Reads an .xvg file and plots it as a PDF."""
    # Disable logging while plotting
    logging.disable(logging.CRITICAL)
    name = strip_extension(svg)
    df = pd.read_csv(svg, sep="\s+", header=None, names=["time", name])
    plt.figure(figsize=(8, 6))
    plt.plot(df["time"], df[name], color="darkblue", linestyle="-", label=name)
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel(name, fontsize=14)
    plt.axhline(0, color="black", linewidth=0.8, linestyle="--")
    plt.tight_layout()
    plt.savefig(f"{name}.pdf", format="pdf")
    plt.close()
    logging.disable(logging.NOTSET)


def set_up_fep_production(
    m_gro,
    MDP_dir,
    base_scratch_dir,
    molecule_name,
    box_type="dodecahedron",
    cation="NA",
    anion="CL",
    d="1.2",
    conc="0.15",
    lambda_values=None,
    exit_screen=True,
):
    """
    Set up and execute a free‐energy perturbation (FEP) simulation using GROMACS.

    The workflow includes:
      1. Preparing the topology by replacing position restraint tags,
         and including force‐field, TIP3P water and ion parameters.
      2. Setting up a scratch directory and defining a simulation box (centered,
         with a minimum distance "d" to the box edge).
      3. Solvating the system and adding ions (to neutralize and achieve the specified concentration).
      4. Energy minimization, followed by NVT and NPT equilibration.
      5. Running production FEP simulations over a series of lambda windows using a template MDP.
      6. Analyzing the free‐energy differences via BAR.
      7. Finalizing the scratch directory and (optionally) closing a screen session.

    Args:
        m_gro (str): Path to the input .gro file.
        MDP_dir (str): Directory containing the MDP files (e.g. em_steep.mdp, NVT.mdp, NPT.mdp, and FEP_template.mdp).
        base_scratch_dir (str): Base directory for the temporary scratch space.
        molecule_name (str): A name for the molecule/simulation (used for screen session and file naming).
        box_type (str, optional): Type of simulation box (default: "dodecahedron").
        cation (str, optional): Name of the cation (default: "NA").
        anion (str, optional): Name of the anion (default: "CL").
        d (str or float, optional): Distance (nm) from the molecule to the box edge (default: "1.2").
        conc (str or float, optional): Ion concentration in mol/L (default: "0.15").
        lambda_values (list, optional): List of lambda values for the FEP simulation.
            If None, defaults to [0.0, 0.1, ..., 1.0].
        exit_screen (bool, optional): Whether to quit the screen session at the end (default: True).

    Raises:
        FileNotFoundError: If any of the required input files are missing.
    """
    if lambda_values is None:
        lambda_values = list(range(13))

    # Set up logging
    fullpath = os.getcwd()
    log_file_path = os.path.join(fullpath, "fep_process.log")
    setup_logging(log_file_path)
    logging.info(f"Logging initialized. Log file: {log_file_path}")

    # Derive base file names from the .gro file
    m_gro_name = strip_extension(m_gro)
    m_itp = f"{m_gro_name}.itp"
    m_top = f"{m_gro_name}.top"
    m_posres = f"posre_{m_gro_name}.itp"

    for f in [m_itp, m_top, m_posres]:
        if not os.path.isfile(f):
            raise FileNotFoundError(f"Required file {f} not found.")

    # Modify the topology file: replace POSRES_LIG with POSRES, include force field, TIP3P, and ion files.
    # replace_posres_in_file(m_top)
    # include_ff(m_top)
    # include_tip3p(m_top)
    # include_ions(m_top)
    # remove_defaults_section_lines(m_top)

    # Set up scratch directory and copy input files
    input_files = [m_gro, m_itp, m_top, m_posres, MDP_dir]
    scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch_dir)
    register_scratch(
    original_dir=original_dir,
    scratch_dir=scratch_dir,
    input_files=input_files,
    # output_files=None,               # optional (default: copy everything else)
    # overwrite=True,                  # optional
    # remove_parent_if_empty=True,     # optional 
)
    os.chdir(scratch_dir)
    logging.info(f"Changed working directory to scratch directory: {scratch_dir}")

    try:
        # 1. Define simulation box (center molecule, set box type and distance)
        run_command(
            f"echo 0 | gmx editconf -f {m_gro} -o {m_gro_name}_box.gro -bt {box_type} -d {d} -c -princ -nobackup",
            cwd=scratch_dir,
        )

        # 2. Solvate the system using TIP3P water
        run_command(
            f"gmx solvate -cp {m_gro_name}_box.gro -cs spc216 -o {m_gro_name}_tip3p.gro -p {m_top} -nobackup",
            cwd=scratch_dir,
        )
        run_command(f"tail {m_top}", cwd=scratch_dir)

        # 3. Generate ions: create an empty ions.mdp and then run grompp and genion
        run_command("touch ions.mdp", cwd=scratch_dir)
        run_command(
            f"gmx grompp -f ions.mdp -c {m_gro_name}_tip3p.gro -p {m_top} -o ions.tpr",
            cwd=scratch_dir,
        )
        run_command(
            f'echo "SOL" | gmx genion -s ions.tpr -o {m_gro_name}_tip3p_ions.gro -conc {conc} -p {m_top} -pname {cation} -nname {anion} -neutral',
            cwd=scratch_dir,
        )
        run_command(f"tail {m_top}", cwd=scratch_dir)

        # 4. Energy Minimization
        run_command(
            f"gmx grompp -f {os.path.join(MDP_dir, 'em_steep.mdp')} -c {m_gro_name}_tip3p_ions.gro -p {m_top} -o em.tpr -nobackup",
            cwd=scratch_dir,
        )
        run_command("gmx mdrun -deffnm em -nobackup", cwd=scratch_dir)
        run_command(
            'echo "Potential\n0\n" | gmx energy -f em.edr -o potential.xvg -xvg none',
            cwd=scratch_dir,
        )
        plot_svg("potential.xvg")

        # 5. NVT Equilibration
        run_command(
            f"gmx grompp -f {os.path.join(MDP_dir, 'NVT.mdp')} -c em.gro -r em.gro -p {m_top} -o nvt.tpr -nobackup",
            cwd=scratch_dir,
        )
        run_command("gmx mdrun -deffnm nvt -nt 16 -pin on -nobackup", cwd=scratch_dir)
        run_command(
            'echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg -xvg none -b 20',
            cwd=scratch_dir,
        )
        plot_svg("temperature.xvg")

        # 6. NPT Equilibration
        run_command(
            f"gmx grompp -f {os.path.join(MDP_dir, 'NPT.mdp')} -c nvt.gro -r nvt.gro -p {m_top} -o npt.tpr -nobackup",
            cwd=scratch_dir,
        )
        run_command("gmx mdrun -deffnm npt -nt 16 -pin on -nobackup", cwd=scratch_dir)
        run_command(
            'echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg -xvg none',
            cwd=scratch_dir,
        )
        plot_svg("pressure.xvg")

        # 7. Production FEP Runs with separate lambda directories:
        for lam in lambda_values:
            lam_dir = os.path.join(scratch_dir, f"lambda_{lam:0>2}")
            os.mkdir(lam_dir)
            shutil.copy(
                os.path.join(scratch_dir, "npt.gro"), os.path.join(lam_dir, "conf.gro")
            )
            shutil.copy(m_top, lam_dir)
            logging.info(
                f"Running FEP simulation at lambda = {lam} in directory {lam_dir}"
            )
            run_command(
                f"sed 's/REPLACE_LAMBDA/{lam}/g' {os.path.join(MDP_dir, 'FEP_template.mdp')} > fep_{lam}.mdp",
                cwd=lam_dir,
            )
            run_command(
                f"gmx grompp -f fep_{lam}.mdp -c conf.gro -p {os.path.basename(m_top)} -o fep_{lam}.tpr -nobackup",
                cwd=lam_dir,
            )
            run_command(
                f"gmx mdrun -deffnm fep_{lam} -nobackup",
                cwd=lam_dir,
            )

        # Copy xvg files back to the scratch directory
        for lam in lambda_values:
            lam_dir = os.path.join(scratch_dir, f"lambda_{lam:0>2}")
            xvg_file = os.path.join(lam_dir, f"fep_{lam}.xvg")
            if os.path.exists(xvg_file):
                shutil.copy(xvg_file, scratch_dir)

        # Run BAR analysis in the scratch directory
        run_command(
            "gmx bar -f fep_*.xvg -o dG.xvg -error dG_error.xvg",
            cwd=scratch_dir,
        )

    except Exception as err:
        logging.error(f"Error occurred during simulation: {err}")
        logging.info("Finalizing scratch directory due to error.")
        finalize_scratch_directory(original_dir, scratch_dir, input_files)
        raise  # Optionally re-raise the error for further handling

    # 8. Finalize the scratch directory (if we reach this point without errors)
    finalize_scratch_directory(original_dir, scratch_dir, input_files)

    if exit_screen:
        os.chdir(original_dir)
        subprocess.run(["screen", "-S", molecule_name, "-X", "quit"])


if __name__ == "__main__":
    # Example call to set_up_fep_production; adjust parameters and paths as needed.
    set_up_fep_production(
        m_gro="IP_111011_111101.gro",
        MDP_dir="MDP",
        base_scratch_dir="/scratch/johannal96/tmp/",
        molecule_name="IP_111011_111101",
        box_type="dodecahedron",
        cation="NA",
        anion="CL",
        d="1.2",
        conc="0.15",
        lambda_values=list(range(13)),
        exit_screen=True,
    )
