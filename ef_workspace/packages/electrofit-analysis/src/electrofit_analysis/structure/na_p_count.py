import logging
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from electrofit.cli.run_commands import run_command
from electrofit.logging import setup_logging

PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
project_path = PROJECT_PATH

# Set the style for seaborn
sns.set_context("talk")


def create_index_file(structure_file, index_file):
    """
    Creates an index file with groups P, P1, P2, P3, P4, P5 using gmx make_ndx.
    """
    logging.info("Starting index file creation with gmx make_ndx...")

    # Define the commands to be sent to gmx make_ndx
    make_ndx_commands = "a P\na P1\na P2\na P3\na P4\na P5\nq\n"

    # Construct the shell command to pipe the commands into gmx make_ndx
    make_ndx_command = f'printf "{make_ndx_commands}" | gmx make_ndx -f {structure_file} -o {index_file}'
    run_command(make_ndx_command)
    logging.info("Index file created successfully.")


def run_pairdist_commands(
    trajectory, topology, index_file, groups, selection_group, output_prefix
):
    """
    Runs gmx pairdist for each group in groups against the selection_group.
    Produces .xvg files with min distance vs time.
    """
    logging.info("Starting pair distance calculations with gmx pairdist...")

    for i, group in enumerate(groups, start=1):
        output_file = f"{output_prefix}{i}.xvg"
        pairdist_command = (
            f"gmx pairdist -f {trajectory} -s {topology} -n {index_file} "
            f'-ref "{group}" -sel "{selection_group}" -pbc no -rmpbc yes -o {output_file}'
        )
        logging.info(f"Running gmx pairdist for group '{group}' -> '{output_file}'...")
        run_command(pairdist_command)
        logging.info(f"Distance calculation for group '{group}' completed.")

    logging.info("All gmx pairdist commands executed successfully.")


def plot_all_distances_subplots(
    output_prefix, num_groups, plot_filename="all_distances_subplots.pdf"
):
    """
    Plots all distance files (.xvg) as subplots (distance vs. time).
    """
    logging.info("Starting to plot all distance data as subplots...")
    # Temporarily suppress logging for cleaner output
    logging.disable(logging.CRITICAL)

    nrows, ncols = 2, 3  # 6 subplots in a 2x3 grid
    fig, axes = plt.subplots(nrows, ncols, figsize=(8, 6), sharex=True, sharey=True)
    axes = axes.flatten()

    for i in range(1, num_groups + 1):
        filename = f"{output_prefix}{i}.xvg"
        ax = axes[i - 1]
        try:
            data = np.loadtxt(filename, comments=("#", "@"))
            time = data[:, 0] / 1000.0  # Convert ps to ns
            distance = data[:, 1]
            group_label = f"P{i}"

            ax.plot(time, distance, color="black", linestyle="-", linewidth=0.5)
            mean_distance = np.mean(distance)
            # Add a red dashed line for the mean
            ax.axhline(y=mean_distance, color="red", linestyle="--", linewidth=1.5)
            ax.text(
                0.95,
                0.95,
                f"Mean: {mean_distance:.3f} nm",
                ha="right",
                va="top",
                transform=ax.transAxes,
                color="red",
                fontsize=12,
                bbox=dict(facecolor="white", alpha=0, edgecolor="none"),
            )

            ax.set_title(f"{group_label} - Na$^+$", fontsize=14)
            ax.set_xlabel("Time (ns)")
            if i % ncols == 1:
                ax.set_ylabel("Distance (nm)")
            ax.grid(False)
        except Exception:
            ax.text(
                0.5,
                0.5,
                "Error loading data",
                ha="center",
                va="center",
                transform=ax.transAxes,
                color="red",
            )
            ax.set_title(f"Distance: P{i} - NA (Error)", fontsize=14)

    plt.tight_layout()
    # fig.suptitle('Minimum Distances Between Phosphorus Groups and Na-Ions', fontsize=18, y=1.05)
    plt.savefig(plot_filename, dpi=300, bbox_inches="tight")
    # plt.show()  # Uncomment if you want to display
    logging.disable(logging.NOTSET)
    logging.info("Subplot plotting completed successfully.")


# --------------------------------------------------------------------
# NEW FEATURE: Counting how many ions are within a certain cutoff
# --------------------------------------------------------------------


def run_ion_count_commands(
    trajectory,
    topology,
    index_file,
    groups,
    selection_group,
    cutoff=0.5,
    output_prefix="ion_count_",
):
    """
    Uses gmx select to count how many ions (selection_group) are within a given
    distance cutoff of each phosphate group. Produces .xvg files with counts over time.
    """
    logging.info(f"Starting ion count within {cutoff} nm for each group in {groups}...")

    for group in groups:
        # We'll name each output file based on the group, e.g. ion_count_P.xvg, ion_count_P1.xvg, etc.
        output_file = f"{output_prefix}{group}.xvg"
        # Construct the 'gmx select' expression
        select_expr = (
            f'group "{selection_group}" and within {cutoff} of group "{group}"'
        )
        command = (
            f"gmx select -f {trajectory} -s {topology} -n {index_file} "
            f"-select '{select_expr}' "
            f"-os {output_file} "  # .xvg with # of selected atoms vs time
            f"-on dummy_index.ndx "  # optional index file (not used further, but needed by gmx)
            f"2>&1"
        )
        logging.info(
            f"Running gmx select for group '{group}', output -> '{output_file}'..."
        )
        run_command(command)

    logging.info("All ion count commands executed successfully.")


# --------------------------------------------------------------------
# NEW FEATURE (volume grid): buried IP6 volume per frame via gmx sasa
# --------------------------------------------------------------------


def run_buried_volume_commands(
    trajectory,
    topology,
    index_file,
    groups,
    ip6_group="Other",
    cutoff=0.5,
    output_prefix="buriedVol_",
    ndots=384,
):
    """
    For every phosphate group make GROMACS evaluate, per frame, how much
    van‑der‑Waals volume of IP6 atoms falls inside radius *cutoff*.

    Implementation:
      gmx sasa -probe 0  -surface '<dynamic selection>' -tv <out.xvg>

    The dynamic selection is:
        group "{ip6_group}" and within {cutoff} of group "{P_i}"
    which is re‑evaluated each frame by the selection engine.

    Output: <output_prefix><group>.xvg, two columns (time [ps], volume [nm^3])
    """
    cmd_common = [
        "gmx",
        "sasa",
        "-f",
        trajectory,
        "-s",
        topology,
        "-n",
        index_file,
        "-probe",
        "0",
        "-ndots",
        str(ndots),  # ↑ accuracy
        "-quiet",  # keep logs succinct
    ]

    logging.info("Starting buried-volume calculations…")
    for grp in groups:
        outfile = f"{output_prefix}{grp}.xvg"
        sel = f'group "{ip6_group}" and within {cutoff} of group "{grp}"'
        cmd = cmd_common + ["-surface", sel, "-tv", outfile]
        logging.info(f"  {outfile}  ←  {sel}")
        run_command(cmd)
    logging.info("Finished.")


def plot_ion_counts_subplots(
    output_prefix, groups, plot_filename="ion_counts_subplots.pdf"
):
    """
    Plots the number of ions within a cutoff for each group as subplots (count vs. time).
    """
    logging.info("Starting to plot ion counts as subplots...")
    logging.disable(logging.CRITICAL)

    num_groups = len(groups)
    nrows, ncols = 2, 3  # Adjust for 6 groups
    fig, axes = plt.subplots(nrows, ncols, figsize=(8, 6), sharex=True, sharey=True)
    axes = axes.flatten()

    for i, group in enumerate(groups):
        filename = f"{output_prefix}{group}.xvg"
        ax = axes[i]
        try:
            data = np.loadtxt(filename, comments=("#", "@"))
            time = data[:, 0] / 1000.0  # time in ns if .xvg is in ps
            ion_count = data[:, 1]
            mean_count = np.mean(ion_count)

            ax.plot(time, ion_count, color="darkblue", linestyle="-", linewidth=0.5)
            ax.axhline(y=mean_count, color="red", linestyle="--", linewidth=1.5)
            ax.text(
                0.95,
                0.95,
                f"Mean: {mean_count:.1f}",
                ha="right",
                va="top",
                transform=ax.transAxes,
                color="red",
                fontsize=12,
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
            )

            ax.set_title(f"{group}", fontsize=14)
            ax.set_xlabel("Time (ns)")
            if i % ncols == 0:
                ax.set_ylabel("Na+ ions")
            ax.grid(False)
        except Exception:
            ax.text(
                0.5,
                0.5,
                "Error loading data",
                ha="center",
                va="center",
                transform=ax.transAxes,
                color="red",
            )
            ax.set_title(f"{group} (Error)", fontsize=14)

    # If we have fewer than nrows*ncols subplots, hide the empty ones
    for j in range(num_groups, nrows * ncols):
        fig.delaxes(axes[j])

    plt.tight_layout()
    # fig.suptitle('Number of Na+ Ions Within Cutoff Distance', fontsize=18, y=1.05)
    plt.savefig(plot_filename, dpi=300, bbox_inches="tight")
    # plt.show()
    logging.disable(logging.NOTSET)
    logging.info("Ion count subplot plotting completed successfully.")


# --------------------------------------------------------------------
#  Third-order (frame-specific) excess definition
#
#  N_excess(t) = N_in_sphere(t) / [ c · ( 4/3·π·r³  −  V_IP6(t) ) ]
#
#  c           = N_total_Na / V_box  (bulk concentration)
#  V_IP6(t)    = n_IP6_atoms_in_sphere(t) · v_atom_IP6  (per-atom excluded volume)
#  r           = cutoff_radius used in the gmx select commands.
#
#  The code logs every piece of the formula so you can reconstruct how
#  each number was obtained when reading the log file.


def _parse_box_volume_from_gro(gro_file):
    """
    Returns simulation box volume (nm³) from the last line of a .gro file.
    For orthorhombic boxes the last line contains three floats lx ly lz.
    For triclinic boxes the first three floats are the box vector lengths
    (off‑diagonal terms are ignored → gives upper bound to actual volume).
    """
    with open(gro_file, "r") as fh:
        lines = fh.readlines()
        # last non‑empty line should contain box vectors
        box_line = lines[-1].strip().split()
        if len(box_line) < 3:
            raise ValueError("Could not read box vectors from .gro file.")
        lx, ly, lz = map(float, box_line[:3])
        return lx * ly * lz  # nm³


def run_ip6_exclusion_count_commands(
    trajectory,
    topology,
    index_file,
    groups,
    ip6_group="Other",
    cutoff=0.5,
    output_prefix="ip6_count_",
):
    """
    Uses gmx select to count how many IP6 atoms (or whatever group name is
    passed via `ip6_group`) are found within `cutoff` nm of each phosphate
    group *for every frame* of the trajectory.

    The output is written to <output_prefix><group>.xvg files in exactly the
    same format that `run_ion_count_commands` produces, enabling frame‑wise
    combination with the Na⁺ counts.

    Parameters
    ----------
    trajectory : str
        Path to the .xtc/.trr trajectory file.
    topology : str
        Path to the .tpr (or any structure readable by GROMACS).
    index_file : str
        .ndx file containing both the phosphate groups and the IP6 group.
    groups : list[str]
        List of phosphate group names (e.g. ["P", "P1", ...]).
    ip6_group : str, optional
        Name of the index group that represents all IP6 atoms. Default "Other".
    cutoff : float, optional
        Cut‑off radius in nm (must be identical to the one used for the Na⁺
        counts). Default 0.5 nm.
    output_prefix : str, optional
        Prefix for the generated .xvg files. Default "ip6_count_".
    """
    logging.info(
        f"Starting IP6 exclusion count within {cutoff} nm for groups: {groups}"
    )
    for group in groups:
        out_file = f"{output_prefix}{group}.xvg"
        select_expr = f'group "{ip6_group}" and within {cutoff} of group "{group}"'
        cmd = (
            f"gmx select -f {trajectory} -s {topology} -n {index_file} "
            f"-select '{select_expr}' "
            f"-os {out_file} "  # writes time‑series of #selected atoms
            f"-on dummy_ip6_index.ndx "  # dummy index output (ignored later)
            f"2>&1"
        )
        logging.info(f'gmx select (IP6) for group "{group}" → "{out_file}"')
        run_command(cmd)
    logging.info("IP6 exclusion counting completed for all groups.")


def _count_na_atoms_in_gro(gro_file):
    """
    Counts how many atoms have the string 'NA' in columns 11‑15 (atom name)
    of a .gro file. Works for standard GROMACS naming.
    """
    count = 0
    with open(gro_file, "r") as fh:
        lines = fh.readlines()
    # skip first 2 title/atom‑number lines, last line is box
    for line in lines[2:-1]:
        atom_name = line[10:15].strip()
        if atom_name.startswith("NA"):
            count += 1
    if count == 0:
        raise ValueError("No Na⁺ atoms found in structure file; check naming.")
    return count


def compute_excess_ion_counts(
    ion_count_prefix,
    buried_prefix,
    groups,
    structure_file,
    cutoff_radius,
    output_prefix="excess_",
):
    """
    Excess‑ion factor using *true* solvent‑accessible volume per frame:

        N_excess(t) = N_in_sphere(t) /
                      [ c_Na · (4/3·π·r³ − V_buried(t)) ]

    where V_buried(t) is obtained from gmx sasa (-tv) and already has units nm³.
    """
    import math

    # Bulk Na⁺ concentration
    n_na_total = _count_na_atoms_in_gro(structure_file)
    v_box = _parse_box_volume_from_gro(structure_file)
    c_na = n_na_total / v_box

    v_sphere = (4.0 / 3.0) * math.pi * cutoff_radius**3
    logging.info("============================================================")
    logging.info("EXCESS‑ION FORMULA (grid/SAV version):")
    logging.info("  N_excess(t) = N_in / [ c_Na · (4/3·π·r³ − V_buried(t)) ]")
    logging.info(f"  r = {cutoff_radius:.3f} nm ; 4/3·π·r³ = {v_sphere:.6f} nm³")
    logging.info(
        f"  c_Na = N_Na_tot / V_Box = {n_na_total}/{v_box:.6f} = {c_na:.6f} nm⁻³"
    )
    logging.info("============================================================")

    for group in groups:
        na_file = f"{ion_count_prefix}{group}.xvg"
        buried_file = f"{buried_prefix}{group}.xvg"
        out_file = f"{output_prefix}{group}.xvg"
        try:
            na_data = np.loadtxt(na_file, comments=("#", "@"))
            buried_data = np.loadtxt(buried_file, comments=("#", "@"))
            if na_data.shape[0] != buried_data.shape[0]:
                raise ValueError("Frame mismatch between Na and buried‑volume files")

            time_ps = na_data[:, 0]
            n_in = na_data[:, 1]
            v_buried = buried_data[:, 1]
            v_eff = v_sphere - v_buried
            # ensure positive effective volume
            v_eff[v_eff <= 1e-9] = 1e-9

            n_solution = c_na * v_eff
            n_excess = n_in / n_solution

            header = (
                '@    title "Excess Na+ ions (grid SAV)"\n'
                '@    xaxis  label "Time (ps)"\n'
                '@    yaxis  label "N_excess"\n'
            )
            np.savetxt(
                out_file,
                np.column_stack((time_ps, n_excess)),
                header=header,
                comments="",
            )

            # log means
            logging.info(
                f"[{group}] mean N_in = {np.mean(n_in):.4f}, "
                f"mean V_buried = {np.mean(v_buried):.4f} nm³, "
                f"mean V_eff = {np.mean(v_eff):.4f} nm³, "
                f"mean N_excess = {np.mean(n_excess):.4f}"
            )
        except Exception as e:
            logging.error(f"Failed for {group}: {e}")


def plot_excess_ion_counts_subplots(
    output_prefix, groups, plot_filename="excess_ion_counts_subplots.pdf"
):
    """
    Generates subplots (N_excess vs time) for each group.
    """
    logging.info("Plotting excess ion counts as subplots...")
    logging.disable(logging.CRITICAL)

    num_groups = len(groups)
    nrows, ncols = 2, 3
    fig, axes = plt.subplots(nrows, ncols, figsize=(8, 6), sharex=True, sharey=True)
    axes = axes.flatten()

    for i, group in enumerate(groups):
        filename = f"{output_prefix}{group}.xvg"
        ax = axes[i]
        try:
            data = np.loadtxt(filename, comments=("#", "@"))
            time = data[:, 0] / 1000.0  # ps → ns
            n_excess = data[:, 1]
            mean_excess = np.mean(n_excess)

            ax.plot(time, n_excess, linestyle="-", linewidth=0.5)
            ax.axhline(y=mean_excess, linestyle="--", linewidth=1.5)
            ax.text(
                0.95,
                0.95,
                f"Mean: {mean_excess:.2f}",
                ha="right",
                va="top",
                transform=ax.transAxes,
                fontsize=12,
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
            )

            ax.set_title(f"{group}", fontsize=14)
            ax.set_xlabel("Time (ns)")
            if i % ncols == 0:
                ax.set_ylabel("N_excess")
            ax.grid(False)
        except Exception:
            ax.text(
                0.5,
                0.5,
                "Error loading data",
                ha="center",
                va="center",
                transform=ax.transAxes,
                color="red",
            )
            ax.set_title(f"{group} (Error)", fontsize=14)

    # hide unused axes
    for j in range(num_groups, nrows * ncols):
        fig.delaxes(axes[j])

    plt.tight_layout()
    plt.savefig(plot_filename, dpi=300, bbox_inches="tight")
    # plt.show()
    logging.disable(logging.NOTSET)
    logging.info("Excess ion count subplot plotting completed successfully.")


# Optional: Count how many ions are near the *entire molecule* (e.g., 'MOL' group)
def run_ion_count_whole_molecule(
    trajectory,
    topology,
    index_file,
    whole_group="Other",
    selection_group="NA",
    cutoff=0.5,
    output_xvg="ion_count_MOL.xvg",
):
    """
    If you have an index group for the entire molecule (named e.g. 'MOL'),
    this function will count how many 'selection_group' atoms are within
    `cutoff` nm of that entire molecule.
    """
    logging.info(f"Starting ion count for entire molecule group '{whole_group}'...")
    select_expr = (
        f'group "{selection_group}" and within {cutoff} of group "{whole_group}"'
    )
    command = (
        f"gmx select -f {trajectory} -s {topology} -n {index_file} "
        f"-select '{select_expr}' "
        f"-os {output_xvg} "
        f"-on dummy_mol_index.ndx "
        f"2>&1"
    )
    run_command(command)
    logging.info("Ion count for entire molecule completed successfully.")


def plot_whole_molecule_ion_count(xvg_file, plot_filename="ion_count_whole_mol.pdf"):
    """
    Plots a single line: # of NA ions within cutoff for the entire molecule vs time.
    """
    logging.info("Plotting ion count for the entire molecule...")
    try:
        data = np.loadtxt(xvg_file, comments=("#", "@"))
        time = data[:, 0] / 1000.0
        ion_count = data[:, 1]
        mean_count = np.mean(ion_count)

        plt.figure(figsize=(6, 4))
        plt.plot(time, ion_count, color="darkblue", linestyle="-", linewidth=0.5)
        plt.axhline(y=mean_count, color="red", linestyle="--", linewidth=1.5)
        plt.text(
            0.95,
            0.1,
            f"Mean: {mean_count:.1f}",
            ha="right",
            va="top",
            transform=plt.gca().transAxes,
            color="red",
            fontsize=12,
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"),
        )
        # plt.title('Na+ Ion Count Near Entire Molecule')
        plt.xlabel("Time (ns)")
        plt.ylabel("# of Na+ ions")
        plt.tight_layout()
        plt.savefig(plot_filename, dpi=300)
        # plt.show()
        logging.info("Ion count for entire molecule plotted successfully.")
    except Exception as e:
        logging.error(f"Error plotting whole molecule ion count: {e}")


# ----------------------------
# Main Execution Flow
# ----------------------------
def main():
    # Define the base process directory
    process_dir = os.path.join(project_path, "process")

    # Loop through each subdirectory in the process directory
    for folder_name in os.listdir(process_dir):
        folder_path = os.path.join(process_dir, folder_name)

        if os.path.isdir(folder_path):
            # Define the 'run_final_gmx_simulation' directory within this folder
            run_final_sim_dir = os.path.join(folder_path, "run_final_gmx_simulation")

            if os.path.isdir(run_final_sim_dir):
                # Define the destination directory 'analyze_final_sim'
                dest_dir = os.path.join(folder_path, "analyze_final_sim")
                dest_dir = os.path.join(dest_dir, "NaP_dist_count")

                os.makedirs(dest_dir, exist_ok=True)
                os.chdir(dest_dir)

                # Define paths (adjust these as necessary)
                structure_file = os.path.join(run_final_sim_dir, "md.gro")
                trajectory_file = os.path.join(run_final_sim_dir, "md_center.xtc")
                topology_file = os.path.join(run_final_sim_dir, "md.tpr")
                index_file = "NA_P_index.ndx"  # We'll create this file
                selection_group = "NA"
                # We'll produce distances like distances_NaP1.xvg -> distances_NaP6.xvg
                output_prefix = "distances_NaP"
                log_file = "distances_NaP_gmx.log"

                # List of phosphorus groups
                p_groups = ["P", "P1", "P2", "P3", "P4", "P5"]

                # Setup logging
                setup_logging(log_file)
                logging.info("Logging is set up.")

                # Check if input files exist
                input_files = [structure_file, trajectory_file, topology_file]
                for file in input_files:
                    if not os.path.isfile(file):
                        logging.error(
                            f"Required input file '{file}' not found. Exiting."
                        )
                        sys.exit(1)
                logging.info("All required input files are present.")

                # Step 1: Create the index file (includes P, P1, P2, P3, P4, P5)
                create_index_file(structure_file, index_file)

                # Step 2: Run gmx pairdist for each phosphorus group
                run_pairdist_commands(
                    trajectory=trajectory_file,
                    topology=topology_file,
                    index_file=index_file,
                    groups=p_groups,
                    selection_group=selection_group,
                    output_prefix=output_prefix,
                )

                # Step 3: Plot the distance data
                plot_all_distances_subplots(
                    output_prefix=output_prefix,
                    num_groups=len(p_groups),
                    plot_filename="all_distances_subplots.pdf",
                )

                # NEW STEPS: Count how many ions are close to each phosphate group
                # Adjust cutoff as desired, e.g., 0.5 nm
                cutoff_distance = 0.5
                ion_count_prefix = "ion_count_"

                run_ion_count_commands(
                    trajectory=trajectory_file,
                    topology=topology_file,
                    index_file=index_file,
                    groups=p_groups,
                    selection_group=selection_group,
                    cutoff=cutoff_distance,
                    output_prefix=ion_count_prefix,
                )

                # Plot the ion count data for each group
                plot_ion_counts_subplots(
                    output_prefix=ion_count_prefix,
                    groups=p_groups,
                    plot_filename="ion_counts_subplots.pdf",
                )

                # ----------------------------------------------------------------
                # NEW FEATURE 2: compute & plot excess Na⁺ counts
                # ----------------------------------------------------------------
                excess_prefix = "excess_NaP_"
                buried_prefix = "buriedVol_"
                run_buried_volume_commands(
                    trajectory=trajectory_file,
                    topology=topology_file,
                    index_file=index_file,
                    groups=p_groups,
                    ip6_group="Other",
                    cutoff=cutoff_distance,
                    output_prefix=buried_prefix,
                )
                compute_excess_ion_counts(
                    ion_count_prefix=ion_count_prefix,
                    buried_prefix=buried_prefix,
                    groups=p_groups,
                    structure_file=structure_file,
                    cutoff_radius=cutoff_distance,
                    output_prefix=excess_prefix,
                )
                plot_excess_ion_counts_subplots(
                    output_prefix=excess_prefix,
                    groups=p_groups,
                    plot_filename="excess_ion_counts_subplots.pdf",
                )

                # OPTIONAL: If you have an index group for the entire molecule
                # (e.g., "MOL") in your index file, you could do:
                run_ion_count_whole_molecule(
                    trajectory=trajectory_file,
                    topology=topology_file,
                    index_file=index_file,
                    whole_group="Other",  # must exist in your index
                    selection_group="NA",
                    cutoff=cutoff_distance,
                    output_xvg="ion_count_MOL.xvg",
                )

                plot_whole_molecule_ion_count(
                    xvg_file="ion_count_MOL.xvg",
                    plot_filename="ion_count_whole_mol.pdf",
                )

                logging.info(
                    "All tasks (distance + ion count) completed successfully.\n"
                )


if __name__ == "__main__":
    main()
