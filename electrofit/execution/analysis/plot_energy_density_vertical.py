#!/usr/bin/env python3
"""
Script: plot_energy_density_vertical.py

This script extracts energy data from GROMACS binary energy files (ener.edr) located
in lambda window subdirectories. It can handle either one or two base directories:

  - If only one directory is provided, it plots a single (top) energy density subplot.
  - If two directories are provided, it plots two vertical subplots that share the x-axis:

     1) Top subplot from the first base directory.
     2) Bottom subplot from the second base directory with the y-axis inverted.

Instead of using a legend for each lambda window, a colorbar represents the lambda index
(or window index) from 0→1 (top) and 1→0 (bottom).

Usage:
    python plot_energy_density_vertical.py <base_dir1> [<base_dir2>] [--energy_term <int>]

Arguments:
    base_dir1      Path to the first base directory containing "lambda_*" subdirectories.
    base_dir2      (Optional) Path to the second base directory containing "lambda_*" subdirectories.
    --energy_term  (Optional) Index of the energy term to extract (default: 10).

Requirements:
    - GROMACS must be installed and accessible (for gmx energy).
    - Python packages: numpy, matplotlib, seaborn.
    - Each base directory should contain subdirectories named "lambda_0", "lambda_1", …,
      each containing an "ener.edr" file.
"""

import os
import sys
import subprocess
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context("talk", font_scale=1.2)

def extract_energy_data(ener_edr, output_xvg, energy_term=10):
    """
    Extracts a specific energy term from a binary ener.edr file using gmx energy,
    and writes the result to an .xvg file.
    """
    cmd = f'echo "{energy_term}\n0" | gmx energy -f {ener_edr} -o {output_xvg} -quiet'
    subprocess.run(cmd, shell=True, check=True)

def parse_xvg(xvg_file):
    """
    Parses a GROMACS .xvg file and returns its numeric data as a NumPy array.
    Lines starting with '#' or '@' are skipped.
    """
    data = []
    with open(xvg_file, 'r') as f:
        for line in f:
            if line.startswith('#') or line.startswith('@'):
                continue
            parts = line.strip().split()
            if parts:
                data.append([float(x) for x in parts])
    return np.array(data)

def plot_single_subdir_set(base_dir, energy_term=10, invert_y=False, ax=None, title="", units="kJ"):
    """
    Plots KDE energy distributions for all lambda_* directories in a single base_dir.

    Parameters
    ----------
    base_dir : str
        Directory containing lambda_* subdirectories.
    energy_term : int
        The energy term index to extract (default: 10).
    invert_y : bool
        Whether to invert the y-axis for the plot (for the second subplot).
    ax : matplotlib.axes.Axes
        The axes on which to plot. Required so we can manage subplots externally.
    title : str
        The title for this subplot.
    """
    # Get sorted list of lambda directories
    lambda_dirs = sorted(
        [os.path.join(base_dir, d) for d in os.listdir(base_dir) if d.startswith("lambda_")],
        key=lambda d: float(d.split("_")[-1])
    )
    n_lam = len(lambda_dirs)
    if n_lam == 0:
        print(f"No lambda_* directories found in {base_dir}.")
        return

    # Use a consistent colormap
    cmap = plt.get_cmap("viridis")

    for idx, lam_dir in enumerate(lambda_dirs):
        ener_edr = os.path.join(lam_dir, "ener.edr")
        output_xvg = os.path.join(lam_dir, "energy.xvg")
        extract_energy_data(ener_edr, output_xvg, energy_term=energy_term)
        data = parse_xvg(output_xvg)
        if data.size == 0:
            print(f"No data found in {output_xvg}. Skipping.")
            continue

        energy = data[:, 1]
        if units.lower() == "kcal":
            energy *= 0.239  # convert energy from kJ/mol to kcal/mol

        # Map index to a 0..1 color value
        if not invert_y:
            color_val = idx / (n_lam - 1) if n_lam > 1 else 0
        else:
            # Invert the color scale from 1..0
            color_val = 1.0 - (idx / (n_lam - 1)) if n_lam > 1 else 1.0

        color = cmap(color_val)
        sns.kdeplot(energy, ax=ax, color=color, lw=1.5)

    ax.set_title(title)
    ax.set_ylabel(r"$P_{fwd}(E)$")
    if invert_y:
        ax.invert_yaxis()

def plot_energy_density_vertical(base_dir1, base_dir2=None, energy_term=10, units="kJ"):
    """
    Plots energy density distributions from 1 or 2 directories:
      - If only base_dir1 is provided, plots a single subplot.
      - If base_dir2 is also provided, plots two vertical subplots.
        The top subplot is from base_dir1, the bottom from base_dir2 (with inverted y-axis).

    Instead of legends for each lambda window, we include a colorbar indicating
    the lambda index range from 0 to 1 (top) or 1 to 0 (bottom).
    """
    # If only one directory is given: single subplot
    if base_dir2 is None:
        fig, ax_top = plt.subplots(1, 1, figsize=(10, 5))
        plot_single_subdir_set(
            base_dir=base_dir1, 
            energy_term=energy_term, 
            invert_y=False, 
            ax=ax_top, 
            title="A → B"
        )
        ax_top.set_xlabel(f"Energy ({units}/mol)")
        
        # Create a colorbar from 0→1
        norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax_top, orientation="vertical", fraction=0.05)
        cbar.set_label(r"$\lambda$")
        plt.tight_layout()
        plt.savefig("plot_single_E.pdf")

    else:
        # Two directories: 2 vertical subplots
        fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        # Top subplot: 0→1 color scale
        plot_single_subdir_set(
            base_dir=base_dir1, 
            energy_term=energy_term, 
            invert_y=False, 
            ax=ax_top,
            title="A → B"
        )
        ax_top.set_ylabel(r"$P_{fwd}(E)$")

        # Bottom subplot: 1→0 color scale
        plot_single_subdir_set(
            base_dir=base_dir2, 
            energy_term=energy_term, 
            invert_y=True, 
            ax=ax_bottom,
            title="B → A"
        )
        ax_bottom.set_ylabel(r"$P_{bwd}(E)$")
        ax_bottom.set_xlabel(f"Energy ({units}/mol)")

        # Create a colorbar for each subplot that indicates the indexing:
        #   top: 0→1, bottom: 1→0
        # We can do it with two colorbars or a single colorbar. 
        # For clarity, let's do a single colorbar showing 0→1:

        norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=[ax_top, ax_bottom], orientation="vertical", fraction=0.05)
        cbar.set_label(r"$\lambda$")

        #plt.tight_layout()
        plt.savefig("plot_double_E.pdf")

def main():
    parser = argparse.ArgumentParser(
        description="Plot vertical energy density distributions for 1 or 2 sets of lambda windows."
    )
    parser.add_argument("base_dir1", type=str,
                        help="First base directory with lambda_* subdirectories.")
    parser.add_argument("base_dir2", nargs="?", default=None,
                        help="(Optional) second base directory with lambda_* subdirectories.")
    parser.add_argument("--energy_term", type=int, default=10,
                        help="Energy term index to extract (default: 10).")
    parser.add_argument("--units", type=str, default="kJ", choices=["kJ", "kcal"],
                        help="Energy units to plot (default: kJ).")
    args = parser.parse_args()

    plot_energy_density_vertical(args.base_dir1, args.base_dir2, energy_term=args.energy_term, units=args.units)

if __name__ == "__main__":
    main()