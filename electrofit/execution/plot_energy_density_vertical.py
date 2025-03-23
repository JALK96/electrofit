#!/usr/bin/env python3
"""
Script: plot_energy_density_vertical.py

This script extracts energy data from GROMACS binary energy files (ener.edr) located
in lambda window subdirectories for two different simulation setups. For each base directory,
the script does the following:
  1. Iterates over subdirectories whose names start with "lambda_" (e.g., lambda_0, lambda_1, …)
     and sorts them in order.
  2. Runs the 'gmx energy' command to extract a specific energy term (e.g., potential energy)
     from the binary energy file, writing the output to an .xvg file.
  3. Parses the .xvg file to obtain the numerical energy data.
  4. Computes a kernel density estimate (KDE) of the energy distribution using Seaborn.

The output figure contains two subplots arranged vertically (top and bottom) that share the x‑axis:
  - The top subplot shows the energy density distributions from the first base directory.
  - The bottom subplot shows the energy density distributions from the second base directory,
    but with the y‑axis inverted.

Usage:
    python plot_energy_density_vertical.py <base_dir1> <base_dir2> [--energy_term <int>]

Arguments:
    base_dir1      Path to the first base directory containing "lambda_*" subdirectories.
    base_dir2      Path to the second base directory containing "lambda_*" subdirectories.
    --energy_term  (Optional) Index of the energy term to extract (default: 10).

Requirements:
    - GROMACS must be installed and accessible in the PATH (to run gmx energy).
    - Python packages: numpy, matplotlib, seaborn.
    - Each base directory should contain subdirectories named "lambda_0", "lambda_1", …,
      and each such directory must contain an "ener.edr" file.
"""

import os
import sys
import subprocess
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def extract_energy_data(ener_edr, output_xvg, energy_term=10):
    """
    Extracts a specific energy term from a binary ener.edr file using gmx energy,
    and writes the result to an .xvg file.

    Parameters
    ----------
    ener_edr : str
        Path to the input binary energy file (ener.edr).
    output_xvg : str
        Path where the output .xvg file will be saved.
    energy_term : int, optional
        The index of the energy term to extract (default: 10).
    """
    # We pipe the energy_term and a terminating "0" into gmx energy.
    cmd = f'echo "{energy_term}\n0" | gmx energy -f {ener_edr} -o {output_xvg} -quiet'
    subprocess.run(cmd, shell=True, check=True)

def parse_xvg(xvg_file):
    """
    Parses a GROMACS .xvg file and returns its numeric data as a NumPy array.
    Lines starting with '#' or '@' are skipped.

    Parameters
    ----------
    xvg_file : str
        Path to the .xvg file.

    Returns
    -------
    numpy.ndarray
        2D array containing the numeric data.
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

def plot_energy_density_vertical(base_dir1, base_dir2, energy_term=10):
    """
    Plots energy density distributions for two sets of lambda windows.
    The figure has two vertically arranged subplots that share the x-axis.
      - Top subplot: Data from base_dir1.
      - Bottom subplot: Data from base_dir2 with the y-axis inverted.

    Parameters
    ----------
    base_dir1 : str
        Directory containing lambda window subdirectories for the first simulation setup.
    base_dir2 : str
        Directory containing lambda window subdirectories for the second simulation setup.
    energy_term : int, optional
        The energy term index to extract (default: 10).
    """
    # Get sorted list of lambda directories for each base directory.
    lambda_dirs1 = sorted(
        [os.path.join(base_dir1, d) for d in os.listdir(base_dir1) if d.startswith("lambda_")],
        key=lambda d: float(d.split("_")[-1])
    )
    lambda_dirs2 = sorted(
        [os.path.join(base_dir2, d) for d in os.listdir(base_dir2) if d.startswith("lambda_")],
        key=lambda d: float(d.split("_")[-1])
    )
    
    n1 = len(lambda_dirs1)
    n2 = len(lambda_dirs2)
    if n1 == 0 or n2 == 0:
        print("No lambda directories found in one or both base directories.")
        return

    # Create subplots arranged vertically and sharing the x-axis.
    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    cmap = plt.get_cmap("viridis")

    # Plot data for base_dir1 (top subplot).
    for idx, lam_dir in enumerate(lambda_dirs1):
        ener_edr = os.path.join(lam_dir, "ener.edr")
        output_xvg = os.path.join(lam_dir, "energy.xvg")
        extract_energy_data(ener_edr, output_xvg, energy_term=energy_term)
        data = parse_xvg(output_xvg)
        if data.size == 0:
            print(f"No data found in {output_xvg}.")
            continue
        # Assume energy is in the second column.
        energy = data[:, 1]
        color = cmap(idx / (n1 - 1)) if n1 > 1 else cmap(0)
        sns.kdeplot(energy, ax=ax_top, color=color, label=os.path.basename(lam_dir))
    ax_top.set_ylabel("Density")
    ax_top.set_title("Energy Density Distribution (Base Dir 1)")

    # Plot data for base_dir2 (bottom subplot) with inverted y-axis.
    for idx, lam_dir in enumerate(lambda_dirs2):
        ener_edr = os.path.join(lam_dir, "ener.edr")
        output_xvg = os.path.join(lam_dir, "energy.xvg")
        extract_energy_data(ener_edr, output_xvg, energy_term=energy_term)
        data = parse_xvg(output_xvg)
        if data.size == 0:
            print(f"No data found in {output_xvg}.")
            continue
        energy = data[:, 1]
        color = cmap(idx / (n2 - 1)) if n2 > 1 else cmap(0)
        sns.kdeplot(energy, ax=ax_bottom, color=color, label=os.path.basename(lam_dir))
    ax_bottom.set_ylabel("Density (Inverted)")
    ax_bottom.invert_yaxis()
    ax_bottom.set_xlabel("Energy (kJ/mol)")
    ax_bottom.set_title("Energy Density Distribution (Base Dir 2)")

    # Add legends to both subplots.
    ax_top.legend(title="Lambda Window", loc="upper right")
    ax_bottom.legend(title="Lambda Window", loc="upper right")
    plt.tight_layout()
    plt.savefig("plot.pdf")

def main():
    parser = argparse.ArgumentParser(
        description="Plot vertical energy density distributions for two sets of lambda windows."
    )
    parser.add_argument("base_dir1", type=str,
                        help="First base directory with lambda_* subdirectories.")
    parser.add_argument("base_dir2", type=str,
                        help="Second base directory with lambda_* subdirectories.")
    parser.add_argument("--energy_term", type=int, default=10,
                        help="Energy term index to extract (default: 10).")
    args = parser.parse_args()
    plot_energy_density_vertical(args.base_dir1, args.base_dir2, energy_term=args.energy_term)

if __name__ == "__main__":
    main()