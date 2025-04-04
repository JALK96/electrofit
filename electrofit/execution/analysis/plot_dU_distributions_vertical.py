#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script to parse GROMACS dhdl.xvg for ΔU distributions, printing debug info
and normalizing each distribution so its maximum peak = 1.

Changes:
1) We now return (current_lambda, next_lambda, deltaU_values, col_index) from parse_single_subdir
   so gather_dU_per_subdir can print those details.
2) In plot_single_dir_dU, we do a manual KDE + scaling so that each distribution's peak = 1.
3) We add print statements showing which columns/first 10 ΔU values, etc.

Usage Example:
    python plot_dU_distributions_vertical.py /path/to/lambda_dir

Requirements:
    - numpy, matplotlib, seaborn, scipy
"""

import os
import re
import sys
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde

sns.set_context("talk")

def parse_dhdl_xvg_current_and_targets(dhdl_file, debug=False):
    """
    Parse a GROMACS dhdl.xvg to find:
      - current_lambda from 'fep-lambda = X' or the 'subtitle' line
      - a map of target_lambda -> column_index (s_index)
        specifically matching lines like:
          @ s2 legend "\\xD\\f{}H \\xl\\f{} to 0.0513"
    Returns:
      current_lambda (float or None)
      targets_map (dict: {float_value: column_index_in_data_arr})
      time_arr, data_arr
      legends_map (dict: sX -> legend_str) for debugging
    """
    current_lambda = None
    targets_map = {}
    legends_map = {}

    time_list = []
    data_list = []

    # Regex to catch e.g.: @ s2 legend "\xD\f{}H \xl\f{} to 0.0513"
    dh_regex = re.compile(
        r'@ s(\d+)\s+legend\s+"\\xD\\f\{\}H\s+\\xl\\f\{\}\s+to\s+([^"]+)"'
    )

    # fallback for e.g. @ subtitle "... fep-lambda = 0.0513"
    subtitle_lambda_regex = re.compile(
        r'@ subtitle ".*fep-lambda\s*=\s*([\d\.]+)"'
    )
    # also lines like: @ s1 legend "dH/d\xl\f{} fep-lambda = 0.0513"
    fep_lambda_line_regex = re.compile(
        r'fep-lambda\s*=\s*([\d\.]+)'
    )

    with open(dhdl_file, 'r', encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith('#'):
                continue

            if line.startswith('@'):
                # 1) Check subtitle with "fep-lambda="
                match_sub = subtitle_lambda_regex.search(line)
                if match_sub:
                    maybe_val = match_sub.group(1)
                    try:
                        current_lambda = float(maybe_val)
                    except ValueError:
                        pass

                # 2) Check "fep-lambda =" in the line
                match_fep = fep_lambda_line_regex.search(line)
                if match_fep:
                    maybe_val = match_fep.group(1)
                    try:
                        current_lambda = float(maybe_val)
                    except ValueError:
                        pass

                # 3) Our special pattern for lines: @ s(\d+) legend "\xD\f{}H \xl\f{} to <float>"
                dh_match = dh_regex.search(line)
                if dh_match:
                    s_index_str = dh_match.group(1)
                    val_str = dh_match.group(2).strip()
                    try:
                        s_index = int(s_index_str)
                        t_val = float(val_str)
                        # time col not counted => col = s_index
                        targets_map[t_val] = s_index
                    except ValueError:
                        pass

                # Also keep all legend lines in legends_map for debug
                leg_match = re.search(r'@ s(\d+)\s+legend\s+"([^"]+)"', line)
                if leg_match:
                    s_i = int(leg_match.group(1))
                    legend_txt = leg_match.group(2)
                    legends_map[s_i] = legend_txt

                continue

            # Otherwise numeric data line
            if line.startswith('@'):
                continue

            parts = line.split()
            floats = [float(x) for x in parts]
            time_list.append(floats[0])
            data_list.append(floats[1:])

    time_arr = np.array(time_list)
    data_arr = np.array(data_list)

    if debug:
        print(f"--- DEBUG parse_dhdl_xvg_current_and_targets({dhdl_file}) ---")
        print(f"Found current_lambda = {current_lambda}")
        print(f"Found targets_map = {targets_map}")
        print("Legends map:", legends_map)
        print(f"Data shape = {data_arr.shape}")

    return current_lambda, targets_map, time_arr, data_arr, legends_map


def find_next_lambda_and_column(current_lambda, targets_map, data_shape, offset_fix=0, debug=False):
    """
    Among the keys in targets_map (i.e. potential next lambdas),
    pick the smallest one that is > current_lambda.
    Then convert sX -> actual data column index = s_index + offset_fix.

    Returns: (next_lambda, col_index) or (None, None) if no bigger target or out of range
    """
    # gather all bigger than current_lambda
    bigger = [x for x in targets_map.keys() if x > current_lambda]
    if not bigger:
        return None, None

    next_lambda = min(bigger)
    s_index = targets_map[next_lambda]
    col_index = s_index + offset_fix

    # check range
    if col_index < 0 or col_index >= data_shape[1]:
        if debug:
            print(f"[DEBUG] col_index {col_index} out of range for data shape {data_shape}")
        return None, None

    return next_lambda, col_index


def parse_single_subdir(dhdl_file, offset_fix=0, debug=False):
    """
    1) parse the file to get current_lambda, target_map
    2) find next_lambda & col_index
    3) return (current_lambda, next_lambda, deltaU_values, col_index)
       or (None, None, None, None) if fails
    """
    (current_lambda,
     targets_map,
     time_arr,
     data_arr,
     legends_map) = parse_dhdl_xvg_current_and_targets(dhdl_file, debug=debug)

    if current_lambda is None or data_arr.size == 0:
        return None, None, None, None

    next_lambda, col_index = find_next_lambda_and_column(
        current_lambda, targets_map, data_arr.shape,
        offset_fix=offset_fix, debug=debug
    )
    if next_lambda is None or col_index is None:
        return current_lambda, None, None, None

    deltaU_values = data_arr[:, col_index]
    return current_lambda, next_lambda, deltaU_values, col_index


def gather_dU_per_subdir(base_dir, offset_fix=0, debug=False):
    """
    For each subdirectory 'lambda_*', parse its dhdl.xvg, returning:
      (cLam, nLam, deltaU, col_index, subdir_index)
    sorted by cLam.

    subdir_index is the integer Y from the folder name "lambda_Y"
    so we can print "state Y".
    """
    results = []
    subdirs = []
    for name in os.listdir(base_dir):
        if name.startswith("lambda_"):
            # parse the index from the directory name
            # e.g. "lambda_2" => subdir_index=2
            # if fails, store None
            try:
                subdir_index = float(name.split("_")[-1])
            except ValueError:
                subdir_index = None

            fpath = os.path.join(base_dir, name, "dhdl.xvg")
            if os.path.isfile(fpath):
                subdirs.append((subdir_index, fpath))

    # sort subdirs by subdir_index if possible
    subdirs.sort(key=lambda x: (x[0] if x[0] is not None else 1e9))

    out_list = []
    for (sdir_idx, dhdl_file) in subdirs:
        cLam, nLam, deltaU, col_idx = parse_single_subdir(dhdl_file, offset_fix=offset_fix, debug=debug)
        if cLam is None:
            print(f"[WARNING] Could not parse current lambda in {dhdl_file}")
            continue
        if nLam is None or deltaU is None:
            print(f"[INFO] No next-lambda found for subdir with current λ={cLam:.4f} "
                  f"in {dhdl_file}. Skipping.")
            continue

        out_list.append((cLam, nLam, deltaU, col_idx, sdir_idx))

    # sort final results by cLam
    out_list.sort(key=lambda x: x[0])
    return out_list


def plot_single_dir_dU(base_dir, invert_y=False, ax=None, offset_fix=0, debug=False, trim=0.1):
    """
    Gather (cLam, nLam, deltaU_values, col_index, subdir_idx), 
    then plot each distribution on 'ax' with color from 0->1 or 1->0.

    We'll do a manual KDE so we can enforce "peak = 1" normalization.
    Optionally, use 'trim' (in percent) to clip extreme outliers.
    """
    dU_list = gather_dU_per_subdir(base_dir, offset_fix=offset_fix, debug=debug)
    if len(dU_list) == 0:
        print(f"[WARNING] No ΔU transitions found in {base_dir}.")
        return

    n_trans = len(dU_list)
    cmap = plt.get_cmap("viridis")

    for i, (cLam, nLam, deltaU, col_idx, sdir_idx) in enumerate(dU_list):
        print(f"Reading transition λ={cLam:.4f} to λ={nLam:.4f} "
              f"from state {sdir_idx} (folder 'lambda_{sdir_idx}'), "
              f"extracting data from column {col_idx}.\n"
              f"First 10 ΔU values: {deltaU[:10]}")

        color_val = i / (n_trans - 1) if n_trans > 1 else 0.0
        color = cmap(color_val)

        label_str = None
        if n_trans <= 10:  # show legend if few lines
            label_str = f"λ={cLam:.4f}->{nLam:.4f}"

        if len(deltaU) < 2:
            continue  # skip trivial cases

        kde = gaussian_kde(deltaU)

        # Determine x-range: use trimmed percentiles if 'trim' is provided,
        # otherwise use min and max of deltaU.
        if trim is not None and trim > 0 and trim < 50:
            left = np.percentile(deltaU, trim)
            right = np.percentile(deltaU, 100 - trim)
            if left == right:
                left -= 0.5
                right += 0.5
        else:
            left = np.min(deltaU)
            right = np.max(deltaU)

        pad = 0.05 * (right - left)
        xvals = np.linspace(left - pad, right + pad, 200)
        yvals = kde(xvals)

        # Normalize so that the peak value is 1.
        peak = np.max(yvals)
        if peak > 0:
            yvals /= peak

        ax.plot(xvals, yvals, color=color, lw=1.5, label=label_str)

    if invert_y:
        ax.invert_yaxis()

    ax.set_ylabel(r"$P(\Delta U)$")
    if n_trans <= 10:
        ax.legend()


def plot_dU_distributions_vertical(base_dir1, base_dir2=None, offset_fix=0, debug=False):
    """
    If only base_dir1 is given -> single subplot
    If base_dir2 is also given -> 2-subplot vertical
    """
    if base_dir2 is None:
        fig, ax_top = plt.subplots(1, 1, figsize=(10, 5))
        plot_single_dir_dU(base_dir1, invert_y=False, ax=ax_top, offset_fix=offset_fix, debug=debug)
        ax_top.set_xlabel(r"ΔU$_{ij}$ (kJ/mol)")
        ax_top.set_title("A → B")

        # colorbar from 0->1
        norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=ax_top, orientation="vertical", fraction=0.05)
        cbar.set_label(r"$\lambda_j$")

        #plt.tight_layout()
        plt.savefig("plot_dU_single.pdf")
        #plt.show()
    else:
        fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        # top
        plot_single_dir_dU(base_dir1, invert_y=False, ax=ax_top, offset_fix=offset_fix, debug=debug)
        ax_top.set_ylabel(r"Norm: $P_{fwd}(\Delta U$")
        ax_top.set_title("A → B")

        # bottom
        plot_single_dir_dU(base_dir2, invert_y=True, ax=ax_bottom, offset_fix=offset_fix, debug=debug)
        ax_bottom.set_ylabel(r"Norm: $P_{bwd}(\Delta U$")
        ax_bottom.set_xlabel(r"ΔU$_{ij}$ (kJ/mol)")
        ax_bottom.set_title("B → A")

        # single colorbar
        norm = matplotlib.colors.Normalize(vmin=0, vmax=1)
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        sm.set_array([])
        cbar = fig.colorbar(sm, ax=[ax_top, ax_bottom], orientation="vertical", fraction=0.05)
        cbar.set_label(r"$\lambda_j$")

        #plt.tight_layout()
        plt.savefig("plot_dU_double.pdf")
        #plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Plot ΔU(λ_i->λ_{i+1}) distributions from dhdl.xvg, forcing each peak=1, printing debug info."
    )
    parser.add_argument("base_dir1", help="First directory with subfolders 'lambda_*'.")
    parser.add_argument("base_dir2", nargs="?", default=None,
                        help="(Optional) second directory for two-subplot layout.")
    parser.add_argument("--offset", type=int, default=0,
                        help="Integer offset to apply to sX->data column (if needed).")
    parser.add_argument("--debug", action="store_true", help="Enable debug prints.")
    args = parser.parse_args()

    base1 = os.path.abspath(args.base_dir1)
    base2 = os.path.abspath(args.base_dir2) if args.base_dir2 else None

    plot_dU_distributions_vertical(base1, base2, offset_fix=args.offset, debug=args.debug)


if __name__ == "__main__":
    main()