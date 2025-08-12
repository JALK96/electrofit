#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
import os
import re

from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
import matplotlib.pyplot as plt
import numpy as np
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
    subtitle_lambda_regex = re.compile(r'@ subtitle ".*fep-lambda\s*=\s*([\d\.]+)"')
    # also lines like: @ s1 legend "dH/d\xl\f{} fep-lambda = 0.0513"
    fep_lambda_line_regex = re.compile(r"fep-lambda\s*=\s*([\d\.]+)")

    with open(dhdl_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            if line.startswith("@"):
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
            if line.startswith("@"):
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


def find_next_lambda_and_column(
    current_lambda, targets_map, data_shape, offset_fix=0, debug=False
):
    """
    Among the keys in targets_map (i.e. potential next lambdas),
    pick the smallest one that is > current_lambda.
    Then convert sX -> actual data column index = s_index + offset_fix.

    Returns: (next_lambda, col_index) or (None, None) if no bigger target or out of range
    """
    bigger = [x for x in targets_map.keys() if x > current_lambda]
    if not bigger:
        return None, None

    next_lambda = min(bigger)
    s_index = targets_map[next_lambda]
    col_index = s_index + offset_fix

    if col_index < 0 or col_index >= data_shape[1]:
        if debug:
            print(
                f"[DEBUG] col_index {col_index} out of range for data shape {data_shape}"
            )
        return None, None

    return next_lambda, col_index


def parse_single_subdir(dhdl_file, offset_fix=0, debug=False):
    """
    1) parse the file to get current_lambda, target_map
    2) find next_lambda & col_index
    3) return (current_lambda, next_lambda, deltaU_values, col_index)
       or (None, None, None, None) if fails
    """
    (current_lambda, targets_map, time_arr, data_arr, legends_map) = (
        parse_dhdl_xvg_current_and_targets(dhdl_file, debug=debug)
    )

    if current_lambda is None or data_arr.size == 0:
        return None, None, None, None

    next_lambda, col_index = find_next_lambda_and_column(
        current_lambda, targets_map, data_arr.shape, offset_fix=offset_fix, debug=debug
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
    """
    subdirs = []
    for name in os.listdir(base_dir):
        if name.startswith("lambda_"):
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
    for sdir_idx, dhdl_file in subdirs:
        cLam, nLam, deltaU, col_idx = parse_single_subdir(
            dhdl_file, offset_fix=offset_fix, debug=debug
        )
        if cLam is None:
            print(f"[WARNING] Could not parse current lambda in {dhdl_file}")
            continue
        if nLam is None or deltaU is None:
            print(
                f"[INFO] No next-lambda found for subdir with current λ={cLam:.4f} "
                f"in {dhdl_file}. Skipping."
            )
            continue

        out_list.append((cLam, nLam, deltaU, col_idx, sdir_idx))

    out_list.sort(key=lambda x: x[0])
    return out_list


def compute_density(deltaU, normalize=False):
    """
    Helper to compute xvals, yvals for a Gaussian KDE.
    If normalize is True, scales the peak so that max(yvals)=1.
    Otherwise, returns the raw density.
    """
    kde = gaussian_kde(deltaU)
    left = np.min(deltaU)
    right = np.max(deltaU)

    # ensure we don't have zero range
    if left == right:
        left -= 0.5
        right += 0.5

    pad = 0.05 * (right - left)
    xvals = np.linspace(left - pad, right + pad, 200)
    yvals = kde(xvals)
    if normalize:
        peak = np.max(yvals)
        if peak > 0:
            yvals /= peak
    return xvals, yvals





def plot_dU_distributions_grid(
    base_dir1, base_dir2=None, offset_fix=0, debug=False, no_normalize=False
):
    """
    Gathers transitions from base_dir1 (forward) and optionally base_dir2 (backward).
    Creates a grid of subplots with 4 columns, enough rows to fit all transitions.
    Each subplot corresponds to one transition index, possibly with two distributions (fwd + bwd).

    no_normalize=True => do NOT scale the peak to 1.
    """
    # Gather forward transitions
    fwd_list = gather_dU_per_subdir(base_dir1, offset_fix=offset_fix, debug=debug)
    # Gather backward transitions if provided
    if base_dir2 is not None:
        bwd_list = gather_dU_per_subdir(base_dir2, offset_fix=offset_fix, debug=debug)
    else:
        bwd_list = []

    n_fwd = len(fwd_list)
    n_bwd = len(bwd_list)
    n_plots = max(n_fwd, n_bwd)

    if n_plots == 0:
        print("[WARNING] No transitions found in the provided directories.")
        return

    # Determine grid size: 4 columns and as many rows as needed
    ncols = 4
    nrows = math.ceil(n_plots / ncols)

    # Create the figure and axes
    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(4 * ncols, 3 * nrows),
        sharex=False,
        sharey=False,
    )

    # Ensure axs is a 2D array for consistent indexing
    axs = np.atleast_2d(axs)

    for i in range(n_plots):
        row = i // ncols
        col = i % ncols
        ax = axs[row, col]

        # Forward distribution if available
        if i < n_fwd:
            cLam_f, nLam_f, deltaU_f, _, sdir_idx_f = fwd_list[i]
            if debug:
                print(
                    f"[DEBUG] FWD i={i} λ={cLam_f:.4f} -> λ={nLam_f:.4f},"
                    f" first 10 ΔU: {deltaU_f[:10]}"
                )
            if len(deltaU_f) > 1:
                xvals_f, yvals_f = compute_density(deltaU_f, normalize=not no_normalize)
                ax.plot(xvals_f, yvals_f, c="k")

        # Backward distribution if available
        if i < n_bwd:
            cLam_b, nLam_b, deltaU_b, _, sdir_idx_b = bwd_list[i]
            if debug:
                print(
                    f"[DEBUG] BWD i={i} λ={cLam_b:.4f} -> λ={nLam_b:.4f},"
                    f" first 10 ΔU: {deltaU_b[:10]}"
                )
            if len(deltaU_b) > 1:
                xvals_b, yvals_b = compute_density(deltaU_b, normalize=not no_normalize)
                ax.plot(xvals_b, yvals_b, "--", c="red")

        # Set title using backward state values if available, otherwise a generic title.
        if i < n_bwd:
            ax.set_title(rf"$\lambda={cLam_b:.3f} \to {nLam_b:.3f}$")
        else:
            ax.set_title(rf"$\lambda={cLam_f:.3f} \to {nLam_f:.3f}$")

    # Hide extra axes if nrows*ncols > n_plots
    total_subplots = nrows * ncols
    for j in range(n_plots, total_subplots):
        row = j // ncols
        col = j % ncols
        axs[row, col].set_visible(False)

    # Remove individual axes x- and y-labels if present
    for ax in axs.flatten():
        ax.set_xlabel("")
        ax.set_ylabel("")

    # Add a common x-axis label at the bottom
    fig.text(0.5, 0.08, r"$\Delta U$ (kJ/mol)", ha="center", va="center", fontsize=30)

    # Create a composite y-axis label using offsetbox
    if no_normalize:
        # Create separate text areas for the two parts of the label
        ybox1 = TextArea(
            r"$P_{fwd}(\Delta U)$",
            textprops=dict(color="k", size=30, rotation=90, ha="center", va="bottom"),
        )
        ybox2 = TextArea(
            ", ",
            textprops=dict(color="k", size=30, rotation=90, ha="center", va="bottom"),
        )
        ybox3 = TextArea(
            r"$P_{bwd}(\Delta U)$",
            textprops=dict(color="red", size=30, rotation=90, ha="center", va="bottom"),
        )
        ybox = VPacker(children=[ybox3, ybox2, ybox1], align="center", pad=0, sep=2)
    else:
        ybox = TextArea(
            "Scaled Density",
            textprops=dict(color="k", size=30, rotation=90, ha="center", va="bottom"),
        )

    # Anchor the composite y-label to the figure.
    anchored_ybox = AnchoredOffsetbox(
        loc="center left",
        child=ybox,
        pad=0.0,
        frameon=False,
        bbox_to_anchor=(0.08, 0.5),
        bbox_transform=fig.transFigure,
        borderpad=0.0,
    )
    fig.add_artist(anchored_ybox)

    # Adjust layout to account for the common labels
    plt.tight_layout(rect=[0.1, 0.1, 1, 0.95])
    plt.savefig("plot_dU_grid.pdf")
    print("Saved Figure to: plot_dU_grid.pdf")
    # plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Plot ΔU(λ_i->λ_{i+1}) distributions in a 4-column grid, each subplot is one transition."
    )
    parser.add_argument(
        "base_dir1", help="First directory with subfolders 'lambda_*' (forward)."
    )
    parser.add_argument(
        "base_dir2",
        nargs="?",
        default=None,
        help="(Optional) second directory with subfolders 'lambda_*' (backward).",
    )
    parser.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Integer offset to apply to sX->data column (if needed).",
    )
    parser.add_argument(
        "--no-normalize",
        action="store_true",
        help="Disable peak=1 normalization and plot raw density instead.",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug prints.")
    args = parser.parse_args()

    base1 = os.path.abspath(args.base_dir1)
    base2 = os.path.abspath(args.base_dir2) if args.base_dir2 else None

    plot_dU_distributions_grid(
        base1,
        base2,
        offset_fix=args.offset,
        debug=args.debug,
        no_normalize=args.no_normalize,
    )


if __name__ == "__main__":
    main()
