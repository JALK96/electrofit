#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import math
import os
import re

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
from scipy.stats import gaussian_kde

sns.set_context("talk", font_scale=1.2)

#####################################################################
#  Parsing & Data Gathering
#####################################################################


def parse_dhdl_xvg_current_and_targets(dhdl_file, debug=False, skip_time=0.0):
    """
    Parse a GROMACS dhdl.xvg to find:
      - current_lambda from 'fep-lambda = X' or the 'subtitle' line
      - a map of target_lambda -> column_index (s_index)
    Also reads time_arr and data_arr from the file.

    If skip_time > 0, discard all rows where time < skip_time.

    Returns:
      current_lambda (float or None)
      targets_map (dict: {float_value: column_index_in_data_arr})
      time_arr, data_arr (filtered if skip_time>0)
      legends_map (dict: sX -> legend_str) for debugging
    """
    current_lambda = None
    targets_map = {}
    legends_map = {}

    time_list = []
    data_list = []

    # Regex patterns
    dh_regex = re.compile(
        r'@ s(\d+)\s+legend\s+"\\xD\\f\{\}H\s+\\xl\\f\{\}\s+to\s+([^"]+)"'
    )
    subtitle_lambda_regex = re.compile(r'@ subtitle ".*fep-lambda\s*=\s*([\d\.]+)"')
    fep_lambda_line_regex = re.compile(r"fep-lambda\s*=\s*([\d\.]+)")

    with open(dhdl_file, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            if line.startswith("@"):
                # 1) "fep-lambda" in subtitle
                match_sub = subtitle_lambda_regex.search(line)
                if match_sub:
                    maybe_val = match_sub.group(1)
                    try:
                        current_lambda = float(maybe_val)
                    except ValueError:
                        pass

                # 2) "fep-lambda =" in the line
                match_fep = fep_lambda_line_regex.search(line)
                if match_fep:
                    maybe_val = match_fep.group(1)
                    try:
                        current_lambda = float(maybe_val)
                    except ValueError:
                        pass

                # 3) Pattern @ s(\d+) legend "\xD\f{}H \xl\f{} to <float>"
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

                # Also keep legend lines in legends_map for debug
                leg_match = re.search(r'@ s(\d+)\s+legend\s+"([^"]+)"', line)
                if leg_match:
                    s_i = int(leg_match.group(1))
                    legend_txt = leg_match.group(2)
                    legends_map[s_i] = legend_txt

                continue

            # Numeric data line
            if line.startswith("@"):
                continue

            parts = line.split()
            floats = [float(x) for x in parts]
            time_list.append(floats[0])
            data_list.append(floats[1:])

    time_arr = np.array(time_list)
    data_arr = np.array(data_list)

    # Skip data prior to skip_time
    if skip_time > 0:
        mask = time_arr >= skip_time
        time_arr = time_arr[mask]
        data_arr = data_arr[mask, :]

    if debug:
        print(f"--- DEBUG parse_dhdl_xvg_current_and_targets({dhdl_file}) ---")
        print(f"Found current_lambda = {current_lambda}")
        print(f"Found targets_map = {targets_map}")
        print("Legends map:", legends_map)
        print(f"Data shape after skip_time={skip_time} filtering: {data_arr.shape}")

    return current_lambda, targets_map, time_arr, data_arr, legends_map


def find_next_lambda_and_column(
    current_lambda, targets_map, data_shape, offset_fix=0, debug=False
):
    """
    Among the keys in targets_map, pick the smallest one that is > current_lambda.
    Then convert sX -> actual data column index = s_index + offset_fix.
    Returns: (next_lambda, col_index) or (None, None).
    """
    bigger = [x for x in targets_map.keys() if x > current_lambda]
    if not bigger:
        return None, None

    next_lambda = min(bigger)
    s_index = targets_map[next_lambda]
    col_index = s_index + offset_fix

    if col_index < 0 or col_index >= data_shape[1]:
        if debug:
            print(f"[DEBUG] col_index {col_index} out of range for shape {data_shape}")
        return None, None

    return next_lambda, col_index


def parse_single_subdir(dhdl_file, offset_fix=0, debug=False, skip_time=0.0):
    """
    Parse the file, find the relevant column for ΔU, return the data.
    Returns (current_lambda, next_lambda, deltaU_values, col_index).
    """
    (current_lambda, targets_map, time_arr, data_arr, legends_map) = (
        parse_dhdl_xvg_current_and_targets(dhdl_file, debug=debug, skip_time=skip_time)
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


def gather_dU_per_subdir(base_dir, offset_fix=0, debug=False, skip_time=0.0):
    """
    For each 'lambda_*' folder in base_dir, parse its dhdl.xvg.
    Return sorted list of (cLam, nLam, deltaU, col_index, subdir_index).
    Skips data prior to skip_time.
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

    # Sort
    subdirs.sort(key=lambda x: (x[0] if x[0] is not None else 1e9))

    out_list = []
    for sdir_idx, dhdl_file in subdirs:
        cLam, nLam, deltaU, col_idx = parse_single_subdir(
            dhdl_file, offset_fix=offset_fix, debug=debug, skip_time=skip_time
        )
        if cLam is None:
            print(f"[WARNING] Could not parse current lambda in {dhdl_file}")
            continue
        if nLam is None or deltaU is None:
            print(
                f"[INFO] No next-lambda found (or no data left) for subdir with λ={cLam:.4f} in {dhdl_file}. Skipping."
            )
            continue

        out_list.append((cLam, nLam, deltaU, col_idx, sdir_idx))

    out_list.sort(key=lambda x: x[0])
    return out_list


#####################################################################
#  Plotting Consecutive Transitions with skip_time + trim + (de)normalization
#####################################################################


def compute_density(deltaU, normalize=False, trim_percentile=None):
    """
    Compute kernel density for the data 'deltaU'.
    - If 'normalize' is True, scale so max(y)=1.
    - If trim_percentile is not None (e.g. 2), we clamp x-range
      to [2nd percentile, 98th percentile], ignoring outliers for display.
    """
    if len(deltaU) < 2:
        return np.array([]), np.array([])

    kde = gaussian_kde(deltaU)

    # 1) Determine left/right
    if trim_percentile is not None and 0 < trim_percentile < 50:
        lower = np.percentile(deltaU, trim_percentile)
        upper = np.percentile(deltaU, 100 - trim_percentile)
        if lower == upper:
            lower -= 0.5
            upper += 0.5
    else:
        lower = np.min(deltaU)
        upper = np.max(deltaU)
        if lower == upper:
            lower -= 0.5
            upper += 0.5

    span = upper - lower
    pad = 0.05 * span
    left = lower - pad
    right = upper + pad

    xvals = np.linspace(left, right, 200)
    yvals = kde(xvals)

    # 2) Normalize if requested
    if normalize:
        peak = np.max(yvals)
        if peak > 0:
            yvals /= peak

    return xvals, yvals


def plot_dU_distributions_consecutive_grid(
    base_dir1,
    base_dir2=None,
    offset_fix=0,
    debug=False,
    normalize=False,
    trim=None,
    skip_time=0.0,
    units="kJ",
):
    """
    Plots each transition i (solid) & i+1 (dashed) for forward (black)
    and backward (red), skipping data before skip_time, optionally trimming
    outliers, and optionally normalizing (peak=1).

    The energy units can be set via the 'units' parameter ("kJ" [default] or "kcal").
    """
    # Gather data
    fwd_list = gather_dU_per_subdir(
        base_dir1, offset_fix=offset_fix, debug=debug, skip_time=skip_time
    )
    bwd_list = []
    if base_dir2 is not None:
        bwd_list = gather_dU_per_subdir(
            base_dir2, offset_fix=offset_fix, debug=debug, skip_time=skip_time
        )

    n_fwd = len(fwd_list)
    n_bwd = len(bwd_list)
    n_plots = max(n_fwd, n_bwd)

    if n_plots == 0:
        print(
            "[WARNING] No transitions found (or no data after skip_time) in the provided directories."
        )
        return

    # 4 columns
    ncols = 5
    nrows = math.ceil(n_plots / ncols)

    fig, axs = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(4 * ncols, 3 * nrows),
        sharex=False,
        sharey=False,
    )
    axs = np.atleast_2d(axs)

    for i in range(n_plots):
        row = i // ncols
        col = i % ncols
        ax = axs[row, col]

        # Forward i => black solid
        if i < n_fwd:
            cLam_f_i, nLam_f_i, deltaU_f_i, _, _ = fwd_list[i]
            if units.lower() == "kcal":
                deltaU_f_i = deltaU_f_i * 0.239
            x_f_i, y_f_i = compute_density(
                deltaU_f_i, normalize=normalize, trim_percentile=trim
            )
            if x_f_i.size > 0:
                ax.plot(x_f_i, y_f_i, "-", color="black", label=f"Fwd i={i}")

        # Forward i+1 => black dashed
        if (i + 1) < n_fwd:
            cLam_f_ip1, nLam_f_ip1, deltaU_f_ip1, _, _ = fwd_list[i + 1]
            if units.lower() == "kcal":
                deltaU_f_ip1 = deltaU_f_ip1 * 0.239
            x_f_ip1, y_f_ip1 = compute_density(
                deltaU_f_ip1, normalize=normalize, trim_percentile=trim
            )
            if x_f_ip1.size > 0:
                ax.plot(x_f_ip1, y_f_ip1, "--", color="black", label=f"Fwd i+1={i + 1}")

        # Backward i => red solid
        if i < n_bwd:
            cLam_b_i, nLam_b_i, deltaU_b_i, _, _ = bwd_list[i]
            if units.lower() == "kcal":
                deltaU_b_i = deltaU_b_i * 0.239
            x_b_i, y_b_i = compute_density(
                deltaU_b_i, normalize=normalize, trim_percentile=trim
            )
            if x_b_i.size > 0:
                ax.plot(x_b_i, y_b_i, "-", color="crimson", label=f"Bwd i={i}")

        # Backward i+1 => red dashed
        if (i + 1) < n_bwd:
            cLam_b_ip1, nLam_b_ip1, deltaU_b_ip1, _, _ = bwd_list[i + 1]
            if units.lower() == "kcal":
                deltaU_b_ip1 = deltaU_b_ip1 * 0.239
            x_b_ip1, y_b_ip1 = compute_density(
                deltaU_b_ip1, normalize=normalize, trim_percentile=trim
            )
            if x_b_ip1.size > 0:
                ax.plot(
                    x_b_ip1, y_b_ip1, "--", color="crimson", label=f"Bwd i+1={i + 1}"
                )

        ax.set_title(rf"{cLam_f_i:.4f}→{nLam_f_i:.4f}")

    # Hide extra subplots
    total_subplots = nrows * ncols
    for j in range(n_plots, total_subplots):
        r = j // ncols
        c = j % ncols
        axs[r, c].set_visible(False)

    # Remove per-axis x,y labels
    for ax in axs.flatten():
        ax.set_xlabel("")
        ax.set_ylabel("")

    # Common X-axis label with units
    fig.text(
        0.5, 0.08, rf"$\Delta U$ ({units}/mol)", ha="center", va="center", fontsize=30
    )

    # Decide Y-axis label depending on whether we have both dirs and whether we are normalizing
    if base_dir2 is None:
        if normalize:
            ybox = TextArea(
                "Scaled Density (Fwd)",
                textprops=dict(
                    color="k", size=30, rotation=90, ha="center", va="bottom"
                ),
            )
        else:
            ybox = TextArea(
                r"$P_{fwd}(\Delta U)$",
                textprops=dict(
                    color="k", size=30, rotation=90, ha="center", va="bottom"
                ),
            )
    else:
        if normalize:
            ybox = TextArea(
                "Scaled Density",
                textprops=dict(
                    color="k", size=30, rotation=90, ha="center", va="bottom"
                ),
            )
        else:
            ybox1 = TextArea(
                r"$P_{fwd}(\Delta U)$",
                textprops=dict(
                    color="k", size=30, rotation=90, ha="center", va="bottom"
                ),
            )
            ybox2 = TextArea(
                ", ",
                textprops=dict(
                    color="k", size=30, rotation=90, ha="center", va="bottom"
                ),
            )
            ybox3 = TextArea(
                r"$P_{bwd}(\Delta U)$",
                textprops=dict(
                    color="crimson", size=30, rotation=90, ha="center", va="bottom"
                ),
            )
            ybox = VPacker(children=[ybox3, ybox2, ybox1], align="center", pad=0, sep=2)

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

    plt.tight_layout(rect=[0.1, 0.1, 1, 0.97])
    outname = "plot_dU_consecutive_trim_skip.pdf"
    plt.savefig(outname)
    print(f"Saved figure to {outname}")
    # plt.show()


#####################################################################
#  Main
#####################################################################


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Plot consecutive ΔU(λ_i→λ_{i+1}) distributions, skipping data before "
            "a specified time, optionally trimming outliers, optionally normalizing."
        )
    )
    parser.add_argument(
        "base_dir1", help="First directory with subfolders 'lambda_*' (forward)."
    )
    parser.add_argument(
        "base_dir2",
        nargs="?",
        default=None,
        help="(Optional) second directory (backward).",
    )
    parser.add_argument(
        "--offset", type=int, default=0, help="Offset to apply to sX→data column."
    )
    parser.add_argument(
        "--normalize",
        action="store_true",
        help="If set, scale each distribution so its peak=1.",
    )
    parser.add_argument(
        "--trim",
        type=float,
        default=None,
        help=(
            "If set, clamp x-range to e.g. [trim%%, 100-trim%%] to "
            "ignore extreme outliers (e.g. '--trim 2')."
        ),
    )
    parser.add_argument(
        "--skip-time",
        type=float,
        default=0.0,
        help="Ignore all data rows with time < skip_time (ps).",
    )
    parser.add_argument("--debug", action="store_true", help="Enable debug prints.")
    parser.add_argument(
        "--units",
        type=str,
        default="kJ",
        choices=["kJ", "kcal"],
        help="Energy units to plot (default: kJ).",
    )
    args = parser.parse_args()

    base1 = os.path.abspath(args.base_dir1)
    base2 = os.path.abspath(args.base_dir2) if args.base_dir2 else None

    plot_dU_distributions_consecutive_grid(
        base_dir1=base1,
        base_dir2=base2,
        offset_fix=args.offset,
        debug=args.debug,
        normalize=args.normalize,
        trim=args.trim,
        skip_time=args.skip_time,
        units=args.units,
    )


if __name__ == "__main__":
    main()
