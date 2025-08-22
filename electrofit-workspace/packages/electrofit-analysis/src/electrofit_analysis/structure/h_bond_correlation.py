#!/usr/bin/env python3

import json
import logging
import sys
from typing import Dict, Tuple

import numpy as np

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def parse_xpm(file_path: str) -> Tuple[np.ndarray, Dict[str, str]]:
    """
    Minimal parse_xpm function:
    Converts an XPM file (from 'gmx hbond -hbm') into a binary (N_bonds x N_frames) matrix.
      1: bond present
      0: bond absent
    Returns (data_matrix, metadata).
    """
    try:
        with open(file_path, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        logging.error(f"XPM file not found: {file_path}")
        return np.array([[]]), {}

    metadata = {}
    color_map = {}
    data_lines = []
    header_found = False
    data_started = False
    width = height = num_colors = chars_per_pixel = 0

    for line in lines:
        line = line.strip()
        if line.startswith("/*"):
            # optional parse for metadata if desired
            continue
        if line.startswith("static char"):
            continue
        if line.startswith('"') and not header_found:
            # header: " width height ncolors chars_per_pixel"
            header = line.strip('",')
            tokens = header.split()
            if len(tokens) >= 4:
                width = int(tokens[0])
                height = int(tokens[1])
                num_colors = int(tokens[2])
                chars_per_pixel = int(tokens[3])
                header_found = True
            continue
        if header_found and not data_started:
            # read color definitions
            color_def = line.strip('",')
            if color_def == "":
                continue
            # parse e.g. "   c #FFFFFF"
            import re

            match = re.match(r"(.{%d})\s+c\s+(\S+)" % chars_per_pixel, color_def)
            if match:
                symbol = match.group(1)
                color = match.group(2)
                color_map[symbol] = color
            if len(color_map) == num_colors:
                data_started = True
            continue
        if data_started:
            if line.startswith('"'):
                data_line = line.strip('",')
                data_lines.append(data_line)

    if not data_lines:
        logging.warning(f"No data lines found in {file_path}")
        return np.array([[]]), {}

    data_matrix = np.zeros((height, width), dtype=int)
    for y, row_str in enumerate(data_lines):
        for x, char in enumerate(row_str):
            if char in color_map and color_map[char] == "#FF0000":
                data_matrix[y, x] = 1
            else:
                data_matrix[y, x] = 0

    return data_matrix, metadata


def compute_hbond_correlation_functions_optimized(
    data_matrix: np.ndarray, time_per_frame_ns: float, max_lag_frames: int = 1000
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Optimized version using vectorization and block detection.
    """
    n_bonds, n_frames = data_matrix.shape
    max_lag_frames = min(max_lag_frames, n_frames)

    # Precompute total bond appearances (for normalization)
    count_origins = data_matrix.sum()

    # 1. Compute S_intermittent using vectorized cross-correlation
    S_inter = np.zeros(max_lag_frames, dtype=float)
    for dt in range(max_lag_frames):
        if dt < n_frames:
            valid_frames = n_frames - dt
            S_inter[dt] = np.sum(
                data_matrix[:, :valid_frames] * data_matrix[:, dt : dt + valid_frames]
            )
    if count_origins > 0:
        S_inter /= count_origins

    # 2. Compute S_continuous using block detection
    S_cont = np.zeros(max_lag_frames, dtype=float)
    for bond_series in data_matrix:
        # Find start/end indices of continuous blocks
        diff = np.diff(bond_series, prepend=0, append=0)
        starts = np.where(diff == 1)[0]
        ends = np.where(diff == -1)[0]
        for start, end in zip(starts, ends):
            L = end - start  # Block length
            D = min(L, max_lag_frames)
            dt_values = np.arange(D)
            contributions = L - dt_values
            np.add.at(S_cont, dt_values, contributions)
    if count_origins > 0:
        S_cont /= count_origins

    times_ns = np.arange(max_lag_frames) * time_per_frame_ns
    return times_ns, S_inter, S_cont


def plot_correlation_functions(
    times_ns: np.ndarray, S_inter: np.ndarray, S_cont: np.ndarray, output_pdf: str
):
    """
    Plots the two correlation functions to a PDF.
    """
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 6))
    plt.plot(times_ns, S_inter, label="Intermittent $S_{HB}(t)$", color="blue")
    plt.plot(times_ns, S_cont, label="Continuous $S_{HB}^{(d)}(t)$", color="red")
    plt.xlabel("Time (ns)")
    plt.ylabel("Correlation")
    plt.title("H-Bond Correlation Functions")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_pdf, format="pdf")
    plt.close()
    logging.info(f"Saved correlation plot to {output_pdf}")


def main():
    """
    Minimal usage:
    python correlation_only.py path_to_matrix.xpm 0.01 1000
      where 0.01 is time_per_frame_ns, 1000 is max_lag_frames
    """
    import argparse

    parser = argparse.ArgumentParser(
        description="Compute and plot H-bond correlation functions from an XPM (presence/absence) file."
    )
    parser.add_argument(
        "xpm_file", type=str, help="Path to the XPM file (from gmx hbond -hbm)."
    )
    parser.add_argument(
        "--time_per_frame_ns",
        type=float,
        default=0.01,
        help="Time step in ns for each frame. Default=0.01 ns (10 ps).",
    )
    parser.add_argument(
        "--max_lag_frames",
        type=int,
        default=1000,
        help="Max number of frames for correlation time-lag. Default=1000.",
    )
    parser.add_argument(
        "--output_prefix",
        type=str,
        default="corr",
        help="Prefix for output PDF and JSON. Default='corr'",
    )
    args = parser.parse_args()

    xpm_file = args.xpm_file
    time_per_frame = args.time_per_frame_ns
    max_lag = args.max_lag_frames
    prefix = args.output_prefix

    logging.info(f"Loading XPM from {xpm_file}")
    data_matrix, _ = parse_xpm(xpm_file)
    if data_matrix.size == 0:
        logging.error("data_matrix is empty. Exiting.")
        sys.exit(1)

    logging.info(f"data_matrix shape = {data_matrix.shape}")

    logging.info("Computing correlation functions...")
    times_ns, S_inter, S_cont = compute_hbond_correlation_functions_optimized(
        data_matrix, time_per_frame_ns=time_per_frame, max_lag_frames=max_lag
    )

    # 1) Plot the correlation
    out_pdf = f"{prefix}_hb_corr_functions.pdf"
    plot_correlation_functions(times_ns, S_inter, S_cont, out_pdf)

    # 2) Dump numeric data to JSON
    corr_json = f"{prefix}_hb_corr_functions.json"
    corr_dict = {
        "times_ns": times_ns.tolist(),
        "S_intermittent": S_inter.tolist(),
        "S_continuous": S_cont.tolist(),
    }
    with open(corr_json, "w") as f:
        json.dump(corr_dict, f, indent=2)
    logging.info(f"Correlation arrays saved to {corr_json}")


if __name__ == "__main__":
    main()
