#!/usr/bin/env python3
"""
plot_free_energy.py

Usage:
  ./plot_free_energy.py <directory>

This command-line tool reads the CSV file (default "free_energy_summary.csv")
from the specified directory and produces a grouped bar plot using the
plotting function. If the directory string contains "stateB", it uses invert_TI=True.

The resulting plot is saved as "breakdown.pdf" in the given directory.
"""

import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_context("talk", font_scale=0.9)


def parse_results_csv(file_path, invert_TI=True):
    """
    Parse a CSV file that contains a results table and extract the data.

    The CSV file is expected to have the following columns:
      state, name, MBAR, MBAR_Error, BAR, BAR_Error, TI, TI_Error

    The 'name' column represents the lambda pair (e.g. "0 -- 1") and will be
    reformatted with an en-dash (–). Rows with a 'state' value of "Stages" or "TOTAL" are ignored.

    Parameters:
      file_path (str): Path to the CSV file.
      invert_TI (bool): If True, multiply the TI values by -1 and invert the lambda pair order.

    Returns:
      pd.DataFrame: DataFrame with selected columns.
    """
    # Read the CSV file.
    df = pd.read_csv(file_path)

    # Filter out rows with 'state' indicating summary entries.
    df = df[~df["state"].str.lower().isin(["stages", "total"])]

    # Rename 'name' column to 'Lambda pair' and replace "--" with an en-dash.
    df.rename(columns={"name": "Lambda pair"}, inplace=True)
    df["Lambda pair"] = df["Lambda pair"].str.replace("--", "–").str.strip()

    if invert_TI:
        df["TI"] = -df["TI"]
        df["TI"] = df["TI"].iloc[::-1].reset_index(drop=True)
        df["TI_Error"] = df["TI_Error"].iloc[::-1].reset_index(drop=True)
        df["Lambda pair"] = df["Lambda pair"].iloc[::-1].reset_index(drop=True)
        df["Lambda pair"] = df["Lambda pair"].apply(
            lambda s: s.split("–")[1] + "–" + s.split("–")[0] if "–" in s else s
        )

    # Keep only desired columns.
    df = df[["Lambda pair", "MBAR", "MBAR_Error", "BAR", "BAR_Error", "TI", "TI_Error"]]

    # Debug: print dataframe summary (can be removed later)
    print(df.head())
    return df


def plot_results_from_csv(file_path, invert_TI=True, group_size=10):
    """
    Parse the free energy results CSV file and produce a grouped bar plot.

    Parameters:
      file_path (str): Path to the CSV file.
      invert_TI (bool): Whether to invert TI values (and their order).
      group_size (int): Maximum number of lambda pairs per row.

    Returns:
      matplotlib.figure.Figure: The generated figure.
    """
    df = parse_results_csv(file_path, invert_TI=invert_TI)

    n = len(df)
    # Group rows for plotting
    groups = []
    for start in range(0, n, group_size):
        groups.append((start, min(start + group_size, n)))
    n_groups = len(groups)

    # Create subplots vertically
    fig, axes = plt.subplots(nrows=n_groups, ncols=1, figsize=(7, 7))
    if n_groups == 1:
        axes = [axes]

    bar_width = 0.25
    for i, (start, end) in enumerate(groups):
        ax = axes[i]
        sub_df = df.iloc[start:end].copy()
        # If the group is smaller than group_size, pad with empty rows.
        if len(sub_df) < group_size:
            pad_count = group_size - len(sub_df)
            pad_data = {
                "Lambda pair": [""] * pad_count,
                "MBAR": [np.nan] * pad_count,
                "MBAR_Error": [0] * pad_count,
                "BAR": [np.nan] * pad_count,
                "BAR_Error": [0] * pad_count,
                "TI": [np.nan] * pad_count,
                "TI_Error": [0] * pad_count,
            }
            pad_df = pd.DataFrame(pad_data)
            sub_df = pd.concat([sub_df, pad_df], ignore_index=True)

        x = np.arange(group_size)
        ax.bar(
            x - bar_width,
            sub_df["MBAR"],
            bar_width,
            yerr=sub_df["MBAR_Error"],
            capsize=3,
            label="MBAR",
            color="#C45AEC",
        )
        ax.bar(
            x + bar_width,
            sub_df["BAR"],
            bar_width,
            yerr=sub_df["BAR_Error"],
            capsize=3,
            label="BAR",
            color="#F9B7FF",
        )
        ax.bar(
            x,
            sub_df["TI"],
            bar_width,
            yerr=sub_df["TI_Error"],
            capsize=3,
            label="TI",
            color="#6698FF",
        )

        ax.set_xticks(x)
        ax.set_xticklabels(sub_df["Lambda pair"], rotation=30, fontsize=14)
        ax.spines["bottom"].set_visible(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(
            axis="x", which="both", bottom=False, top=False, labelbottom=True
        )
        if i == 0:
            ax.legend(fontsize=14, frameon=False)

    fig.supylabel("ΔG (kcal/mol)", fontsize=20)
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    out_file = os.path.join(os.getcwd(), "breakdown.pdf")
    print(f"Saving figure to {out_file}")
    plt.savefig(out_file)
    # plt.show()  # Uncomment to display interactively

    return fig


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Plot free energy summary from a CSV file in a given directory."
    )
    parser.add_argument(
        "directory",
        help="Directory where the CSV file (free_energy_summary.csv) is located.",
    )
    parser.add_argument(
        "--csv",
        default="free_energy_summary.csv",
        help="Name of the CSV file (default: free_energy_summary.csv)",
    )
    parser.add_argument(
        "--invert-ti", action="store_true", help="Invert TI values (default: False)"
    )
    parser.add_argument(
        "--group-size",
        type=int,
        default=10,
        help="Number of lambda pairs per subplot (default: 10)",
    )
    args = parser.parse_args()

    # Resolve directory
    input_dir = Path(args.directory).resolve()
    if not input_dir.is_dir():
        print(f"Error: {input_dir} is not a valid directory.")
        sys.exit(1)

    csv_file = input_dir / args.csv
    if not csv_file.is_file():
        print(f"Error: CSV file {csv_file} not found in {input_dir}.")
        sys.exit(1)

    # Determine invert_TI based on command-line flag or content in directory name
    invert_TI = args.invert_ti or ("stateB" in str(input_dir).lower())

    print(f"Processing CSV file {csv_file} with invert_TI = {invert_TI}")
    plot_results_from_csv(
        str(csv_file), invert_TI=invert_TI, group_size=args.group_size
    )


if __name__ == "__main__":
    main()
