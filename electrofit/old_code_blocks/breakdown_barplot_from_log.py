import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
sns.set_context("talk", font_scale=0.9)

import os



def parse_results_log(file_path, invert_TI=True):
    """
    Parse a log file that contains a results table and extract the data.
    
    The function searches for a header line containing "MBAR", "BAR", and "TI" and then
    parses subsequent lines that either start with "States" or with a numeric token (for rows
    after the first) that indicate a lambda pair. Lines starting with "Stages" or "TOTAL" are ignored.
    
    Parameters
    ----------
    file_path : str
        Path to the log file.
    invert_TI : bool, optional
        If True, multiply the TI values by -1. (Default is True)
    
    Returns
    -------
    pd.DataFrame
        A DataFrame with columns: Lambda pair, MBAR, MBAR_Error, BAR, BAR_Error, TI, TI_Error.
        The Lambda pair is formatted with an en-dash (–) between the two numbers.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # Find the header line of the table by looking for a line that contains all three keywords.
    table_start = None
    for i, line in enumerate(lines):
        if "MBAR" in line and "BAR" in line and "TI" in line:
            table_start = i
            break
    if table_start is None:
        raise ValueError("Could not find the table header in the log file.")
    
    # Process subsequent lines to extract data rows.
    data_rows = []
    for line in lines[table_start+1:]:
        if not line.strip():
            continue
        parts = line.strip().split()
        # Skip rows that start with unwanted tokens (e.g. "Stages" or "TOTAL")
        if parts[0].lower() in ["stages", "total"]:
            continue
        
        # Check if the row begins with "States"
        if parts[0].lower() == "states":
            # Expected format: "States", state1, "--", state2, value1, value2, ..., value6
            if len(parts) < 8:
                continue
            lambda_pair = parts[1] + "–" + parts[3]  # using an en-dash
            numeric_parts = parts[4:]
        # Otherwise, if the first token is numeric, assume it's a continuation row.
        elif parts[0].replace('.', '', 1).isdigit():
            # Expected format: state1, "--", state2, value1, value2, ..., value6
            if len(parts) < 7:
                continue
            lambda_pair = parts[0] + "–" + parts[2]
            numeric_parts = parts[3:]
        else:
            continue

        if len(numeric_parts) < 6:
            continue
        
        try:
            mbar = float(numeric_parts[0])
            mbar_err = float(numeric_parts[1])
            bar = float(numeric_parts[2])
            bar_err = float(numeric_parts[3])
            ti = float(numeric_parts[4])
            ti_err = float(numeric_parts[5])
        except ValueError:
            continue
        
        # Multiply TI by -1 if requested (for instance, if TI needs to be inverted for plotting)
        ti_val = -ti if invert_TI else ti
        
        data_rows.append({
            "Lambda pair": lambda_pair,
            "MBAR": mbar,
            "MBAR_Error": mbar_err,
            "BAR": bar,
            "BAR_Error": bar_err,
            "TI": ti_val,
            "TI_Error": ti_err
        })
    
    if not data_rows:
        raise ValueError("No data rows were parsed from the log file.")
    
    df = pd.DataFrame(data_rows)
        # Invert the order of the TI and TI_Error columns if requested
    if invert_TI:
        df["TI"] = df["TI"].iloc[::-1].reset_index(drop=True)
        df["TI_Error"] = df["TI_Error"].iloc[::-1].reset_index(drop=True)
        df["Lambda pair"] = df["Lambda pair"].iloc[::-1].reset_index(drop=True)
        df["Lambda pair"] = df["Lambda pair"].apply(
            lambda s: s.split("–")[1] + "–" + s.split("–")[0] if "–" in s else s
        )
    print(df)
    return df

def plot_results_from_log(file_path, invert_TI=True, group_size=10):
    """
    Parse the results from a log file and produce a grouped bar plot.
    
    Parameters
    ----------
    file_path : str
        Path to the log file containing the results table.
    invert_TI : bool, optional
        If True, the TI values are multiplied by -1. (Default is True)
    group_size : int, optional
        Number of lambda pairs to display per subplot row. (Default is 10)
    """
    df = parse_results_log(file_path, invert_TI=invert_TI)
    
    n = len(df)
    # Create groups of size 'group_size'
    groups = []
    for start in range(0, n, group_size):
        groups.append((start, min(start+group_size, n)))
    n_groups = len(groups)
    
    # Create vertically stacked subplots
    fig, axes = plt.subplots(nrows=n_groups, ncols=1, figsize=(7, 7))
    if n_groups == 1:
        axes = [axes]
    
    bar_width = 0.25
    for i, (start, end) in enumerate(groups):
        ax = axes[i]
        sub_df = df.iloc[start:end].copy()
        
        # If the group is not full, pad with empty rows to keep consistent alignment
        if len(sub_df) < group_size:
            pad_count = group_size - len(sub_df)
            pad_data = {
                'MBAR': [np.nan] * pad_count,
                'MBAR_Error': [0] * pad_count,
                'BAR': [np.nan] * pad_count,
                'BAR_Error': [0] * pad_count,
                'TI': [np.nan] * pad_count,
                'TI_Error': [0] * pad_count,
                'Lambda pair': [''] * pad_count
            }
            pad_df = pd.DataFrame(pad_data)
            sub_df = pd.concat([sub_df, pad_df], ignore_index=True)
        
        x = np.arange(group_size)
        
        # Plot the three sets of bars: MBAR, BAR, and TI (already inverted if requested)
        ax.bar(x - bar_width, sub_df['MBAR'], bar_width,
               yerr=sub_df['MBAR_Error'], capsize=3, label='MBAR')
        ax.bar(x, sub_df['BAR'], bar_width,
               yerr=sub_df['BAR_Error'], capsize=3, label='BAR')
        ax.bar(x + bar_width, sub_df['TI'], bar_width,
               yerr=sub_df['TI_Error'], capsize=3, label='TI')
        
        # Set the x-ticks to the lambda pair labels (empty labels for padded entries)
        ax.set_xticks(x)
        ax.set_xticklabels(sub_df['Lambda pair'], rotation=30)
        
        # Hide the top, bottom, and right spines while leaving x-tick labels visible
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
        
        # Only add the legend to the first subplot to avoid repetition.
        if i == 0:
            ax.legend(fontsize=14, frameon=False)
    
    # Add one common y-axis label on the left side
    fig.supylabel('ΔG (kcal/mol)', fontsize=20)
    plt.tight_layout(rect=[0, 0, 1, 0.94])
    plt.savefig("breakdown.pdf")
    #plt.show()


paths = [
    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/011111-101111/edges/edge_011111_101111/water/analyse_r1/stateA",
    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/011111-101111/edges/edge_011111_101111/water/analyse_r1/stateB",
    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/011111-111011/edges/edge_011111_111011/water/analyse_r1/stateA",
    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/011111-111011/edges/edge_011111_111011/water/analyse_r1/stateB",
    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/011111-111101/edges/edge_011111_111101/water/analyse_r1/stateA",
    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/011111-111101/edges/edge_011111_111101/water/analyse_r1/stateB",


    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/101111-111011_com_cos/edges/edge_101111_111011/water/analyse_r1/stateA",
    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/101111-111011_com_cos/edges/edge_101111_111011/water/analyse_r1/stateB",

    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/101111-111101_com_cos/edges/edge_101111_111101/water/analyse_r1/stateA",
    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/101111-111101_com_cos/edges/edge_101111_111101/water/analyse_r1/stateB",

    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/111011-111101_com_cos/edges/edge_111011_111101/water/analyse_r1/stateA",
    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/111011-111101_com_cos/edges/edge_111011_111101/water/analyse_r1/stateB",

    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/3ns/011111-101111/workpath/edge_011111_101111/water/analyse_r1/stateA",
    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/3ns/011111-101111/workpath/edge_011111_101111/water/analyse_r1/stateB",
]
for path in paths:
    # Change directory to the current analysis directory.
    os.chdir(path)
    
    # Assume the log file is named "analysis.log" within each directory.
    log_file = "analysis.log"
    
    # If the path contains "stateB", invert the TI data.
    if "stateB" in path:
        plot_results_from_log(log_file, invert_TI=True)
    else:
        plot_results_from_log(log_file, invert_TI=False)