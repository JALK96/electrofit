import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import seaborn as sns

sns.set_context("talk", font_scale=0.9)

def parse_results_csv(file_path, invert_TI=True):
    """
    Parse a CSV file that contains a results table and extract the data.
    
    The CSV file is expected to have the following columns:
      state, name, MBAR, MBAR_Error, BAR, BAR_Error, TI, TI_Error
    
    The 'name' column represents the lambda pair (e.g. "0 -- 1") and will be reformatted 
    with an en-dash (–). Rows with a 'state' value of "Stages" or "TOTAL" are ignored.
    
    Parameters
    ----------
    file_path : str
        Path to the CSV file.
    invert_TI : bool, optional
        If True, multiply the TI values by -1 and invert the lambda pair order.
    
    Returns
    -------
    pd.DataFrame
        A DataFrame with columns: Lambda pair, MBAR, MBAR_Error, BAR, BAR_Error, TI, TI_Error.
    """
    # Read the CSV file.
    df = pd.read_csv(file_path)
    
    # Filter out rows where the first column indicates "Stages" or "TOTAL"
    # (Note: The CSV has a column "state" that has "States" for data and "Stages" for summary rows).
    df = df[~df['state'].str.lower().isin(["stages", "total"])]
    
    # Convert the 'name' column into the lambda pair column with an en-dash
    df.rename(columns={"name": "Lambda pair"}, inplace=True)
    df["Lambda pair"] = df["Lambda pair"].str.replace("--", "–").str.strip()
    
    # Multiply TI values by -1 if requested.
    if invert_TI:
        df["TI"] = -df["TI"]
    
        # Invert the order of the rows.
        df["TI"] = df["TI"].iloc[::-1].reset_index(drop=True)
        df["TI_Error"] = df["TI_Error"].iloc[::-1].reset_index(drop=True)

        df["Lambda pair"] = df["Lambda pair"].iloc[::-1].reset_index(drop=True)
        df["Lambda pair"] = df["Lambda pair"].apply(
            lambda s: s.split("–")[1] + "–" + s.split("–")[0] if "–" in s else s
        )
    
    # Ensure the DataFrame contains only the desired columns.
    df = df[["Lambda pair", "MBAR", "MBAR_Error", "BAR", "BAR_Error", "TI", "TI_Error"]]
    
    print(df)
    return df

def plot_results_from_csv(file_path, invert_TI=True, group_size=10):
    """
    Parse the results from a CSV file and produce a grouped bar plot.
    
    Parameters
    ----------
    file_path : str
        Path to the CSV file containing the results table.
    invert_TI : bool, optional
        If True, the TI values are multiplied by -1 and the lambda pair order is inverted.
    group_size : int, optional
        Number of lambda pairs to display per subplot row.
    """
    df = parse_results_csv(file_path, invert_TI=invert_TI)
    
    n = len(df)
    # Create groups of size 'group_size'
    groups = []
    for start in range(0, n, group_size):
        groups.append((start, min(start+group_size, n)))
    n_groups = len(groups)
    
    # Create vertically stacked subplots (same aesthetics as your original code)
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
                'Lambda pair': [''] * pad_count,
                'MBAR': [np.nan] * pad_count,
                'MBAR_Error': [0] * pad_count,
                'BAR': [np.nan] * pad_count,
                'BAR_Error': [0] * pad_count,
                'TI': [np.nan] * pad_count,
                'TI_Error': [0] * pad_count
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
    # plt.show()
    
    return fig

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


paths = [
    "/home/johannal96/MA/electrofit/FEP/equilibrium-FEP/intermediat/edge_011111-101111/stateB/run1/analysis",
]

for path in paths:
    os.chdir(path)
    csv_file = "free_energy_summary.csv"
    if "stateB" in path:
        plot_results_from_csv(csv_file, invert_TI=True)
    else:
        plot_results_from_csv(csv_file, invert_TI=False)