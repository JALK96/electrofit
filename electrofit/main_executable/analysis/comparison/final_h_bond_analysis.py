#!/usr/bin/env python3
"""
Combined Script for Intra Hydrogen Bond Analysis and Oxygen Occurrence Mapping
with Plotting Logic Matching Your Provided Snippet for the Oxygen Summary

Generates:
 1) A 3-column PDF showing lifetimes (box), distances (box), and H-bond counts (violin).
 2) A summary PDF with horizontal bars (occurrence times) and a heatmap (existence fraction),
    matching the axis directions and ordering you specified in your snippet,
    AND using the same color scheme as the 3-column figure.

Logs: analysis.log
Author: Your Name
Date:   2025-01-27
"""

import os
import re
import sys
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

logging.basicConfig(
    filename='analysis.log',
    filemode='w',  # overwrite on each run
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)

sns.set_context("talk")

# ---------------------------------------------------------------------------
# Common / Shared Utilities
# ---------------------------------------------------------------------------

def refine_atom_name(atom):
    """
    Refines the atom name based on specified rules:
      - If the atom is named 'O', change it to 'O1'.
      - If the atom is named 'O<number>', increment the number by 1 (e.g., O10 -> O11).
    """
    if re.fullmatch(r'[A-Za-z]+', atom):
        return atom + '1'
    match = re.fullmatch(r'([A-Za-z]+)(\d+)', atom)
    if match:
        name = match.group(1)
        number = int(match.group(2)) + 1
        return f"{name}{number}"
    return atom

def find_project_root(current_dir, project_name="electrofit"):
    """
    Find the root directory named `project_name` by ascending from `current_dir`.
    Returns the outermost match or raises FileNotFoundError if not found.
    """
    logger.info(f"Locating project root '{project_name}' from: {current_dir}")
    root = None
    while True:
        parent = os.path.dirname(current_dir)
        if os.path.basename(current_dir) == project_name:
            root = current_dir
        if parent == current_dir:  # reached filesystem root
            if root is None:
                raise FileNotFoundError(f"Project root '{project_name}' not found.")
            logger.info(f"Found project root: {root}")
            return root
        current_dir = parent

# ---------------------------------------------------------------------------
# Functions For 3-Column (Box/Violin) Plots
# ---------------------------------------------------------------------------

def load_hb_num_xvg(filename):
    """
    Loads (# H-bonds vs. time) from a GROMACS .xvg file, skipping @/# lines.
    Returns a numpy array with shape (N, 2): time in first column, count in second.
    """
    logger.info(f"Loading H-bond number data from {filename}")
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('@','#')):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    floats = list(map(float, parts[:2]))
                    data.append(floats)
                except ValueError:
                    pass
    arr = np.array(data)
    logger.info(f"H-bond number data shape: {arr.shape}")
    return arr

def parse_xpm(file_path):
    """
    Parses an XPM file and converts it to a binary NumPy array.

    Assumes that '#FF0000' in the file indicates a "present bond" (1),
    and everything else => 0.

    Returns
    -------
    data_matrix : np.ndarray
        2D array of shape (height, width)
    metadata : dict
        Parsed metadata (title etc.) from commented lines
    """
    logger.info(f"Parsing XPM file: {file_path}")
    with open(file_path, 'r') as file:
        lines = file.readlines()

    metadata = {}
    data_lines = []
    color_map = {}
    header_found = False
    data_started = False

    # Initialize placeholders. Will be assigned once we parse the header line:
    width = height = num_colors = chars_per_pixel = None

    for line_index, line_raw in enumerate(lines, start=1):
        line = line_raw.strip()

        # 1) Grab metadata from lines like:
        # /* title: "XYZ" */
        if line.startswith("/*") and not data_started:
            comment = line.strip("/* ").strip(" */")
            if ':' in comment:
                key, value = comment.split(":", 1)
                metadata[key.strip().lower()] = value.strip().strip('"')
            continue

        # 2) Skip 'static char' line
        if line.startswith('static char'):
            continue

        # 3) Look for the header line containing width, height, num_colors, chars_per_pixel
        #    Typically looks like: "64  64    10  1",
        if (not header_found) and line.startswith('"'):
            # Safely remove leading/trailing quote(s) and comma
            # e.g.  "64 64 10 1",
            raw_header = line.strip().rstrip(',')
            if raw_header.startswith('"'):
                raw_header = raw_header[1:]
            if raw_header.endswith('"'):
                raw_header = raw_header[:-1]

            tokens = raw_header.split()
            logger.debug(f"Header tokens = {tokens}")

            if len(tokens) >= 4:
                try:
                    width, height, num_colors, chars_per_pixel = map(int, tokens[:4])
                    header_found = True
                    logger.info(f"XPM header => width={width}, height={height}, "
                                f"colors={num_colors}, chars/pixel={chars_per_pixel}")
                    continue
                except ValueError:
                    logger.error(f"Error parsing header line {line_index}: {line_raw}")
                    raise
            else:
                logger.error(f"Header line {line_index} doesn't have enough tokens: {line_raw}")
                raise ValueError(f"Invalid XPM header line: {line_raw}")
        
        # 4) If we've parsed the header but not started data, we must parse color definitions
        if header_found and not data_started:
            cdef = line.strip().rstrip(',')
            if cdef.startswith('"'):
                cdef = cdef[1:]
            if cdef.endswith('"'):
                cdef = cdef[:-1]

            # Example pattern: (chars_per_pixel) c #XYZ
            pattern = rf'(.{{{chars_per_pixel}}})\s+c\s+(\S+)'
            match = re.match(pattern, cdef)
            if match:
                symbol, color_val = match.groups()
                color_map[symbol] = color_val
            if len(color_map) == num_colors:
                data_started = True
            continue

        # 5) Reading the actual pixel rows
        if data_started and line.startswith('"'):
            row_data = line.strip().rstrip(',')
            if row_data.startswith('"'):
                row_data = row_data[1:]
            if row_data.endswith('"'):
                row_data = row_data[:-1]
            data_lines.append(row_data)

    if width is None or height is None:
        logger.error("Could not find a valid XPM header in the file. Aborting.")
        raise ValueError("No valid XPM header found.")
    if len(data_lines) != height:
        logger.warning(f"Expected {height} data lines, found {len(data_lines)} in {file_path}.")

    # 6) Build the binary matrix
    data_matrix = np.zeros((height, width), dtype=int)
    for y, row_str in enumerate(data_lines):
        for x, char in enumerate(row_str):
            color = color_map.get(char, None)
            if color == '#FF0000':
                data_matrix[y, x] = 1
            else:
                data_matrix[y, x] = 0

    logger.info(f"Parsed XPM: final matrix shape = {data_matrix.shape}")
    return data_matrix, metadata

def analyze_hydrogen_bonds(data_matrix, metadata):
    """
    For a binary matrix (rows=bonds, columns=frames),
    returns dict of: hbonds_over_time, hbonds_per_index, lifetimes.
    """
    logger.info("Analyzing hydrogen bonds from the matrix.")
    hbonds_over_time = np.sum(data_matrix, axis=0)
    hbonds_per_index = np.sum(data_matrix, axis=1)

    # Lifetimes
    lifetimes = []
    for row in data_matrix:
        current = 0
        bond_lifetimes = []
        for state in row:
            if state == 1:
                current += 1
            else:
                if current>0:
                    bond_lifetimes.append(current)
                    current = 0
        if current>0:
            bond_lifetimes.append(current)
        lifetimes.append(bond_lifetimes)

    return {
        'hbonds_over_time': hbonds_over_time,
        'hbonds_per_index': hbonds_per_index,
        'lifetimes': lifetimes
    }

def parse_hbond_log_to_dataframe(file_path):
    """
    Parses a GROMACS .log to retrieve donor-acceptor pairs with refined atom names.
    Returns a DataFrame [idx, donor, acceptor].
    """
    logger.info(f"Parsing H-bond log: {file_path}")
    hbond_pairs = []
    line_pattern = re.compile(r'^\s*(\S+)\s+-\s+(\S+)\s*$')
    atom_pattern = re.compile(r'^[A-Za-z]+\d+([A-Za-z]+\d*)$')

    with open(file_path, 'r') as f:
        for line_number, line in enumerate(f, start=1):
            line=line.strip()
            if (not line) or line.startswith('#') or line.startswith('"""') or line.startswith('*'):
                continue
            match = line_pattern.match(line)
            if match:
                donor_full, acceptor_full = match.groups()
                d_match = atom_pattern.match(donor_full)
                a_match = atom_pattern.match(acceptor_full)
                if d_match and a_match:
                    donor_atom = refine_atom_name(d_match.group(1))
                    acceptor_atom = refine_atom_name(a_match.group(1))
                    hbond_pairs.append({'donor': donor_atom, 'acceptor': acceptor_atom})
                else:
                    logger.warning(f"Line {line_number}: Could not parse donor/acceptor => {line}")
            else:
                logger.warning(f"Line {line_number}: Did not match pattern => {line}")

    df = pd.DataFrame(hbond_pairs)
    df.reset_index(inplace=True)
    df.rename(columns={'index':'idx'}, inplace=True)
    df['idx'] = df.index
    logger.info(f"Found {len(df)} bonds in log.")
    return df

def one_row_three_columns_box_violin(
    df_lifetimes,    # columns: [species, num_ones, lifetime_ns]
    df_distance,     # columns: [species, num_ones, distance_nm, frequency]
    df_hbonds,       # columns: [species, num_ones, hbond_count]
    color_dict=None, # e.g. { "species_id": "#RRGGBB" }
    output_prefix="intra_hbonds",
    figure_size=(20, 8),
    spacer_prefix="~space_"
):
    """
    Creates a single figure with 3 columns:
      1) Lifetime (Box),
      2) Distance (Box, expanded by frequency),
      3) H-bond Count (Violin).
    Ensures grouping by 5->4->3 (with spacer rows).
    Saves the figure as PDF with the given prefix.
    """
    logger.info("Generating 3-column figure (lifetime, distance, hbond_count)...")

    # Check empties
    if (df_lifetimes is None or df_lifetimes.empty) and \
       (df_distance is None or df_distance.empty) and \
       (df_hbonds is None or df_hbonds.empty):
        logger.warning("All input DataFrames empty. Skipping 3-col plot.")
        return

    # Collect species in any DF
    species_in_any = set()
    for df in [df_lifetimes, df_distance, df_hbonds]:
        if df is not None and not df.empty:
            species_in_any.update(df["species"].unique())

    # Filter only species that have 3,4,5 ones
    def get_num_ones(sp, dfs):
        for d in dfs:
            if d is not None and not d.empty:
                row = d.loc[d["species"] == sp]
                if not row.empty:
                    return row.iloc[0]["num_ones"]
        return None

    all_with_num = []
    for sp in species_in_any:
        val = get_num_ones(sp, [df_lifetimes, df_distance, df_hbonds])
        if val in [3,4,5]:
            all_with_num.append((sp, val))

    if not all_with_num:
        logger.warning("No species with num_ones in [3,4,5]. Skipping 3-col plot.")
        return

    # Sort so 5->top, 4->middle, 3->bottom, alpha tiebreak
    all_with_num.sort(key=lambda x: (-x[1], x[0]))

    # Insert spacer labels
    grouped_species_order = []
    prev_n_ones = None
    for i, (sp, n_ones) in enumerate(all_with_num):
        if i>0 and n_ones != prev_n_ones:
            grouped_species_order.append((f"{spacer_prefix}{prev_n_ones}to{n_ones}", None))
        grouped_species_order.append((sp, n_ones))
        prev_n_ones = n_ones

    full_order = [t[0] for t in grouped_species_order]

    # Expand distance by frequency
    if df_distance is not None and not df_distance.empty:
        df_distance = df_distance.copy()
        df_distance['frequency'] = df_distance['frequency'].astype(int)
        df_dist_expanded = df_distance.loc[
            df_distance.index.repeat(df_distance['frequency'])
        ].reset_index(drop=True)
    else:
        df_dist_expanded = pd.DataFrame()

    def inject_spacers(df, measure_col):
        """Add spacer rows for dummy species so Seaborn sees them but draws blank."""
        if df is None or df.empty:
            return pd.DataFrame(columns=["species", measure_col, "num_ones"])
        new_df = df.copy()
        spacer_rows = []
        for sp, val in grouped_species_order:
            if val is None:
                spacer_rows.append({
                    "species": sp,
                    measure_col: np.nan,
                    "num_ones": -1
                })
        if spacer_rows:
            spacer_df = pd.DataFrame(spacer_rows)
            combined = pd.concat([new_df, spacer_df], ignore_index=True)
            return combined
        return new_df

    df_life2 = inject_spacers(df_lifetimes, "lifetime_ns")
    df_dist2 = inject_spacers(df_dist_expanded, "distance_nm")
    df_hbond2 = inject_spacers(df_hbonds, "hbond_count")

    # Ensure color_dict has entries for spacer
    if color_dict is None:
        color_dict = {}
    else:
        color_dict = color_dict.copy()

    spacer_color = "#D3D3D3"
    for sp, val in grouped_species_order:
        if val is None and sp not in color_dict:
            color_dict[sp] = spacer_color

    # Create figure
    fig, axes = plt.subplots(1,3, figsize=figure_size)
    ax_life, ax_dist, ax_hbond = axes

    # Plot 1: lifetimes
    if df_life2 is not None and not df_life2.empty:
        sns.boxplot(
            data=df_life2, x="lifetime_ns", y="species",
            order=full_order, orient='h', palette=color_dict,
            showmeans=True, showfliers=False,
            meanprops={
                'marker':'o','markerfacecolor':'black','markeredgecolor':'black','markersize':8
            },
            ax=ax_life
        )
        ax_life.set_title("Lifetime (Box)")
        ax_life.set_xlabel("Lifetime (ns)")
        ax_life.set_ylabel("Species")

        # Hide spacer labels
        new_labels = []
        for lbl in ax_life.get_yticklabels():
            txt = lbl.get_text()
            if txt.startswith(spacer_prefix): new_labels.append("")
            else: new_labels.append(txt)
        ax_life.set_yticklabels(new_labels)
    else:
        ax_life.set_title("No Lifetime Data")
        ax_life.set_xlabel("")
        ax_life.set_ylabel("")

    # Plot 2: distance
    if not df_dist2.empty:
        sns.boxplot(
            data=df_dist2, x="distance_nm", y="species",
            order=full_order, orient='h', palette=color_dict,
            showmeans=True, showfliers=False,
            meanprops={
                'marker':'o','markerfacecolor':'black','markeredgecolor':'black','markersize':8
            },
            ax=ax_dist
        )
        ax_dist.set_title("Distance (Box)")
        ax_dist.set_xlabel("Distance (nm)")
        ax_dist.set_ylabel("")
        ax_dist.set_yticklabels([])
    else:
        ax_dist.set_title("No Distance Data")
        ax_dist.set_xlabel("")
        ax_dist.set_ylabel("")
        ax_dist.set_yticklabels([])

    # Plot 3: hbond count
    if df_hbond2 is not None and not df_hbond2.empty:
        sns.violinplot(
            data=df_hbond2, x="hbond_count", y="species",
            order=full_order, orient='h', palette=color_dict,
            inner=None, cut=0, ax=ax_hbond
        )
        ax_hbond.set_title("H-bond Count (Violin)")
        ax_hbond.set_xlabel("H-bond Count")
        ax_hbond.set_ylabel("")
        ax_hbond.set_yticklabels([])

        # overlay means
        df_means = df_hbond2.groupby('species')['hbond_count'].mean().reset_index()
        # exclude spacer
        df_means = df_means[~df_means['species'].str.startswith(spacer_prefix)]
        ax_hbond.scatter(
            df_means['hbond_count'], df_means['species'],
            color='black', marker='o', s=80, label='Mean', zorder=5
        )
        ax_hbond.legend(loc='lower right')
    else:
        ax_hbond.set_title("No H-bond Data")
        ax_hbond.set_xlabel("")
        ax_hbond.set_ylabel("")
        ax_hbond.set_yticklabels([])

    plt.tight_layout()
    out_fname = f"{output_prefix}_3col_naive_space.pdf"
    plt.savefig(out_fname, dpi=300)
    plt.close()
    logger.info(f"Saved 3-column figure to '{out_fname}'")


# ---------------------------------------------------------------------------
# Functions For Oxygen Occurrence Summary (Matching the Snippet's Logic)
# ---------------------------------------------------------------------------

def count_oxygen_occurrences_from_matrix(hbond_df, hbonds_per_index):
    """
    Based on hbond_df['idx'] -> row index in data_matrix,
    sums total occurrences (hbonds_per_index) for each oxygen.
    """
    logger.info("Counting oxygen occurrences from the matrix.")
    # If 'count' already exists, remove or rename to avoid confusion
    if 'count' in hbond_df.columns:
        hbond_df = hbond_df.drop(columns=['count'])
    hbond_df['count'] = hbonds_per_index[hbond_df['idx'].values]

    melted = hbond_df.melt(
        id_vars=['idx','count'], value_vars=['donor','acceptor'],
        value_name='oxygen'
    ).dropna(subset=['oxygen'])
    # Sum over oxygen
    counts = melted.groupby('oxygen')['count'].sum()
    logger.info(f"Found {len(counts)} oxygen atoms in the occurrences.")
    return counts

def generate_summary_plot(
    occurrence_data, 
    existence_data,
    color_dict=None,  # <== we pass the same color_dict from main
    time_per_frame_ns=0.01,
    output_file='oxygen_occurrences_summary.pdf',
    folder_order=None
):
    """
    EXACT LOGIC from your snippet. Produces bar plots of oxygen occurrence (time) on the left,
    and existence heatmap on the right. The x-axis is reversed to match your snippet.
    Now uses color_dict[species_id] for the bar color, ensuring consistency with 3-col figure.

    occurrence_data : dict{ species_id -> pd.Series(index=oxygen, data=counts) }
    existence_data  : dict{ species_id -> 2D np.array (#oxygens, #bins) }
    color_dict      : dict{ species_id -> color string }
    """
    logger.info("Generating summary plot with existence maps & occurrence bars (per snippet).")
    num_species = len(occurrence_data)
    if num_species == 0:
        logger.warning("No data in occurrence_data; skipping summary plot.")
        return

    # Collect all unique oxygens to fix a universal y-axis
    all_oxygens = set()
    for series_data in occurrence_data.values():
        all_oxygens.update(series_data.index.tolist())
    # Sort O1..On by numeric portion
    all_oxygens = sorted(
        all_oxygens,
        key=lambda x: int(re.findall(r'\d+', x)[0]) if re.findall(r'\d+', x) else 0
    )

    # Possibly reorder species by folder_order, else alphabetical
    if folder_order:
        sorted_species_ids = [sp for sp in folder_order if sp in occurrence_data]
    else:
        sorted_species_ids = sorted(occurrence_data.keys())

    # Determine max occurrence => x-limit
    max_occurrence = max([ser.max() for ser in occurrence_data.values()]) if occurrence_data else 0
    x_limit = max_occurrence * time_per_frame_ns * 1.1

    # Figure
    fig_height = max(6, 5.5*num_species)
    fig, axes = plt.subplots(
        nrows=num_species, ncols=2,
        figsize=(16, fig_height),
        constrained_layout=True,
        gridspec_kw={'width_ratios':[1,1]}
    )
    if num_species==1:
        axes = [axes]

    for i, species_id in enumerate(sorted_species_ids):
        counts = occurrence_data[species_id].reindex(all_oxygens, fill_value=0)
        time_vals = counts * time_per_frame_ns

        ax_occ, ax_exist = axes[i]

        # We'll look up color from color_dict if present
        color_for_bars = color_dict.get(species_id, "gray") if color_dict else "gray"

        # Subplot 1: horizontal bar plot
        sns.barplot(
            x=time_vals.values,
            y=all_oxygens,
            ax=ax_occ,
            color=color_for_bars
        )
        ax_occ.set_title(f"H-Bond Occurrence Time - {species_id}", fontweight='bold')
        ax_occ.set_xlabel("Time (ns)")
        # Reverse x-axis
        ax_occ.set_xlim(x_limit, 0)
        # Annotate each bar
        for p in ax_occ.patches:
            width = p.get_width()
            if width>0:
                ax_occ.annotate(
                    f"{width:.2f}",
                    (width, p.get_y() + p.get_height()/2),
                    ha='right', va='center',
                    xytext=(-5,0),
                    textcoords='offset points',
                    fontsize=13
                )
        ax_occ.grid(True, linestyle='--', alpha=0.5, axis='x')
        # Hide y tick labels for the occurrence plot
        ax_occ.set_yticklabels([])
        ax_occ.yaxis.set_ticks_position('right')

        # Subplot 2: existence heatmap
        aggregated_binned_matrix = existence_data[species_id]
        logger.info(f"Species '{species_id}' existence map shape: {aggregated_binned_matrix.shape}")
        if aggregated_binned_matrix.shape[0] != len(all_oxygens):
            logger.warning(
                f"Mismatch: species '{species_id}' existence matrix has {aggregated_binned_matrix.shape[0]} rows, "
                f"but we have {len(all_oxygens)} oxygens. Skipping heatmap."
            )
            ax_exist.set_visible(False)
            continue

        sns.heatmap(
            aggregated_binned_matrix,
            cmap="Reds",
            cbar=True,
            ax=ax_exist,
            linewidths=0,
            linecolor='white',
            cbar_kws={"label":"Fraction of Bond Presence"},
            vmin=0,
            vmax=1
        )
        ax_exist.set_title(f"H-Bond Existence Map - {species_id}", fontweight='bold')
        ax_exist.set_xlabel("Time (ns)")

        # Set x-ticks based on num_bins
        num_bins = aggregated_binned_matrix.shape[1]
        bin_size_ns = 0.2  # same assumption as snippet
        time_bins = np.arange(num_bins)*bin_size_ns
        # limit # of ticks
        n_ticks = min(6, num_bins)
        tick_positions = np.linspace(0, num_bins-1, n_ticks, dtype=int)
        tick_labels = [f"{int(time_bins[pos])}" for pos in tick_positions]

        ax_exist.set_xticks(tick_positions+0.5)
        ax_exist.set_xticklabels(tick_labels, rotation=0, ha='right')

        # y-ticks => all_oxygens
        ax_exist.set_yticks(np.arange(len(all_oxygens)) + 0.5)
        ax_exist.set_yticklabels(all_oxygens, rotation=0)

        # Make spines visible
        for spine in ax_exist.spines.values():
            spine.set_visible(True)

    plt.savefig(output_file, dpi=300)
    plt.close()
    logger.info(f"Summary plot saved as '{output_file}'")

# ---------------------------------------------------------------------------
# Main Routine
# ---------------------------------------------------------------------------

def main():
    logger.info("Starting combined analysis script...")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    try:
        project_path = find_project_root(script_dir)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)

    process_dir = os.path.join(project_path, "process.nobackup")
    if not os.path.isdir(process_dir):
        logger.error(f"Not a valid directory: {process_dir}")
        sys.exit(1)
    logger.info(f"Using process directory: {process_dir}")

    # ----------- Data for the 3-col figure -----------
    time_per_frame_ns = 0.01
    hbonds_data = []
    lifetimes_data = []
    distance_data = []
    all_species_ids = []

    process_path = Path(process_dir)
    for folder_name in sorted(os.listdir(process_path)):
        if not folder_name.startswith("IP_"):
            continue
        folder_path = process_path / folder_name
        if not folder_path.is_dir():
            continue

        species_id = folder_name.replace("IP_","")
        num_ones = species_id.count('1')
        if num_ones not in [3,4,5]:
            continue

        all_species_ids.append(species_id)

        # relevant files
        hb_num_file = folder_path / "analyze_final_sim" / "h_bonds" / "intra_hb_num.xvg"
        xpm_file    = folder_path / "analyze_final_sim" / "h_bonds" / "intra_hb_matrix.xpm"
        dist_file   = folder_path / "analyze_final_sim" / "h_bonds" / "intra_hb_dist.xvg"

        if not (hb_num_file.is_file() and xpm_file.is_file() and dist_file.is_file()):
            logger.warning(f"Skipping {species_id}: missing one of the hbond files in {folder_path}")
            continue

        # (1) H-bond counts
        data_num = load_hb_num_xvg(str(hb_num_file))
        if data_num.size < 2:
            logger.warning(f"Empty hbond num data for {species_id}")
            continue
        hbond_counts = data_num[:,1]
        for val in hbond_counts:
            hbonds_data.append( (species_id, num_ones, val) )

        # (2) Lifetime
        matrix_data, meta = parse_xpm(str(xpm_file))
        analysis = analyze_hydrogen_bonds(matrix_data, meta)
        lifetime_frames = [lf for bond_lf in analysis['lifetimes'] for lf in bond_lf]
        lifetimes_ns = [f * time_per_frame_ns for f in lifetime_frames]
        for lf in lifetimes_ns:
            lifetimes_data.append( (species_id, num_ones, lf) )

        # (3) Distances
        d_list = []
        with open(dist_file, 'r') as fdist:
            for line in fdist:
                if not line.strip() or line.startswith(('#','@')):
                    continue
                parts = line.split()
                if len(parts)==2:
                    try:
                        x = float(parts[0])
                        y = float(parts[1])
                        d_list.append( (x,y) )
                    except ValueError:
                        pass
        d_list = np.array(d_list)
        if d_list.size<2:
            logger.warning(f"No valid distance data for {species_id}")
            continue

        # filter out distances < 0.2 nm
        d_list = d_list[d_list[:,0]>=0.2]
        for row in d_list:
            distance_data.append((species_id,num_ones,row[0],row[1]))

    # Build DataFrames
    df_hbonds = pd.DataFrame(hbonds_data, columns=["species","num_ones","hbond_count"])
    df_lifetimes = pd.DataFrame(lifetimes_data, columns=["species","num_ones","lifetime_ns"])
    df_distance = pd.DataFrame(distance_data, columns=["species","num_ones","distance_nm","frequency"])

    logger.info(f"Data for 3-col figure: H-bonds {df_hbonds.shape}, lifetimes {df_lifetimes.shape}, distances {df_distance.shape}")

    # assign color dict
    all_species_ids = list(set(all_species_ids))
    all_species_ids.sort()
    palette = sns.color_palette("tab10", len(all_species_ids))
    color_dict = {}
    for i, sp in enumerate(all_species_ids):
        color_dict[sp] = palette[i % len(palette)]

    # produce 3-col figure
    one_row_three_columns_box_violin(
        df_lifetimes=df_lifetimes,
        df_distance=df_distance,
        df_hbonds=df_hbonds,
        color_dict=color_dict,
        output_prefix="intra_hbonds",
        figure_size=(20,8)
    )

    # ----------- Data for the oxygen summary -----------
    occurrence_data = {}
    existence_data = {}
    all_oxygens = set()

    folder_list = []
    for folder_name in sorted(os.listdir(process_path)):
        if not folder_name.startswith("IP_"):
            continue
        species_id = folder_name.replace("IP_","")
        num_ones = species_id.count('1')
        if num_ones not in [3,4,5]:
            continue
        folder_list.append((folder_name, species_id, num_ones))

    # sort: 5->4->3 descending, then alphabetical species_id
    folder_list.sort(key=lambda x: (-x[2], x[1]))
    sorted_species_ids = [x[1] for x in folder_list]

    # first pass: gather all oxygens
    for folder_name, species_id, n_ones in folder_list:
        folder_path = process_path / folder_name
        hb_log_file = folder_path / "analyze_final_sim" / "h_bonds" / "intra_hb.log"
        if not hb_log_file.is_file():
            continue
        hbond_df = parse_hbond_log_to_dataframe(str(hb_log_file))
        if hbond_df.empty:
            continue
        all_oxygens.update(hbond_df['donor'].unique())
        all_oxygens.update(hbond_df['acceptor'].unique())

    # second pass: build occurrence + existence data
    logger.info(f"Total unique oxygens found: {len(all_oxygens)}")
    for folder_name, species_id, n_ones in folder_list:
        folder_path = process_path / folder_name
        xpm_file = folder_path / "analyze_final_sim" / "h_bonds" / "intra_hb_matrix.xpm"
        hb_log_file = folder_path / "analyze_final_sim" / "h_bonds" / "intra_hb.log"

        if not (xpm_file.is_file() and hb_log_file.is_file()):
            continue
        hbond_df = parse_hbond_log_to_dataframe(str(hb_log_file))
        if hbond_df.empty:
            continue
        data_matrix, meta = parse_xpm(str(xpm_file))
        results = analyze_hydrogen_bonds(data_matrix, meta)
        hbonds_per_index = results['hbonds_per_index']
        # count occurrences
        occ_series = count_oxygen_occurrences_from_matrix(hbond_df, hbonds_per_index)
        occurrence_data[species_id] = occ_series

        # existence data
        bin_size_ns = 0.2
        frames_per_bin = int(bin_size_ns / time_per_frame_ns)
        if frames_per_bin < 1:
            logger.warning(f"Invalid bin size for {species_id}, skipping existence map.")
            continue

        n_frames = data_matrix.shape[1]
        num_bins = n_frames // frames_per_bin
        if n_frames % frames_per_bin != 0:
            num_bins += 1
        binned_matrix = np.zeros((data_matrix.shape[0], num_bins))

        # fill binned matrix
        for i in range(num_bins):
            start_idx = i*frames_per_bin
            end_idx = min((i+1)*frames_per_bin, n_frames)
            chunk = data_matrix[:,start_idx:end_idx]
            binned_matrix[:,i] = np.mean(chunk, axis=1)

        # convert bond-centric to oxygen-centric
        oxygen_to_indices = {}
        for idx, row in hbond_df.iterrows():
            d = row['donor']
            a = row['acceptor']
            oxygen_to_indices.setdefault(d, []).append(idx)
            oxygen_to_indices.setdefault(a, []).append(idx)

        aggregated_binned_matrix = np.zeros((len(all_oxygens), num_bins))
        # Sort all_oxygens the same way in the final plot function
        oxy_sorted = sorted(
            all_oxygens,
            key=lambda x: int(re.findall(r'\d+', x)[0]) if re.findall(r'\d+', x) else 0
        )
        # Build row by row
        for o_idx, oxygen in enumerate(oxy_sorted):
            bonds_for_o = oxygen_to_indices.get(oxygen, [])
            if bonds_for_o:
                for b_idx in range(num_bins):
                    aggregated_binned_matrix[o_idx, b_idx] = np.max(binned_matrix[bonds_for_o, b_idx])

        existence_data[species_id] = aggregated_binned_matrix
        logger.info(f"Built existence map for {species_id}, shape {aggregated_binned_matrix.shape}")

    # Produce summary plot with color_dict to ensure same species colors
    if occurrence_data:
        generate_summary_plot(
            occurrence_data, 
            existence_data,
            color_dict=color_dict,  # <== pass the same color_dict we used for 3-col figure
            time_per_frame_ns=time_per_frame_ns,
            output_file='oxygen_occurrences_summary.pdf',
            folder_order=sorted_species_ids
        )
    else:
        logger.warning("No occurrence data to plot. Skipping summary plot.")

    logger.info("All analysis completed.")

if __name__ == "__main__":
    main()