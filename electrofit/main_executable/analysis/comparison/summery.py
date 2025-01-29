#!/usr/bin/env python3
"""
Combined Script for H-Bond Analysis with a 2×2 Figure:
  (Top-Left)    Lifetime (Box)
  (Top-Right)   Distance (Box)
  (Bottom-Left) Engaged Time / Atom (Bar) [Donor-based]
  (Bottom-Right)H-bond Count (Violin)

Additionally, a summary plot is produced:
 - Left: Oxygen-level occurrence bars (reversed x-axis)
 - Right: Existence map heatmap
All using the same color dictionary.

We force Times New Roman font, enlarge bar-annotation font in "Engaged Time" plot,
and color donor oxygens in red in the summary heatmap's y-ticks.

Logs: analysis.log
Author: Arthur Laux
Date:   2025-01-27
"""

import os
import re
import sys
import logging
import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# ---------------------------------------------------------------------------
# Global Logging & Font Setup
# ---------------------------------------------------------------------------
logging.basicConfig(
    filename='analysis.log',
    filemode='w',  # overwrite each run
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)


sns.set_context("talk")
# Force Times New Roman as the default font
mpl.rcParams['font.family'] = 'Times New Roman'

# ---------------------------------------------------------------------------
# Common / Shared Utilities
# ---------------------------------------------------------------------------
def refine_atom_name(atom):
    """
    Refine GROMACS atom names:
      - 'O' -> 'O1'
      - 'O10' -> 'O11'
      (Increment trailing integer by 1)
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
    Traverse upward from current_dir to find a folder named project_name.
    Returns the path if found, otherwise raises FileNotFoundError.
    """
    logger.info(f"Locating project root '{project_name}' from: {current_dir}")
    root = None
    while True:
        parent = os.path.dirname(current_dir)
        if os.path.basename(current_dir) == project_name:
            root = current_dir
        if parent == current_dir:
            # Reached filesystem root
            if root is None:
                raise FileNotFoundError(f"Project root '{project_name}' not found.")
            logger.info(f"Found project root: {root}")
            return root
        current_dir = parent


def load_hb_num_xvg(filename):
    """
    Loads (# H-bonds vs. time) from a GROMACS .xvg file,
    skipping lines starting with @ or #.
    Returns Nx2 float array: [time, hbond_count].
    """
    logger.info(f"Loading H-bond number data from: {filename}")
    data = []
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(('@', '#')):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    floats = list(map(float, parts[:2]))
                    data.append(floats)
                except ValueError:
                    pass
    arr = np.array(data)
    logger.info(f"H-bond num data shape: {arr.shape}")
    return arr


def parse_xpm(file_path):
    """
    Parses an XPM file and converts #FF0000 => 1, everything else => 0.
    Returns (data_matrix, metadata).
    """
    logger.info(f"Parsing XPM: {file_path}")
    with open(file_path, 'r') as file:
        lines = file.readlines()

    metadata = {}
    data_lines = []
    color_map = {}
    header_found = False
    data_started = False

    width = height = num_colors = chars_per_pixel = None

    for line_idx, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        # Potential metadata lines: /* title: "some text" */
        if line.startswith("/*") and not data_started:
            comment = line.strip("/* ").strip(" */")
            if ':' in comment:
                key, value = comment.split(":", 1)
                metadata[key.strip().lower()] = value.strip().strip('"')
            continue

        if line.startswith('static char'):
            continue  # skip

        if (not header_found) and line.startswith('"'):
            # parse the header line "64 64 10 1",
            line_str = line.strip(',').lstrip('"').rstrip('"')
            tokens = line_str.split()
            if len(tokens) >= 4:
                try:
                    width, height, num_colors, chars_per_pixel = map(int, tokens[:4])
                    header_found = True
                    logger.info(f"XPM header => w={width}, h={height}, c={num_colors}, cpp={chars_per_pixel}")
                    continue
                except ValueError:
                    logger.error(f"Invalid header line {line_idx}: {raw_line}")
                    raise
            else:
                logger.error(f"Header line {line_idx} missing tokens: {raw_line}")
                raise ValueError(f"Bad XPM header line: {raw_line}")

        if header_found and not data_started:
            # color definitions
            cdef = line.strip(',').lstrip('"').rstrip('"')
            pattern = rf'(.{{{chars_per_pixel}}})\s+c\s+(\S+)'
            match = re.match(pattern, cdef)
            if match:
                symbol, color_val = match.groups()
                color_map[symbol] = color_val
            if len(color_map) == num_colors:
                data_started = True
            continue

        if data_started and line.startswith('"'):
            row_str = line.strip(',').lstrip('"').rstrip('"')
            data_lines.append(row_str)

    if width is None or height is None:
        logger.error("No valid XPM header found.")
        raise ValueError("XPM header missing.")

    if len(data_lines) != height:
        logger.warning(f"Expected {height} lines, got {len(data_lines)} in {file_path}.")

    arr = np.zeros((height, width), dtype=int)
    for y, row_str in enumerate(data_lines):
        for x, ch in enumerate(row_str):
            color = color_map.get(ch, None)
            if color == '#FF0000':
                arr[y,x] = 1

    logger.info(f"XPM matrix shape: {arr.shape}")
    return arr, metadata


def analyze_hydrogen_bonds(data_matrix, metadata):
    """
    Summarize a binary matrix (rows = bond indices, columns = frames):
      - hbonds_over_time (sum along rows)
      - hbonds_per_index (sum along columns)
      - lifetimes (runs of 1's in each row)
    """
    logger.info("Analyzing hydrogen bond matrix.")
    hbonds_over_time = np.sum(data_matrix, axis=0)
    hbonds_per_index = np.sum(data_matrix, axis=1)

    lifetimes = []
    for row in data_matrix:
        cur_run = 0
        run_list = []
        for val in row:
            if val == 1:
                cur_run += 1
            else:
                if cur_run > 0:
                    run_list.append(cur_run)
                    cur_run = 0
        if cur_run>0:
            run_list.append(cur_run)
        lifetimes.append(run_list)

    return {
        'hbonds_over_time': hbonds_over_time,
        'hbonds_per_index': hbonds_per_index,
        'lifetimes': lifetimes
    }


def parse_hbond_log_to_dataframe(file_path):
    """
    GROMACS .log -> DataFrame of shape [idx, donor, acceptor].
    We refine the final portion of each name (the numeric bits).
    """
    logger.info(f"Parsing H-bond log: {file_path}")
    line_pattern = re.compile(r'^\s*(\S+)\s+-\s+(\S+)\s*$')
    atom_pattern = re.compile(r'^[A-Za-z]+\d+([A-Za-z]+\d*)$')

    pairs = []
    with open(file_path, 'r') as f:
        for line_no, line_raw in enumerate(f, start=1):
            line = line_raw.strip()
            if not line or line.startswith('#') or line.startswith('"""') or line.startswith('*'):
                continue
            match = line_pattern.match(line)
            if match:
                donor_full, acceptor_full = match.groups()
                d_match = atom_pattern.match(donor_full)
                a_match = atom_pattern.match(acceptor_full)
                if d_match and a_match:
                    donor_atom = refine_atom_name(d_match.group(1))
                    acceptor_atom = refine_atom_name(a_match.group(1))
                    pairs.append({'donor': donor_atom, 'acceptor': acceptor_atom})
                else:
                    logger.warning(f"Line {line_no}: cannot parse => {line}")
            else:
                logger.warning(f"Line {line_no}: no match => {line}")

    df = pd.DataFrame(pairs)
    df.reset_index(inplace=True)
    df.rename(columns={'index':'idx'}, inplace=True)
    df['idx'] = df.index
    logger.info(f"Found {len(df)} H-bond pairs in log.")
    return df


# ------------------------------------------------------------------------------
# Summation of Donor-only occurrences
# ------------------------------------------------------------------------------

def count_donor_occurrences_from_matrix(hbond_df, hbonds_per_index):
    """
    Sums total frames for *donor* oxygens only (ignoring acceptors).
    This avoids double-counting a single H-bond from both sides.

    Steps:
      1) For each row of hbond_df => bond index => hbonds_per_index (frames).
      2) Only meltdown 'donor' atoms, ignoring 'acceptor'.
      3) Group by 'donor' oxygen => sum the frames.
      4) Return a Series of shape [oxygen -> frames].
    """
    logger.info("Counting DONOR-only occurrences from the matrix.")
    df = hbond_df.copy()
    # If 'count' already exists, remove it to avoid collisions
    if 'count' in df.columns:
        df = df.drop(columns=['count'])

    # Put the sum of frames (hbonds_per_index) into a new 'count' col
    df['count'] = hbonds_per_index[df['idx'].values]

    # Melt only the donor column
    melted = df.melt(
        id_vars=['idx','count'],
        value_vars=['donor'],  # <--- ignoring 'acceptor'
        value_name='oxygen'
    ).dropna(subset=['oxygen'])

    # Now sum over the donors
    counts = melted.groupby('oxygen')['count'].sum()
    logger.info(f"Found {len(counts)} donor atoms with nonzero occurrences.")
    return counts


# ------------------------------------------------------------------------------
# The 2×2 Figure Layout
# ------------------------------------------------------------------------------
def two_by_two_plots_box_violin(
    df_lifetimes,   # [species, num_ones, lifetime_ns]
    df_distance,    # [species, num_ones, distance_nm, frequency]
    df_engaged,     # [species, num_ones, engaged_time_ns] (Bar, bottom-left)
    df_hbonds,      # [species, num_ones, hbond_count]     (Violin, bottom-right)
    color_dict=None,
    output_prefix="intra_hbonds",
    figure_size=(16, 12),
    spacer_prefix="~space_"
):
    """
    2×2 layout:

       (0,0) Lifetime (Box)   | (0,1) Distance (Box)
       (1,0) Engaged Time/H-Atom| (1,1) H-bond Count (Violin)

    We group species by # ones => 5->4->3, with spacer labels in each subplot.
    """
    logger.info("Generating 2×2 figure: top=Lifetime/Distance, bottom=Time/H-Atom & H-bond Count.")

    # Check if everything is empty
    all_empty = True
    for df in [df_lifetimes, df_distance, df_engaged, df_hbonds]:
        if df is not None and not df.empty:
            all_empty = False
            break
    if all_empty:
        logger.warning("All DataFrames empty; skipping 2×2 plot.")
        return

    # Collect all species in any DF
    species_in_any = set()
    for df in [df_lifetimes, df_distance, df_engaged, df_hbonds]:
        if df is not None and not df.empty:
            species_in_any.update(df["species"].unique())

    def get_num_ones(sp, df_list):
        for d in df_list:
            if d is not None and not d.empty:
                row = d.loc[d["species"] == sp]
                if not row.empty:
                    return row.iloc[0]["num_ones"]
        return None

    # We only plot species with #ones in [3,4,5]
    all_with_num = []
    for sp in species_in_any:
        val = get_num_ones(sp, [df_lifetimes, df_distance, df_engaged, df_hbonds])
        if val in [3,4,5]:
            all_with_num.append((sp, val))
    if not all_with_num:
        logger.warning("No species with num_ones in [3,4,5]. Skipping 2×2 plot.")
        return

    # Sort: 5->top, 4->middle, 3->bottom
    all_with_num.sort(key=lambda x: (-x[1], x[0]))

    # Insert spacer labels
    grouped_species_order = []
    prev_n = None
    for i, (sp, n_ones) in enumerate(all_with_num):
        if i>0 and n_ones != prev_n:
            grouped_species_order.append((f"{spacer_prefix}{prev_n}to{n_ones}", None))
        grouped_species_order.append((sp, n_ones))
        prev_n = n_ones

    full_order = [t[0] for t in grouped_species_order]

    def inject_spacers(df, measure_col):
        if df is None or df.empty:
            return pd.DataFrame(columns=["species", measure_col, "num_ones"])
        new_df = df.copy()
        row_spacers = []
        for sp, val in grouped_species_order:
            if val is None:
                row_spacers.append({
                    "species": sp,
                    measure_col: np.nan,
                    "num_ones": -1
                })
        if row_spacers:
            dummy_df = pd.DataFrame(row_spacers)
            return pd.concat([new_df, dummy_df], ignore_index=True)
        return new_df

    # Expand distance by frequency if needed
    if df_distance is not None and not df_distance.empty:
        df_distance = df_distance.copy()
        df_distance['frequency'] = df_distance['frequency'].astype(int)
        df_distance_expanded = df_distance.loc[
            df_distance.index.repeat(df_distance['frequency'])
        ].reset_index(drop=True)
    else:
        df_distance_expanded = pd.DataFrame()

    df_life2   = inject_spacers(df_lifetimes, "lifetime_ns")
    df_dist2   = inject_spacers(df_distance_expanded, "distance_nm")
    df_eng2    = inject_spacers(df_engaged, "engaged_time_ns")
    df_hbond2  = inject_spacers(df_hbonds, "hbond_count")

    # Extend color_dict for spacer rows
    if color_dict is None:
        color_dict = {}
    else:
        color_dict = color_dict.copy()
    spacer_color = "#D3D3D3"
    for sp, val in grouped_species_order:
        if val is None and sp not in color_dict:
            color_dict[sp] = spacer_color

    # Create the 2x2 figure
    fig, axes = plt.subplots(2, 2, figsize=figure_size)
    # axes[0][0] => Lifetime, axes[0][1] => Distance
    # axes[1][0] => Engaged Time, axes[1][1] => H-bond Count

    ax_life = axes[0][0]
    ax_dist = axes[0][1]
    ax_eng  = axes[1][0]  # bottom-left
    ax_hbond= axes[1][1]  # bottom-right

    # (1) Lifetime (top-left)
    if not df_life2.empty:
        sns.boxplot(
            data=df_life2,
            x="lifetime_ns",
            y="species",
            order=full_order,
            orient='h',
            palette=color_dict,
            showmeans=True,
            meanprops={'marker':'o','markerfacecolor':'black','markeredgecolor':'black','markersize':8},
            showfliers=False,
            ax=ax_life
        )
        ax_life.set_title("Lifetime (Box)")
        ax_life.set_xlabel("Lifetime (ns)")
        ax_life.set_ylabel("(micro) Protonationstate")
        # Hide spacer labels
        new_lbls = []
        for lbl in ax_life.get_yticklabels():
            txt = lbl.get_text()
            if txt.startswith(spacer_prefix):
                new_lbls.append("")
            else:
                new_lbls.append(txt)
        ax_life.set_yticklabels(new_lbls)
    else:
        ax_life.set_title("No Lifetime Data")
        ax_life.set_xlabel("")
        ax_life.set_ylabel("")

    # (2) Distance (top-right)
    if not df_dist2.empty:
        sns.boxplot(
            data=df_dist2,
            x="distance_nm",
            y="species",
            order=full_order,
            orient='h',
            palette=color_dict,
            showmeans=True,
            meanprops={'marker':'o','markerfacecolor':'black','markeredgecolor':'black','markersize':8},
            showfliers=False,
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

    # (3) Engaged Time / Atom (bottom-left)
    if not df_eng2.empty:
        sns.barplot(
            data=df_eng2,
            x="engaged_time_ns",
            y="species",
            order=full_order,
            orient='h',
            palette=color_dict,
            ax=ax_eng
        )
        ax_eng.set_title("Time / Atom (Bar)")
        ax_eng.set_xlabel("Engaged Time / Atom (ns)")
        ax_eng.set_ylabel("(micro) Protonationstate")
        ax_eng.set_yticklabels(new_lbls)

        # Larger annotation font
        for patch in ax_eng.patches:
            width = patch.get_width()
            y_center = patch.get_y() + patch.get_height()/2
            if not np.isnan(width) and width>0:
                ax_eng.annotate(
                    f"{width:.2f}",
                    (width, y_center),
                    ha='right', va='center',
                    xytext=(-5,0),
                    textcoords='offset points',
                    fontsize=14  # Enlarge
                )
    else:
        ax_eng.set_title("No Engaged-Time Data")
        ax_eng.set_xlabel("")
        ax_eng.set_ylabel("(micro) Protonationstate")
        ax_eng.set_yticklabels(new_lbls)

    # (4) H-bond Count (bottom-right)
    if not df_hbond2.empty:
        sns.violinplot(
            data=df_hbond2,
            x="hbond_count",
            y="species",
            order=full_order,
            orient='h',
            palette=color_dict,
            cut=0,
            inner=None,
            ax=ax_hbond
        )
        ax_hbond.set_title("H-bond Count (Violin)")
        ax_hbond.set_xlabel("H-bond Count")
        ax_hbond.set_ylabel("")
        ax_hbond.set_yticklabels([])

        # Means
        df_means = df_hbond2.groupby('species')['hbond_count'].mean().reset_index()
        df_means = df_means[~df_means['species'].str.startswith(spacer_prefix)]
        ax_hbond.scatter(
            df_means['hbond_count'], df_means['species'],
            color='black', marker='o', s=80,  zorder=5 # ,label='Mean'
        )
        #ax_hbond.legend(loc='lower right')
    else:
        ax_hbond.set_title("No H-bond Data")
        ax_hbond.set_xlabel("")
        ax_hbond.set_ylabel("")
        ax_hbond.set_yticklabels([])

    plt.tight_layout()
    out_pdf = f"{output_prefix}_2x2_with_engaged.pdf"
    plt.savefig(out_pdf, dpi=300)
    plt.close()
    logger.info(f"Saved 2x2 figure to '{out_pdf}'")


# ------------------------------------------------------------------------------
# Summaries: Occurrence Bars + Existence Heatmap
# ------------------------------------------------------------------------------

def count_oxygen_occurrences_from_matrix(hbond_df, hbonds_per_index):
    """
    Original function that sums total frames from both donor & acceptor sides.
    Not used for "Time / Atom," but used for the final summary existence plot.
    """
    logger.info("Counting oxygen occurrences (both donor & acceptor).")
    df = hbond_df.copy()
    if 'count' in df.columns:
        df = df.drop(columns=['count'])
    df['count'] = hbonds_per_index[df['idx'].values]

    melted = df.melt(
        id_vars=['idx','count'],
        value_vars=['donor','acceptor'],
        value_name='oxygen'
    ).dropna(subset=['oxygen'])
    counts = melted.groupby('oxygen')['count'].sum()
    logger.info(f"Found {len(counts)} total oxygen atoms.")
    return counts


def generate_summary_plot(
    occurrence_data,
    existence_data,
    donor_atoms_by_species=None,
    color_dict=None,
    time_per_frame_ns=0.01,
    output_file='oxygen_occurrences_summary.pdf',
    folder_order=None
):
    """
    Snippet-like summary plot:
      Left: barplot of occurrence times (reversed x-axis),
      Right: existence heatmap.

    If donor_atoms_by_species is given, any oxygen in that set is labeled red on y-axis.
    """
    logger.info("Generating summary plot with existence maps & occurrence bars.")
    num_species = len(occurrence_data)
    if num_species == 0:
        logger.warning("No data in occurrence_data; skipping summary plot.")
        return

    if donor_atoms_by_species is None:
        donor_atoms_by_species = {}

    # Gather all possible oxygens
    all_oxygens = set()
    for ser in occurrence_data.values():
        all_oxygens.update(ser.index)
    all_oxygens = sorted(
        all_oxygens,
        key=lambda x: int(re.findall(r'\d+', x)[0]) if re.findall(r'\d+', x) else 0
    )

    if folder_order:
        sorted_species_ids = [sp for sp in folder_order if sp in occurrence_data]
    else:
        sorted_species_ids = sorted(occurrence_data.keys())

    max_occ = max([s.max() for s in occurrence_data.values()]) if occurrence_data else 0
    x_limit = max_occ * time_per_frame_ns * 1.1

    fig_height = max(6, 5.5*num_species)
    fig, axes = plt.subplots(
        nrows=num_species, ncols=2,
        figsize=(16, fig_height),
        constrained_layout=True,
        gridspec_kw={'width_ratios':[1,1]}
    )
    if num_species == 1:
        axes = [axes]

    for i, sp in enumerate(sorted_species_ids):
        ax_occ, ax_exist = axes[i]
        occ_series = occurrence_data[sp].reindex(all_oxygens, fill_value=0)
        time_vals = occ_series * time_per_frame_ns

        bar_color = color_dict.get(sp, "gray") if color_dict else "gray"

        # Left subplot: bars
        sns.barplot(
            x=time_vals.values,
            y=all_oxygens,
            color=bar_color,
            ax=ax_occ
        )
        ax_occ.set_title(f"H-Bond Occurrence Time - {sp}", fontweight='bold')
        ax_occ.set_xlabel("Time (ns)")
        ax_occ.set_xlim(x_limit, 0)

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
        ax_occ.set_yticklabels([])
        ax_occ.yaxis.set_ticks_position('right')

        # Right subplot: existence heatmap
        mat = existence_data[sp]
        if mat.shape[0] != len(all_oxygens):
            logger.warning(f"Mismatch for {sp}: existence map {mat.shape[0]} rows vs {len(all_oxygens)} oxygens.")
            ax_exist.set_visible(False)
            continue

        sns.heatmap(
            mat, cmap="Reds", cbar=True, ax=ax_exist,
            linewidths=0, linecolor='white',
            vmin=0, vmax=1,
            cbar_kws={"label":"Fraction of Bond Presence"}
        )
        ax_exist.set_title(f"H-Bond Existence Map - {sp}", fontweight='bold')
        ax_exist.set_xlabel("Time (ns)")

        nbins = mat.shape[1]
        bin_size_ns = 0.2
        tvals = np.arange(nbins)*bin_size_ns
        ticks_count = min(6, nbins)
        tick_positions = np.linspace(0, nbins-1, ticks_count, dtype=int)
        tick_labels = [f"{int(tvals[pos])}" for pos in tick_positions]

        ax_exist.set_xticks(tick_positions + 0.5)
        ax_exist.set_xticklabels(tick_labels, rotation=0, ha='right')
        ax_exist.set_yticks(np.arange(len(all_oxygens)) + 0.5)
        ax_exist.set_yticklabels(all_oxygens, rotation=0)

        # Color donor y-labels red
        donorset = donor_atoms_by_species.get(sp, set())
        for txt_label in ax_exist.get_yticklabels():
            oxy_txt = txt_label.get_text()
            if oxy_txt in donorset:
                txt_label.set_color("red")

        for spine in ax_exist.spines.values():
            spine.set_visible(True)

    plt.savefig(output_file, dpi=300)
    plt.close()
    logger.info(f"Summary plot saved as '{output_file}'")


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

def main():
    """
    Main routine orchestrating:
      1) Scan 'process.nobackup' for species with 3/4/5 '1'-bits
      2) Load H-bond data => produce 2×2 figure:
         - (Lifetime box, Distance box, Time/H-Atom bar [donors only], H-bond count violin)
      3) Also produce a summary snippet plot (occurrence bars + existence heatmap),
         coloring donor oxygens red, same color dictionary.
    """
    logger.info("Starting combined analysis with donor-based engaged time, 2×2 layout...")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    try:
        project_path = find_project_root(script_dir)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)

    process_dir = os.path.join(project_path, "process.nobackup")
    if not os.path.isdir(process_dir):
        logger.error(f"Invalid process directory: {process_dir}")
        sys.exit(1)
    logger.info(f"Using process directory: {process_dir}")

    process_path = Path(process_dir)
    time_per_frame_ns = 0.01

    # We gather data for 2x2 figure + a final snippet plot
    hbonds_data    = []
    lifetimes_data = []
    distance_data  = []

    # For computing Time/H-Atom, we sum only donor occurrences
    occurrence_data_for_time_atom = {}  # { species_id -> sum of donor frames }
    # We'll also store donor sets per species => used in final summary to color them red
    donor_atoms_by_species = {}
    # We'll keep track of species IDs as well
    all_species_ids = []

    # (1) Collect subfolders that start with "IP_"
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

        # required files
        hb_num_file = folder_path / "analyze_final_sim" / "h_bonds" / "intra_hb_num.xvg"
        xpm_file    = folder_path / "analyze_final_sim" / "h_bonds" / "intra_hb_matrix.xpm"
        dist_file   = folder_path / "analyze_final_sim" / "h_bonds" / "intra_hb_dist.xvg"
        hb_log_file = folder_path / "analyze_final_sim" / "h_bonds" / "intra_hb.log"

        if not (hb_num_file.is_file() and xpm_file.is_file() and dist_file.is_file() and hb_log_file.is_file()):
            logger.warning(f"Skipping {species_id}: missing input files.")
            continue

        # -- A) H-bond Count vs Time --
        data_num = load_hb_num_xvg(str(hb_num_file))
        if data_num.size < 2:
            logger.warning(f"No hbond count data for {species_id}")
            continue
        hbond_counts = data_num[:,1]
        for val in hbond_counts:
            hbonds_data.append((species_id, num_ones, val))

        # -- B) XPM => matrix => analyze => lifetimes
        xpm_mat, meta = parse_xpm(str(xpm_file))
        analysis_res  = analyze_hydrogen_bonds(xpm_mat, meta)

        # lifetimes => frames => ns
        lifetimes_frames = [lf for bond_lf in analysis_res['lifetimes'] for lf in bond_lf]
        lifetimes_ns = [f * time_per_frame_ns for f in lifetimes_frames]
        for lf_ns in lifetimes_ns:
            lifetimes_data.append((species_id, num_ones, lf_ns))

        # -- C) Distance
        d_list = []
        with open(dist_file, 'r') as dfile:
            for line in dfile:
                if not line.strip() or line.startswith(('#','@')):
                    continue
                parts = line.split()
                if len(parts)==2:
                    try:
                        dx = float(parts[0])
                        dy = float(parts[1])
                        d_list.append((dx, dy))
                    except ValueError:
                        pass
        arr_dist = np.array(d_list)
        if arr_dist.size<2:
            logger.warning(f"No distance data for {species_id}")
            continue
        arr_dist = arr_dist[arr_dist[:,0]>=0.2]
        for row in arr_dist:
            distance_data.append((species_id, num_ones, row[0], row[1]))

        # -- D) Donor-based sum for engaged time
        # parse .log => DataFrame => sum frames only for donors
        hbond_df = parse_hbond_log_to_dataframe(str(hb_log_file))
        if hbond_df.empty:
            logger.warning(f"No hbond pairs for {species_id}")
            continue
        # track which oxygens are donors => used for color in summary
        local_donors = set(hbond_df['donor'].unique())
        donor_atoms_by_species[species_id] = local_donors

        hpi = analysis_res['hbonds_per_index']
        donor_counts = count_donor_occurrences_from_matrix(hbond_df, hpi)
        occurrence_data_for_time_atom[species_id] = donor_counts

    # (2) Build DataFrames for the 2x2 figure
    df_hbonds = pd.DataFrame(hbonds_data, columns=["species","num_ones","hbond_count"])
    df_lifetimes = pd.DataFrame(lifetimes_data, columns=["species","num_ones","lifetime_ns"])
    df_distance = pd.DataFrame(distance_data, columns=["species","num_ones","distance_nm","frequency"])

    logger.info(f"Hbond DF: {df_hbonds.shape}, Lifetime DF: {df_lifetimes.shape}, Dist DF: {df_distance.shape}")

    # (3) Build DataFrame => Time/H-Atom from the donor-based occurrences
    engaged_list = []
    for sp, donor_occ_series in occurrence_data_for_time_atom.items():
        # sum frames across all donor oxygens
        total_frames = donor_occ_series.sum()
        # convert frames => ns
        total_time_ns = total_frames * time_per_frame_ns
        # # of donors = # of '1' bits in species ID
        n_ones = sp.count('1')
        if n_ones<=0:
            continue
        engaged_time_ns = total_time_ns / n_ones
        engaged_list.append((sp, n_ones, engaged_time_ns))

    df_engaged = pd.DataFrame(engaged_list, columns=["species","num_ones","engaged_time_ns"])
    logger.info(f"Engaged DF: {df_engaged.shape}")

    # (4) Color dict (tab10)
    all_species_ids = sorted(set(all_species_ids))
    palette = sns.color_palette("tab10", len(all_species_ids))
    color_dict = {}
    for i, sp in enumerate(all_species_ids):
        color_dict[sp] = palette[i % len(palette)]

    # (5) 2×2 figure
    two_by_two_plots_box_violin(
        df_lifetimes=df_lifetimes,
        df_distance=df_distance,
        df_engaged=df_engaged,      # bottom-left
        df_hbonds=df_hbonds,        # bottom-right
        color_dict=color_dict,
        output_prefix="intra_hbonds",
        figure_size=(16,12)
    )

    # (6) Final snippet summary (occurrence + existence map), 
    #     but with full donors/acceptors (not just donors).
    occurrence_data = {}
    existence_data = {}
    folder_list = []

    for folder_name in sorted(os.listdir(process_path)):
        if not folder_name.startswith("IP_"):
            continue
        species_id = folder_name.replace("IP_","")
        num_ones = species_id.count('1')
        if num_ones not in [3,4,5]:
            continue
        folder_list.append((folder_name, species_id, num_ones))

    # sort 5->4->3
    folder_list.sort(key=lambda x: (-x[2], x[1]))
    sorted_ids = [x[1] for x in folder_list]
    all_oxygens = set()

    # first pass => collect all oxygens
    for fname, sp, _ in folder_list:
        fpath = process_path / fname / "analyze_final_sim" / "h_bonds" / "intra_hb.log"
        if not fpath.is_file():
            continue
        bdf = parse_hbond_log_to_dataframe(str(fpath))
        if bdf.empty:
            continue
        all_oxygens.update(bdf['donor'].unique())
        all_oxygens.update(bdf['acceptor'].unique())

    # second pass => existence data
    for fname, sp, _ in folder_list:
        folder_path = process_path / fname
        xpm_file = folder_path / "analyze_final_sim" / "h_bonds" / "intra_hb_matrix.xpm"
        hb_log_file = folder_path / "analyze_final_sim" / "h_bonds" / "intra_hb.log"
        if not (xpm_file.is_file() and hb_log_file.is_file()):
            continue

        bond_df = parse_hbond_log_to_dataframe(str(hb_log_file))
        if bond_df.empty:
            continue

        mat, meta = parse_xpm(str(xpm_file))
        an = analyze_hydrogen_bonds(mat, meta)
        hpi = an['hbonds_per_index']

        # store occurrence data (donor+acceptor) for the snippet
        c_ser = count_oxygen_occurrences_from_matrix(bond_df, hpi)
        occurrence_data[sp] = c_ser

        # existence map
        bin_size_ns = 0.2
        frames_per_bin = int(bin_size_ns / time_per_frame_ns)
        if frames_per_bin<1:
            continue
        n_frames = mat.shape[1]
        nbins = n_frames // frames_per_bin
        if n_frames % frames_per_bin != 0:
            nbins +=1

        binned = np.zeros((mat.shape[0], nbins))
        for i in range(nbins):
            st = i*frames_per_bin
            en = min((i+1)*frames_per_bin, n_frames)
            chunk = mat[:, st:en]
            binned[:,i] = np.mean(chunk, axis=1)

        # bond->oxygen
        oxy_map = {}
        for idx, row in bond_df.iterrows():
            d = row['donor']
            a = row['acceptor']
            oxy_map.setdefault(d, []).append(idx)
            oxy_map.setdefault(a, []).append(idx)

        # build aggregated existence
        oxy_sorted = sorted(all_oxygens, key=lambda x: int(re.findall(r'\d+', x)[0]) if re.findall(r'\d+', x) else 0)
        agg = np.zeros((len(oxy_sorted), nbins))
        for o_i, oxy in enumerate(oxy_sorted):
            bond_indices = oxy_map.get(oxy, [])
            if bond_indices:
                for bn in range(nbins):
                    agg[o_i, bn] = np.max(binned[bond_indices, bn])
        existence_data[sp] = agg

    # (7) Final snippet-based summary plot, coloring donor atoms red
    generate_summary_plot(
        occurrence_data=occurrence_data,
        existence_data=existence_data,
        donor_atoms_by_species=donor_atoms_by_species,
        color_dict=color_dict,
        time_per_frame_ns=time_per_frame_ns,
        output_file='oxygen_occurrences_summary.pdf',
        folder_order=sorted_ids
    )

    logger.info("All analysis completed successfully!")


if __name__ == "__main__":
    main()