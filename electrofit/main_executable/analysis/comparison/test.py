#!/usr/bin/env python3
"""
Combined Script for H-Bond Analysis with:
  1) A 2×2 Figure (Lifetime, Distance, Time/Atom, H-Bond Count)
  2) A Summary Plot (Occurrence bars + Existence Heatmap)
  3) A NEW Phosphate Group Violin Plot:
     - Aggregating donor occurrences by phosphate label (P1..P6)
     - Each species (micro protonation state) shown as separate distribution

We focus on IP6 microprotonation states. 
If a donor is in [O7,O8,O9] => P1, [O10,O11,O12] => P2, ..., [O22,O23,O24] => P6.

We force Times New Roman font, color donor oxygens red in the summary plot, 
and produce a new file 'phosphate_group_violin.pdf'.

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

# ------------------------------------------------------------------------------
# Global Logging & Font Setup
# ------------------------------------------------------------------------------
logging.basicConfig(
    filename='analysis.log',
    filemode='w',  # overwrite each run each time
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)

# Use Times New Roman globally
mpl.rcParams['font.family'] = 'Times New Roman'
sns.set_context("talk")


# ------------------------------------------------------------------------------
# Donor & Acceptor → Phosphate Mapping
# ------------------------------------------------------------------------------
donor_to_phosphate = {
    # P1
    'O1': 'P1','O7': 'P1','O8': 'P1','O9': 'P1',
    # P2
    'O2': 'P2','O10': 'P2','O11': 'P2','O12': 'P2',
    # P3
    'O3': 'P3','O13': 'P3','O14': 'P3','O15': 'P3',
    # P4
    'O4': 'P4','O16': 'P4','O17': 'P4','O18': 'P4',
    # P5
    'O5': 'P5','O19': 'P5','O20': 'P5','O21': 'P5',
    # P6
    'O6': 'P6','O22': 'P6','O23': 'P6','O24': 'P6'
}

# If acceptors also appear in [O1..O24], we can reuse the same mapping:
acceptor_to_phosphate = donor_to_phosphate




# ------------------------------------------------------------------------------
# Common Utilities
# ------------------------------------------------------------------------------
def refine_atom_name(atom):
    """
    In GROMACS logs:
      - 'O' => 'O1'
      - 'O10' => 'O11'
    (Increment the trailing digit by 1 if present.)
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
    """Ascend from current_dir until we find a folder named project_name."""
    logger.info(f"Locating project root '{project_name}' from: {current_dir}")
    root = None
    while True:
        parent = os.path.dirname(current_dir)
        if os.path.basename(current_dir) == project_name:
            root = current_dir
        if parent == current_dir:  # Reached filesystem root
            if root is None:
                raise FileNotFoundError(f"Project root '{project_name}' not found.")
            logger.info(f"Found project root: {root}")
            return root
        current_dir = parent


def load_hb_num_xvg(filename):
    """
    Loads (#H-bonds vs. time) from .xvg (skips lines with @/#).
    Returns Nx2 array => [time, hbond_count].
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
    logger.info(f"XVG shape: {arr.shape}")
    return arr


def parse_xpm(file_path):
    """
    Parse XPM => binary array (#FF0000 => 1). Return (data_matrix, metadata).
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

    for idx, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()

        if line.startswith("/*") and not data_started:
            comment = line.strip("/* ").strip(" */")
            if ':' in comment:
                key, value = comment.split(":", 1)
                metadata[key.strip().lower()] = value.strip().strip('"')
            continue

        if line.startswith('static char'):
            continue

        if (not header_found) and line.startswith('"'):
            # something like: "64 64 10 1",
            line_str = line.strip(',').lstrip('"').rstrip('"')
            tokens = line_str.split()
            if len(tokens) >= 4:
                try:
                    width, height, num_colors, chars_per_pixel = map(int, tokens[:4])
                    header_found = True
                    logger.info(f"XPM header => w={width}, h={height}, c={num_colors}, cpp={chars_per_pixel}")
                    continue
                except ValueError:
                    logger.error(f"Invalid header line {idx}: {raw_line}")
                    raise
            else:
                logger.error(f"Header line {idx} missing tokens: {raw_line}")
                raise ValueError(f"Bad XPM header line: {raw_line}")

        if header_found and not data_started:
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
            row_data = line.strip(',').lstrip('"').rstrip('"')
            data_lines.append(row_data)

    if width is None or height is None:
        logger.error("No valid XPM header found.")
        raise ValueError("XPM header missing.")

    if len(data_lines) != height:
        logger.warning(f"XPM mismatch: expected {height} lines, found {len(data_lines)}")

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
    Summarize a binary matrix => over_time, per_index, lifetimes.
    """
    logger.info("Analyzing hydrogen bond matrix.")
    hbonds_over_time = np.sum(data_matrix, axis=0)
    hbonds_per_index = np.sum(data_matrix, axis=1)

    lifetimes = []
    for row in data_matrix:
        current = 0
        bond_lf = []
        for val in row:
            if val == 1:
                current += 1
            else:
                if current>0:
                    bond_lf.append(current)
                    current = 0
        if current>0:
            bond_lf.append(current)
        lifetimes.append(bond_lf)

    return {
        'hbonds_over_time': hbonds_over_time,
        'hbonds_per_index': hbonds_per_index,
        'lifetimes': lifetimes
    }


def parse_hbond_log_to_dataframe(file_path):
    """
    .log => DataFrame [idx, donor, acceptor], refining final numeric bits.
    """
    logger.info(f"Parsing H-bond log: {file_path}")
    pairs = []
    line_pattern = re.compile(r'^\s*(\S+)\s+-\s+(\S+)\s*$')
    atom_pattern = re.compile(r'^[A-Za-z]+\d+([A-Za-z]+\d*)$')

    with open(file_path, 'r') as f:
        for line_no, raw_line in enumerate(f, start=1):
            line = raw_line.strip()
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
                    logger.warning(f"Line {line_no}: parse fail => {line}")
            else:
                logger.warning(f"Line {line_no}: no match => {line}")

    df = pd.DataFrame(pairs)
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'idx'}, inplace=True)
    df['idx'] = df.index
    logger.info(f"Found {len(df)} H-bond pairs.")
    return df





# ------------------------------------------------------------------------------
# Donor-based Summation (Donor Only)
# ------------------------------------------------------------------------------
def count_donor_occurrences_from_matrix(hbond_df, hbonds_per_index):
    """
    Only sums frames for the 'donor' column in hbond_df. 
    Returns a Series: donor_atom -> total frames.
    """
    logger.info("Counting DONOR-based frames only (avoid double-counting).")
    df = hbond_df.copy()
    if 'count' in df.columns:
        df.drop(columns=['count'], inplace=True)
    df['count'] = hbonds_per_index[df['idx'].values]

    melted = df.melt(
        id_vars=['idx','count'],
        value_vars=['donor'],  # ignoring acceptor
        value_name='oxygen'
    ).dropna(subset=['oxygen'])

    counts = melted.groupby('oxygen')['count'].sum()
    logger.info(f"Donor-based occurrences => {len(counts)} donors.")
    return counts

# ------------------------------------------------------------------------------
# Donor-Acceptor-based Summation
# ------------------------------------------------------------------------------

def count_oxygen_occurrences_from_matrix(hbond_df, hbonds_per_index):
    """
    Summation for both donor & acceptor => used for final existence heatmap.
    """
    logger.info("Counting oxygen occurrences (both donor & acceptor).")
    df=hbond_df.copy()
    if 'count' in df.columns:
        df.drop(columns=['count'],inplace=True)
    df['count']=hbonds_per_index[df['idx'].values]

    melted=df.melt(
        id_vars=['idx','count'],
        value_vars=['donor','acceptor'],
        value_name='oxygen'
    ).dropna(subset=['oxygen'])
    grouped=melted.groupby('oxygen')['count'].sum()
    logger.info(f"Found {len(grouped)} total oxygen atoms.")
    return grouped

# ------------------------------------------------------------------------------
# Summaries: Occurrence Bars + Existence Heatmap
# ------------------------------------------------------------------------------


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
# Build/Print Table: Donor Phosphate => top acceptor targets
# ------------------------------------------------------------------------------

def build_donor_acceptor_summary(
    species_id,
    bond_df,            
    hbonds_per_index,   
    time_per_frame_ns,
    donor_map,
    acceptor_map
):
    """
    For each bond in bond_df => (donor_atom, acceptor_atom).
    - Convert donor_atom -> donor_pgroup
    - Convert acceptor_atom -> acceptor_pgroup
    - engaged_ns = frames * time_per_frame_ns
    Store in a nested dict => data[donor_pgroup]['sum_donor_time'] and data[donor_pgroup]['targets'][(acc_pgroup, acc_atom)].
    Returns that nested structure for printing.
    """
    result = {}

    for i, row in bond_df.iterrows():
        b_idx = row['idx']
        d_atom= row['donor']
        a_atom= row['acceptor']

        frames = hbonds_per_index[b_idx]
        engaged_ns= frames*time_per_frame_ns

        donor_pgroup = donor_map.get(d_atom,None)
        acceptor_pgroup = acceptor_map.get(a_atom,None)
        if donor_pgroup is None or acceptor_pgroup is None:
            # skip if unknown
            continue

        if donor_pgroup not in result:
            result[donor_pgroup] = {
                "sum_donor_time": 0.0,
                "targets": {}
            }
        # accumulate total
        result[donor_pgroup]["sum_donor_time"] += engaged_ns

        # accumulate target details => (acceptor_pgroup, a_atom)
        key_t=(acceptor_pgroup, a_atom)
        if key_t not in result[donor_pgroup]["targets"]:
            result[donor_pgroup]["targets"][key_t]=0.0
        result[donor_pgroup]["targets"][key_t]+= engaged_ns

    return result


def print_donor_acceptor_table(species_id, data):
    """
    data => nested dict: data[donor_pgroup]['sum_donor_time'] + .targets => (acc_pg,acc_atom)->ns
    Print lines:

      Species: 010101

        Donor P2 => total time: XX ns
          main target Pi => P3(O3) = 6.43 ns
          second target Pi => ...
          third target Pi => ...
    """
    print(f"Species: {species_id}\n")

    # sort donor pgroups => P1..P6
    donor_sorted = sorted(data.keys())
    for dpgroup in donor_sorted:
        sum_time = data[dpgroup]["sum_donor_time"]
        print(f"  Donor {dpgroup} => total time: {sum_time:.2f} ns")

        # top 3 targets
        t_map = data[dpgroup]["targets"]
        t_list = sorted(t_map.items(), key=lambda x: x[1], reverse=True)  # sorted by time desc
        top_names=["main target Pi","second target Pi","third target Pi"]
        for i,( (acc_pg,acc_atom), val_ns ) in enumerate(t_list):
            if i<3:
                print(f"    {top_names[i]:20s} => {acc_pg}({acc_atom}) = {val_ns:.2f} ns")
        print()
    print("-"*70)
    print()

# ------------------------------------------------------------------------------
# Two-by-Two Figure Layout
# ------------------------------------------------------------------------------
def two_by_two_plots_box_violin(
    df_lifetimes,   # [species, num_ones, lifetime_ns]
    df_distance,    # [species, num_ones, distance_nm, frequency]
    df_engaged,     # [species, num_ones, engaged_time_ns]
    df_hbonds,      # [species, num_ones, hbond_count]
    color_dict=None,
    output_prefix="intra_hbonds",
    figure_size=(16, 12),
    spacer_prefix="~space_"
):
    """
    2×2 layout:
      (0,0) => Lifetime (Box)
      (0,1) => Distance (Box)
      (1,0) => Time/Atom (Bar)
      (1,1) => H-bond Count (Violin)

    Insert spacers for species grouping (5->4->3).
    """
    logger.info("Generating 2×2 figure with swapped bottom row (Time/Atom, HbondCount).")

    # Quick check if all are empty
    all_empty = True
    for df in [df_lifetimes, df_distance, df_engaged, df_hbonds]:
        if df is not None and not df.empty:
            all_empty = False
            break
    if all_empty:
        logger.warning("All DataFrames empty => skipping 2×2 plot.")
        return

    # Gather species
    species_in_any = set()
    for df in [df_lifetimes, df_distance, df_engaged, df_hbonds]:
        if df is not None and not df.empty:
            species_in_any.update(df["species"].unique())

    def get_num_ones(sp, df_list):
        for d in df_list:
            if d is not None and not d.empty:
                row = d.loc[d["species"]==sp]
                if not row.empty:
                    return row.iloc[0]["num_ones"]
        return None

    # Only keep those with num_ones in [3,4,5]
    all_with_num = []
    for sp in species_in_any:
        val = get_num_ones(sp, [df_lifetimes, df_distance, df_engaged, df_hbonds])
        if val in [3,4,5]:
            all_with_num.append((sp,val))
    if not all_with_num:
        logger.warning("No species with 3,4,5 ones => skipping 2×2.")
        return

    # Sort: 5->top, then 4->3
    all_with_num.sort(key=lambda x: (-x[1], x[0]))
    # Insert spacers
    grouped_species_order = []
    prev_n_ones = None
    for i, (sp,n) in enumerate(all_with_num):
        if i>0 and n!=prev_n_ones:
            grouped_species_order.append((f"{spacer_prefix}{prev_n_ones}to{n}", None))
        grouped_species_order.append((sp,n))
        prev_n_ones = n
    full_order = [t[0] for t in grouped_species_order]

    def inject_spacers(df, measure_col):
        if df is None or df.empty:
            return pd.DataFrame(columns=["species", measure_col, "num_ones"])
        new_df = df.copy()
        row_spacers = []
        for sp, val in grouped_species_order:
            if val is None:  # indicates dummy label
                row_spacers.append({
                    "species": sp,
                    measure_col: np.nan,
                    "num_ones": -1
                })
        if row_spacers:
            dummy_df = pd.DataFrame(row_spacers)
            return pd.concat([new_df, dummy_df], ignore_index=True)
        return new_df

    # Expand distance if needed
    if df_distance is not None and not df_distance.empty:
        df_distance = df_distance.copy()
        df_distance['frequency'] = df_distance['frequency'].astype(int)
        df_dist_expanded = df_distance.loc[
            df_distance.index.repeat(df_distance['frequency'])
        ].reset_index(drop=True)
    else:
        df_dist_expanded = pd.DataFrame()

    df_life2 = inject_spacers(df_lifetimes, "lifetime_ns")
    df_dist2 = inject_spacers(df_dist_expanded, "distance_nm")
    df_eng2  = inject_spacers(df_engaged, "engaged_time_ns")
    df_hbond2= inject_spacers(df_hbonds, "hbond_count")

    # Ensure color_dict has spacer color
    if color_dict is None:
        color_dict = {}
    else:
        color_dict = color_dict.copy()
    spacer_color = "#D3D3D3"
    for sp,val in grouped_species_order:
        if val is None and sp not in color_dict:
            color_dict[sp] = spacer_color

    fig, axes = plt.subplots(2,2, figsize=figure_size)
    ax_life = axes[0][0]
    ax_dist = axes[0][1]
    ax_eng  = axes[1][0]  # Time/Atom
    ax_hbond= axes[1][1]  # H-bond count

    # (1) Lifetime
    if not df_life2.empty:
        sns.boxplot(
            data=df_life2,
            x="lifetime_ns", y="species",
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
        ax_life.set_ylabel("Species")
        # Hide spacer labels
        new_lbls=[]
        for lbl in ax_life.get_yticklabels():
            txt=lbl.get_text()
            if txt.startswith(spacer_prefix):
                new_lbls.append("")
            else:
                new_lbls.append(txt)
        ax_life.set_yticklabels(new_lbls)
    else:
        ax_life.set_title("No Lifetime Data")
        ax_life.set_xlabel("")
        ax_life.set_ylabel("")

    # (2) Distance
    if not df_dist2.empty:
        sns.boxplot(
            data=df_dist2,
            x="distance_nm", y="species",
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

    # (3) Time/Atom (bottom-left)
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
        ax_eng.set_ylabel("")
        ax_eng.set_yticklabels([])

        # Larger annotation font
        for patch in ax_eng.patches:
            width = patch.get_width()
            y_center = patch.get_y() + patch.get_height()/2
            if not np.isnan(width) and width>0:
                ax_eng.annotate(
                    f"{width:.2f}",
                    (width, y_center),
                    ha='right',
                    va='center',
                    xytext=(-5,0),
                    textcoords='offset points',
                    fontsize=14  # bigger font
                )
    else:
        ax_eng.set_title("No Engaged-Time Data")
        ax_eng.set_xlabel("")
        ax_eng.set_ylabel("")
        ax_eng.set_yticklabels([])

    # (4) H-bond count (bottom-right)
    if not df_hbond2.empty:
        sns.violinplot(
            data=df_hbond2,
            x="hbond_count", y="species",
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
            color='black', marker='o', s=80, label='Mean', zorder=5
        )
        ax_hbond.legend(loc='lower right')
    else:
        ax_hbond.set_title("No H-bond Data")
        ax_hbond.set_xlabel("")
        ax_hbond.set_ylabel("")
        ax_hbond.set_yticklabels([])

    plt.tight_layout()
    out_pdf = f"{output_prefix}_2x2_with_engaged.pdf"
    plt.savefig(out_pdf, dpi=300)
    plt.close()
    logger.info(f"Saved 2×2 figure to '{out_pdf}'")


# ------------------------------------------------------------------------------
# New Plot: Donor → Phosphate Group Distribution (Violin)
# ------------------------------------------------------------------------------
def plot_phosphate_group_violin(
    df_phosphate,
    color_dict=None,
    output_file="phosphate_group_violin.pdf"
):
    """
    Creates a violin plot showing:
      x-axis = phosphate group label (P1..P6),
      y-axis = engaged_ns (the total time in ns for each donor),
      hue    = species (microprotonation state).

    This compares how donors in each phosphate group 
    are distributed across the 11 species.
    """
    logger.info("Generating phosphate group violin plot...")

    if df_phosphate is None or df_phosphate.empty:
        logger.warning("No phosphate-group donor data => skipping violin plot.")
        return
    
    plt.figure(figsize=(10,6))

    # We'll do x='phosphate_group', y='engaged_ns', hue='species'
    sns.violinplot(
        data=df_phosphate,
        x='phosphate_group',
        y='engaged_ns',
        hue='species',
        palette=color_dict if color_dict else "tab10",
        cut=0
    )
    plt.title("Donor Engagement by Phosphate Group (Violin)")
    plt.xlabel("Phosphate Group")
    plt.ylabel("Engaged Time (ns)")

    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    logger.info(f"Saved phosphate group violin plot => '{output_file}'")


# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
def main():
    """
    Main routine that:
      1) Finds 'process.nobackup'
      2) Gathers species with 3..5 ones
      3) Creates:
         (a) A 2×2 figure: Lifetime, Distance, Time/Atom, H-bond Count
         (b) A summary plot (occurrence bars + existence heatmap)
         (c) A new phosphate-group violin plot using donor_to_phosphate
    """
    logger.info("Starting combined analysis with new phosphate group violin plot...")

    script_dir = os.path.dirname(os.path.abspath(__file__))
    try:
        project_path = find_project_root(script_dir)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)

    process_dir = os.path.join(project_path, "process.nobackup")
    if not os.path.isdir(process_dir):
        logger.error(f"Invalid process dir => {process_dir}")
        sys.exit(1)
    logger.info(f"Using process directory => {process_dir}")

    time_per_frame_ns = 0.01
    process_path = Path(process_dir)

    # We'll collect data for:
    #   - (2x2) => hbonds_data, lifetimes_data, distance_data, engaged_data
    #   - (summary) => occurrence_data, existence_data
    #   - (phosphate group) => df_phosphate (rows = donors, columns = [species, phosphate_group, engaged_ns])
    hbonds_data     = []
    lifetimes_data  = []
    distance_data   = []
    engaged_data    = []  # for the 2x2 time/atom bar

    # For the new phosphate-group distribution
    df_phosphate_rows = []  # Will store [species, phosphate_group, engaged_ns]

    # We'll also hold onto donor sets for each species => color them red in summary
    donor_atoms_by_species = {}
    # For final summary existence data
    occurrence_data = {}
    existence_data  = {}
    donor_acceptor_summaries = {}

    # We'll gather species from subfolders
    all_species_ids = []

    folder_list = []
    for fname in sorted(os.listdir(process_path)):
        if not fname.startswith("IP_"):
            continue
        folder_path = process_path / fname
        if not folder_path.is_dir():
            continue

        species_id = fname.replace("IP_","")
        n_ones = species_id.count('1')
        if n_ones not in [3,4,5]:
            continue
        folder_list.append((fname,species_id,n_ones))
        all_species_ids.append(species_id)

    # Sort 5->4->3
    folder_list.sort(key=lambda x: (-x[2], x[1]))
    sorted_species_list = [x[1] for x in folder_list]

    for folder_name, species_id, n_ones in folder_list:
        # Check required files
        fpath = process_path / folder_name
        hb_num_file = fpath / "analyze_final_sim" / "h_bonds" / "intra_hb_num.xvg"
        xpm_file    = fpath / "analyze_final_sim" / "h_bonds" / "intra_hb_matrix.xpm"
        dist_file   = fpath / "analyze_final_sim" / "h_bonds" / "intra_hb_dist.xvg"
        log_file    = fpath / "analyze_final_sim" / "h_bonds" / "intra_hb.log"

        if not (hb_num_file.is_file() and xpm_file.is_file() and dist_file.is_file() and log_file.is_file()):
            logger.warning(f"Skipping {species_id}: missing required hbond files.")
            continue

        # 1) H-bond counts vs time
        data_num = load_hb_num_xvg(str(hb_num_file))
        if data_num.size<2:
            logger.warning(f"No hbond count data for {species_id}")
            continue
        hbond_counts = data_num[:,1]
        for val in hbond_counts:
            hbonds_data.append((species_id,n_ones,val))

        # 2) Parse XPM => analyze => lifetimes
        xpm_mat, meta = parse_xpm(str(xpm_file))
        analysis_res = analyze_hydrogen_bonds(xpm_mat, meta)
        # lifetimes in frames => convert => ns
        lf_frames = [lf for bond_lf in analysis_res['lifetimes'] for lf in bond_lf]
        lf_ns = [f * time_per_frame_ns for f in lf_frames]
        for val in lf_ns:
            lifetimes_data.append((species_id,n_ones,val))

        # 3) Distances
        d_list=[]
        with open(dist_file, 'r') as df_:
            for line in df_:
                if not line.strip() or line.startswith(('#','@')):
                    continue
                parts=line.split()
                if len(parts)==2:
                    try:
                        dx=float(parts[0])
                        dy=float(parts[1])
                        d_list.append((dx,dy))
                    except ValueError:
                        pass
        arr_d = np.array(d_list)
        if arr_d.size<2:
            logger.warning(f"No distance data => {species_id}")
            continue
        arr_d = arr_d[arr_d[:,0]>=0.2]
        for row in arr_d:
            distance_data.append((species_id,n_ones,row[0], row[1]))

        # 4) Parse .log => donor-based sum => time/atom
        hbond_df = parse_hbond_log_to_dataframe(str(log_file))
        if hbond_df.empty:
            logger.warning(f"No hbond pairs => {species_id}")
            continue
        donor_atoms_by_species[species_id] = set(hbond_df['donor'].unique())

        # frames per bond index
        hpi = analysis_res['hbonds_per_index']
        donor_counts = count_donor_occurrences_from_matrix(hbond_df, hpi)
        total_donor_frames = donor_counts.sum()
        engaged_ns = (total_donor_frames * time_per_frame_ns) / n_ones
        engaged_data.append((species_id, n_ones, engaged_ns))

        # 4b) Build the new phosphate-group data
        # For each donor in donor_counts => find phosphate group => store row
        for donor_atom, frames in donor_counts.items():
            # Convert frames => ns
            engaged_time_ns = frames * time_per_frame_ns
            # Map donor => phosphate label (if found)
            phosphate_label = donor_to_phosphate.get(donor_atom, None)
            if phosphate_label is not None:
                df_phosphate_rows.append({
                    'species': species_id,
                    'num_ones': n_ones,
                    'donor_atom': donor_atom,
                    'phosphate_group': phosphate_label,
                    'engaged_ns': engaged_time_ns
                })
            else:
                # Possibly a donor not in O7..O24 => skip
                pass

        summary_data = build_donor_acceptor_summary(species_id, bond_df, hpi, time_per_frame_ns, donor_to_phosphate, acceptor_to_phosphate)
        donor_acceptor_summaries[species_id] = summary_data

    # Convert into DataFrames for the 2x2 figure
    df_hbonds = pd.DataFrame(hbonds_data, columns=["species","num_ones","hbond_count"])
    df_lifetimes = pd.DataFrame(lifetimes_data, columns=["species","num_ones","lifetime_ns"])
    df_distance = pd.DataFrame(distance_data, columns=["species","num_ones","distance_nm","frequency"])

    df_engaged = pd.DataFrame(engaged_data, columns=["species","num_ones","engaged_time_ns"])
    logger.info(f"Prepared data => hbond={df_hbonds.shape}, life={df_lifetimes.shape}, dist={df_distance.shape}, engaged={df_engaged.shape}")

    # Build color dict
    all_species_ids = sorted(set(all_species_ids))
    palette = sns.color_palette("tab10", len(all_species_ids))
    color_dict = {}
    for i, sp in enumerate(all_species_ids):
        color_dict[sp] = palette[i % len(palette)]

    # 2×2 figure
    two_by_two_plots_box_violin(
        df_lifetimes=df_lifetimes,
        df_distance=df_distance,
        df_engaged=df_engaged,
        df_hbonds=df_hbonds,
        color_dict=color_dict,
        output_prefix="intra_hbonds",
        figure_size=(16,12)
    )

    # ---------------------------------------------
    # Now we do the summary plot => existence heatmap
    # Then the new phosphate group violin
    # ---------------------------------------------
    # For the summary, we build occurrence_data + existence_data
    occurrence_data.clear()
    existence_data.clear()

    # We'll also collect all oxygens across these species
    all_oxys = set()
    for fname, species_id, n_ones in folder_list:
        log_path = process_path / fname / "analyze_final_sim" / "h_bonds" / "intra_hb.log"
        if not log_path.is_file():
            continue
        bdf = parse_hbond_log_to_dataframe(str(log_path))
        if bdf.empty:
            continue
        all_oxys.update(bdf['donor'].unique())
        all_oxys.update(bdf['acceptor'].unique())

    for fname, species_id, n_ones in folder_list:
        fpath = process_path / fname
        xpm_file = fpath / "analyze_final_sim" / "h_bonds" / "intra_hb_matrix.xpm"
        log_file = fpath / "analyze_final_sim" / "h_bonds" / "intra_hb.log"
        if not (xpm_file.is_file() and log_file.is_file()):
            continue
        bond_df = parse_hbond_log_to_dataframe(str(log_file))
        if bond_df.empty:
            continue
        mat, meta = parse_xpm(str(xpm_file))
        an = analyze_hydrogen_bonds(mat, meta)
        hpi = an['hbonds_per_index']

        # This time we sum donor+acceptor => occurrence_data
        from_both = count_oxygen_occurrences_from_matrix(bond_df, hpi)
        occurrence_data[species_id] = from_both

        # Build existence map
        bin_size_ns = 0.2
        frames_per_bin = int(bin_size_ns / time_per_frame_ns)
        if frames_per_bin<1:
            continue
        n_frames = mat.shape[1]
        nbins = n_frames // frames_per_bin
        if n_frames % frames_per_bin != 0:
            nbins+=1

        binned = np.zeros((mat.shape[0], nbins))
        for i in range(nbins):
            st = i*frames_per_bin
            en = min((i+1)*frames_per_bin, n_frames)
            chunk = mat[:, st:en]
            binned[:, i] = np.mean(chunk, axis=1)

        # Bond -> oxygen
        oxy_map = {}
        for idx, row in bond_df.iterrows():
            d = row['donor']
            a = row['acceptor']
            oxy_map.setdefault(d, []).append(idx)
            oxy_map.setdefault(a, []).append(idx)

        # aggregator
        sorted_oxy = sorted(all_oxys, key=lambda x: int(re.findall(r'\d+', x)[0]) if re.findall(r'\d+', x) else 0)
        agg = np.zeros((len(sorted_oxy), nbins))
        for o_i, oxy in enumerate(sorted_oxy):
            b_list = oxy_map.get(oxy, [])
            if b_list:
                for bn in range(nbins):
                    agg[o_i,bn] = np.max(binned[b_list,bn])
        existence_data[species_id] = agg


        # Summary plot
    
    from functools import partial
    # pass donor_atoms_by_species => color red in y-axis
    generate_summary_plot(
        occurrence_data, 
        existence_data,
        donor_atoms_by_species=donor_atoms_by_species,
        color_dict=color_dict,
        time_per_frame_ns=time_per_frame_ns,
        output_file='oxygen_occurrences_summary.pdf',
        folder_order=sorted_species_list
    )

    # ---------------------------------------------
    # Finally, create the new phosphate-group violin
    # ---------------------------------------------
    df_phosphate = pd.DataFrame(df_phosphate_rows)
    print(df_phosphate)
    # df_phosphate columns => [species, num_ones, donor_atom, phosphate_group, engaged_ns]
    # We do a single figure with x=phosphate_group, y=engaged_ns, hue=species
    plot_phosphate_group_violin(
        df_phosphate=df_phosphate,
        color_dict=color_dict,
        output_file="phosphate_group_violin.pdf"
    )

    # Finally, print the textual table
    for species_id in donor_acceptor_summaries:
        data = donor_acceptor_summaries[species_id]
        print_donor_acceptor_table(species_id, data)


    logger.info("All analysis steps completed successfully!")


if __name__ == "__main__":
    main()


