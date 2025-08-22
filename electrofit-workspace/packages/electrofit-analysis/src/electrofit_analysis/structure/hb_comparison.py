#!/usr/bin/env python3
"""
Combined Script for H-Bond Analysis with:
 1) A 2×2 Figure (Lifetime, Distance, Time/Atom, H-Bond Count)
 2) A Summary Plot (Occurrence bars + Existence Heatmap)
 3) A Phosphate Group Violin Plot (Donor occurrences by P1..P6)
 4) A textual table summarizing donor→acceptor phosphate usage.
    - additionally, a directed graph (network) showing donor→target relationships.

We assume IP6 microprotonation states, with donor→phosphate mapping:
  - P1: O1, O7, O8, O9
  - P2: O2, O10, O11, O12
  - P3: O3, O13, O14, O15
  - P4: O4, O16, O17, O18
  - P5: O5, O19, O20, O21
  - P6: O6, O22, O23, O24

Author: Arthur Laux
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
import networkx as nx

# ------------------------------------------------------------------------------
# Logging & Font Setup
# ------------------------------------------------------------------------------
logging.basicConfig(
    filename='analysis.log',
    filemode='w',  
    level=logging.INFO,
    format='%(asctime)s [%(levelname)s] %(message)s'
)
logger = logging.getLogger(__name__)

# Use Times New Roman + "talk" context from Seaborn
#mpl.rcParams['font.family'] = 'Times New Roman'
sns.set_context("talk")
sns.set_style({'font.family':'serif', 'font.serif':'Times New Roman'})
# ------------------------------------------------------------------------------
# Donor & Acceptor → Phosphate Mapping
# ------------------------------------------------------------------------------
donor_to_phosphate = {
    'O1': 'P1','O7': 'P1','O8': 'P1','O9': 'P1',
    'O2': 'P2','O10': 'P2','O11': 'P2','O12': 'P2',
    'O3': 'P3','O13': 'P3','O14': 'P3','O15': 'P3',
    'O4': 'P4','O16': 'P4','O17': 'P4','O18': 'P4',
    'O5': 'P5','O19': 'P5','O20': 'P5','O21': 'P5',
    'O6': 'P6','O22': 'P6','O23': 'P6','O24': 'P6'
}
acceptor_to_phosphate = donor_to_phosphate  # same mapping if acceptors are in [O1..O24]


# ------------------------------------------------------------------------------
# Common Utility Functions (PLACEHOLDER)
# ------------------------------------------------------------------------------
def refine_atom_name(atom):
    """
    If we see e.g. 'O' => 'O1', or 'O10' => 'O11'.
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
    """Ascend from current_dir until we find folder named project_name."""
    logger.info(f"Locating project root '{project_name}' from: {current_dir}")
    root = None
    while True:
        parent = os.path.dirname(current_dir)
        if os.path.basename(current_dir) == project_name:
            root = current_dir
        if parent == current_dir:  
            if root is None:
                raise FileNotFoundError(f"Project root '{project_name}' not found.")
            logger.info(f"Found project root: {root}")
            return root
        current_dir = parent

def load_hb_num_xvg(filename):
    """
    (#H-bonds vs time) from .xvg => Nx2 array.
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
    """Parses XPM => binary array (#FF0000 =>1). Returns (arr, metadata)."""
    logger.info(f"Parsing XPM: {file_path}")
    with open(file_path, 'r') as file:
        lines = file.readlines()

    metadata={}
    data_lines=[]
    color_map={}
    header_found=False
    data_started=False
    width=height=num_colors=chars_per_pixel=None

    for idx, raw_line in enumerate(lines, start=1):
        line=raw_line.strip()
        if line.startswith("/*") and not data_started:
            comment = line.strip("/* ").strip(" */")
            if ':' in comment:
                key, value = comment.split(":",1)
                metadata[key.strip().lower()] = value.strip().strip('"')
            continue

        if line.startswith('static char'):
            continue

        if (not header_found) and line.startswith('"'):
            line_str=line.strip(',').lstrip('"').rstrip('"')
            tokens=line_str.split()
            if len(tokens)>=4:
                try:
                    width,height,num_colors,chars_per_pixel=map(int,tokens[:4])
                    header_found=True
                    logger.info(f"XPM header => w={width},h={height},colors={num_colors},cpp={chars_per_pixel}")
                    continue
                except ValueError:
                    logger.error(f"Invalid header line {idx}: {raw_line}")
                    raise
            else:
                logger.error(f"Header line {idx} missing tokens: {raw_line}")
                raise ValueError(f"Bad XPM header line: {raw_line}")
        if header_found and not data_started:
            cdef=line.strip(',').lstrip('"').rstrip('"')
            pattern=rf'(.{{{chars_per_pixel}}})\s+c\s+(\S+)'
            match=re.match(pattern, cdef)
            if match:
                symbol,color_val=match.groups()
                color_map[symbol]=color_val
            if len(color_map)==num_colors:
                data_started=True
            continue
        if data_started and line.startswith('"'):
            row_data=line.strip(',').lstrip('"').rstrip('"')
            data_lines.append(row_data)

    if width is None or height is None:
        logger.error("No valid XPM header found.")
        raise ValueError("XPM header missing.")
    if len(data_lines)!=height:
        logger.warning(f"XPM mismatch: expected {height} lines,found {len(data_lines)}")

    arr=np.zeros((height,width), dtype=int)
    for y,row_str in enumerate(data_lines):
        for x,ch in enumerate(row_str):
            color=color_map.get(ch,None)
            if color=='#FF0000':
                arr[y,x]=1
    logger.info(f"XPM matrix => shape={arr.shape}")
    return arr,metadata

def analyze_hydrogen_bonds(data_matrix, metadata):
    """Summaries => over_time, per_index, lifetimes."""
    logger.info("Analyzing hydrogen bond matrix.")
    hbonds_over_time=np.sum(data_matrix,axis=0)
    hbonds_per_index=np.sum(data_matrix,axis=1)
    lifetimes=[]
    for row in data_matrix:
        cur=0
        run_list=[]
        for val in row:
            if val==1:
                cur+=1
            else:
                if cur>0:
                    run_list.append(cur)
                    cur=0
        if cur>0:
            run_list.append(cur)
        lifetimes.append(run_list)
    return {
        'hbonds_over_time':hbonds_over_time,
        'hbonds_per_index':hbonds_per_index,
        'lifetimes':lifetimes
    }

def parse_hbond_log_to_dataframe(file_path):
    """
    GROMACS .log => DataFrame [idx, donor, acceptor], refining numeric bits
    """
    logger.info(f"Parsing H-bond log: {file_path}")
    pairs=[]
    line_pattern=re.compile(r'^\s*(\S+)\s+-\s+(\S+)\s*$')
    atom_pattern=re.compile(r'^[A-Za-z]+\d+([A-Za-z]+\d*)$')
    with open(file_path,'r') as f:
        for line_no,raw_line in enumerate(f,start=1):
            line=raw_line.strip()
            if not line or line.startswith('#') or line.startswith('"""') or line.startswith('*'):
                continue
            match=line_pattern.match(line)
            if match:
                d_full,a_full=match.groups()
                d_match=atom_pattern.match(d_full)
                a_match=atom_pattern.match(a_full)
                if d_match and a_match:
                    donor_atom=refine_atom_name(d_match.group(1))
                    acceptor_atom=refine_atom_name(a_match.group(1))
                    pairs.append({'donor':donor_atom,'acceptor':acceptor_atom})
                else:
                    logger.warning(f"Line {line_no}: parse fail => {line}")
            else:
                logger.warning(f"Line {line_no}: no match => {line}")
    df=pd.DataFrame(pairs)
    df.reset_index(inplace=True)
    df.rename(columns={'index':'idx'},inplace=True)
    df['idx']=df.index
    logger.info(f"Found {len(df)} H-bond pairs.")
    return df

def count_donor_occurrences_from_matrix(hbond_df, hbonds_per_index):
    """
    Sums frames for 'donor' => avoid double-count. Returns Series: donor->frames
    """
    logger.info("Counting DONOR-based frames only.")
    df=hbond_df.copy()
    if 'count' in df.columns:
        df.drop(columns=['count'],inplace=True)
    df['count']=hbonds_per_index[df['idx'].values]
    melted=df.melt(
        id_vars=['idx','count'],
        value_vars=['donor'],
        value_name='oxygen'
    ).dropna(subset=['oxygen'])
    grouped=melted.groupby('oxygen')['count'].sum()
    logger.info(f"Donor-based occurrences => {len(grouped)} donors.")
    return grouped

def count_oxygen_occurrences_from_matrix(hbond_df, hbonds_per_index):
    """
    Summation for both donor & acceptor => used for final existence heatmap
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

def expand_phosphate_data_for_violins(df_phosphate):
    """
    Expand each row of df_phosphate (which has columns e.g.:
       [species, phosphate_group, engaged_ns]
    into multiple rows, effectively simulating a distribution
    for each (species, phosphate_group).

    Returns a new DataFrame with columns [species, pgroup_index].
    'pgroup_index' is repeated ~ engaged_ns times for each row.
    """
    # 1) Map each phosphate group to an integer index
    #    e.g. P1->0, P2->1, ... P6->5
    pgroup_index_map = {
        'P1': 0,
        'P2': 1,
        'P3': 2,
        'P4': 3,
        'P5': 4,
        'P6': 5
    }

    expanded_rows = []

    for row in df_phosphate.itertuples(index=False):
        species_id = row.species
        pgroup     = row.phosphate_group
        engaged_ns = row.engaged_ns * 1000

        # choose how many times to repeat. E.g. round it:
        repeat_count = int(round(engaged_ns))  
        # If engaged_ns is large, repeat_count could be big => large DataFrame.

        idx_val = pgroup_index_map.get(pgroup, None)
        if idx_val is None:
            continue  # skip unknown group

        # expand
        for _ in range(repeat_count):
            expanded_rows.append({
                'species': species_id,
                'pgroup_index': idx_val
            })

    expanded_df = pd.DataFrame(expanded_rows)
    return expanded_df

import matplotlib.pyplot as plt
import networkx as nx

def draw_phosphorus_diagram(species_id, top_three_map, edge_weights=None):
    """
    Draws a directed graph (P1..P6) showing donor->target relationships,
    with the species label centered and pure‐acceptor nodes in white.
    """
    if edge_weights is None:
        edge_weights = {}

    # 1) Build graph
    G = nx.DiGraph()
    for i in range(1, 7):
        G.add_node(i)
    for donor, ranks in top_three_map.items():
        for rank, target in ranks.items():
            G.add_edge(donor, target, rank=rank)

    # 2) Custom circular layout (rotated so P1 sits where P3 normally is)
    oldpos = nx.circular_layout(G)
    pos = {1: oldpos[3], 2: oldpos[4], 3: oldpos[5],
           4: oldpos[6], 5: oldpos[1], 6: oldpos[2]}

    # 3) Determine node colors: 
    #    white for nodes with in_degree>0 and out_degree==0 ("pure acceptors"), else lightgray
    node_colors = []
    for n in G.nodes():
        if G.in_degree(n) > 0 and G.out_degree(n) == 0:
            node_colors.append("white")
        else:
            node_colors.append("lightgray")

    # 4) Prepare figure
    fig, ax = plt.subplots(figsize=(8,6))

    # 5) Draw nodes
    nx.draw_networkx_nodes(
        G, pos, node_color=node_colors,
        edgecolors="black", node_size=1300, ax=ax
    )

    # 6) Split edges by rank
    edges_first  = [(u,v) for u,v,d in G.edges(data=True) if d['rank']=="first"]
    edges_second = [(u,v) for u,v,d in G.edges(data=True) if d['rank']=="second"]
    edges_third  = [(u,v) for u,v,d in G.edges(data=True) if d['rank']=="third"]

    #  – edge-width scaling
    min_w, max_w = 1.0, 6.0
    if edge_weights:
        max_val = max(edge_weights.values()) or 1.0
        def scale(t):
            return min_w + (max_w - min_w) * (t / max_val)
    else:
        scale = lambda t: min_w

    # 7) Draw first‐tier (red), handling reciprocals
    mutual = [(u,v) for (u,v) in edges_first if (v,u) in edges_first and u<v]
    mutual_set = set(mutual) | set((v,u) for (u,v) in mutual)
    nonmutual = [e for e in edges_first if e not in mutual_set]

    nx.draw_networkx_edges(
        G, pos, edgelist=nonmutual,
        width=[scale(edge_weights.get(e,0.0)) for e in nonmutual],
        edge_color="red", arrows=True, arrowstyle='-|>', arrowsize=20,
        connectionstyle="arc3,rad=0.0",
        min_source_margin=45, min_target_margin=45, ax=ax
    )
    for u,v in mutual:
        w_uv = scale(edge_weights.get((u,v),0.0))
        nx.draw_networkx_edges(
            G, pos, edgelist=[(u,v)],
            width=[w_uv], edge_color="red", arrows=True,
            arrowstyle='-|>', arrowsize=20,
            connectionstyle="arc3,rad=+0.2",
            min_source_margin=45, min_target_margin=45, ax=ax
        )
        w_vu = scale(edge_weights.get((v,u),0.0))
        nx.draw_networkx_edges(
            G, pos, edgelist=[(v,u)],
            width=[w_vu], edge_color="red", arrows=True,
            arrowstyle='-|>', arrowsize=20,
            connectionstyle="arc3,rad=+0.2",
            min_source_margin=45, min_target_margin=45, ax=ax
        )

    # 8) Second‐tier (black)
    nx.draw_networkx_edges(
        G, pos, edgelist=edges_second,
        width=[scale(edge_weights.get(e,0.0)) for e in edges_second],
        edge_color="black", arrows=True, arrowstyle='-|>', arrowsize=20,
        connectionstyle="arc3,rad=+0.2",
        min_source_margin=45, min_target_margin=45, ax=ax
    )

    # 9) Third‐tier (blue)
    nx.draw_networkx_edges(
        G, pos, edgelist=edges_third,
        width=[scale(edge_weights.get(e,0.0)) for e in edges_third],
        edge_color="blue", arrows=True, arrowstyle='-|>', arrowsize=20,
        connectionstyle="arc3,rad=-0.2",
        min_source_margin=45, min_target_margin=45, ax=ax
    )

    # 10) Draw node labels
    nx.draw_networkx_labels(
        G, pos, labels={i: f"P{i}" for i in G.nodes()}, font_size=14, ax=ax
    )

    # 11) Add centered species label in the middle of the plot
    #    Use axis‐coordinates (0.5,0.5) so it's always at the visual center
    ax.text(
        0.5, 0.5, f"({species_id})",
        fontsize=18, fontweight="bold",
        ha="center", va="center",
        transform=ax.transAxes,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="black", alpha=0)
    )

    ax.set_axis_off()
    plt.tight_layout()
    plt.savefig(f"species_{species_id}.pdf")
    plt.close()


# ------------------------------------------------------------------------------
# 1) 2×2 Figure (Lifetime, Distance, Time/Atom, H-Bond Count)
# ------------------------------------------------------------------------------
def three_by_two_plots_box_violin(
    df_lifetimes,   
    df_distance,    
    df_engaged,     
    df_hbonds,      
    df_phosphate,   # [species, phosphate_group, engaged_ns, ...]
    color_dict=None,
    output_prefix="intra_hbonds",
    figure_size=(16, 22),   
    spacer_prefix="~space_"
):
    """
    Creates a 3×2 layout:

      Row 0: (0,0) => (a) Lifetime (Box)  
              (0,1) => (b) Distance (Box)

      Row 1: (1,0) => (c) Time/Atom (Stacked)
              (1,1) => (d) H-bond Count (Violin)

      Row 2: (2,0) => (e) Horizontal Phosphate Violin
              (2,1) => Empty (white)

    We do not alter the 4 subplots from your 2×2 figure, 
    but add a 5th plot (phosphate violin) in row=2, col=0.
    Spacers are injected so all y-axis labels match across the subplots.
    """

    logger.info("Generating a 3×2 figure: the original 4 subplots + a 5th horizontal violin below them.")

    # ------------------------------------------------------
    # 0) Check empties for the top-4 data
    # ------------------------------------------------------
    all_empty = True
    for df_check in [df_lifetimes, df_distance, df_engaged, df_hbonds]:
        if df_check is not None and not df_check.empty:
            all_empty = False
            break
    if all_empty:
        logger.warning("All top-4 DataFrames empty => skipping top-4 plots.")
        # We still allow the 5th plot if df_phosphate is non-empty

    # ------------------------------------------------------
    # 1) Gather species from the top-4 DataFrames
    # ------------------------------------------------------
    species_in_any = set()
    for df_check in [df_lifetimes, df_distance, df_engaged, df_hbonds]:
        if df_check is not None and not df_check.empty:
            species_in_any.update(df_check["species"].unique())

    # Helper: figure out how many '1' bits for a species
    def get_num_ones(sp, df_list):
        for d in df_list:
            if d is not None and not d.empty:
                row = d.loc[d["species"] == sp]
                if not row.empty:
                    return row.iloc[0]["num_ones"]
        return None

    # Keep only species with 3,4,5
    all_with_num = []
    for sp in species_in_any:
        val = get_num_ones(sp, [df_lifetimes, df_distance, df_engaged, df_hbonds])
        if val in [3, 4, 5]:
            all_with_num.append((sp, val))

    if not all_with_num:
        logger.warning("No species with num_ones in [3,4,5]. Skipping all top-4 plots + 5th.")
        return

    # Sort => 5->top, 4->middle, 3->bottom
    all_with_num.sort(key=lambda x: (-x[1], x[0]))

    # Insert spacers for transitions (5->4->3)
    grouped_species_order = []
    prev_n_ones = None
    for i, (sp, n_ones) in enumerate(all_with_num):
        if i > 0 and n_ones != prev_n_ones:
            grouped_species_order.append((f"{spacer_prefix}{prev_n_ones}to{n_ones}", None))
        grouped_species_order.append((sp, n_ones))
        prev_n_ones = n_ones
    # final list => species + placeholders
    full_order = [t[0] for t in grouped_species_order]

    # A helper to inject spacers into a DataFrame
    def inject_spacers(df, measure_col):
        """
        For each spacer in grouped_species_order, add a row with measure_col=NaN.
        """
        if df is None or df.empty:
            return pd.DataFrame(columns=["species", measure_col, "num_ones"])
        new_df = df.copy()
        row_spacers = []
        for sp, val in grouped_species_order:
            if val is None:  # a dummy spacer
                row_spacers.append({
                    "species": sp,
                    measure_col: np.nan,
                    "num_ones": -1
                })
        if row_spacers:
            dummy_df = pd.DataFrame(row_spacers)
            combined = pd.concat([new_df, dummy_df], ignore_index=True)
            return combined
        return new_df

    # ------------------------------------------------------
    # 2) For df_distance, expand by frequency if needed
    # ------------------------------------------------------
    if df_distance is not None and not df_distance.empty:
        df_distance = df_distance.copy()
        df_distance['frequency'] = df_distance['frequency'].astype(int)
        df_dist_expanded = df_distance.loc[df_distance.index.repeat(df_distance['frequency'])].reset_index(drop=True)
    else:
        df_dist_expanded = pd.DataFrame()

    # Now inject spacers for each measure
    df_life2  = inject_spacers(df_lifetimes, "lifetime_ns")
    df_dist2  = inject_spacers(df_dist_expanded, "distance_nm")
    df_eng2   = inject_spacers(df_engaged, "engaged_time_ns")
    df_hbond2 = inject_spacers(df_hbonds, "hbond_count")

    # ------------------------------------------------------
    # 3) Summarize df_phosphate => row=species, col=P1..P6, so we can do time/atom stacked
    #    (unchanged from your code)
    # ------------------------------------------------------
    grouped_phos = df_phosphate.groupby(['species','phosphate_group'])['engaged_ns'].sum().reset_index()
    pivot_phos = grouped_phos.pivot(index='species', columns='phosphate_group', values='engaged_ns').fillna(0)

    phosphate_cols = ['P1','P2','P3','P4','P5','P6']
    for c in phosphate_cols:
        if c not in pivot_phos.columns:
            pivot_phos[c] = 0.0
    pivot_phos = pivot_phos[phosphate_cols]

    # If no color_dict => create an empty one
    if color_dict is None:
        color_dict = {}
    else:
        color_dict = dict(color_dict)  # copy to avoid mutating the original

    # Ensure spacers have a color
    spacer_color = "#D3D3D3"
    for sp, val in grouped_species_order:
        if val is None and sp not in color_dict:
            color_dict[sp] = spacer_color

    # ------------------------------------------------------
    # 4) Create a 3×2 figure 
    # ------------------------------------------------------
    fig, axes = plt.subplots(nrows=3, ncols=2, figsize=figure_size)

    # Ax references
    ax_life   = axes[0][0]  # (a)
    ax_dist   = axes[0][1]  # (b)
    ax_stack  = axes[1][0]  # (c)
    ax_hbond  = axes[1][1]  # (d)
    ax_violin = axes[2][0]  # (e)
    ax_blank  = axes[2][1]  
    ax_blank.axis('off')   # White cell

    # ~~~~~~~~~~~~~ (a) Lifetime (Box) ~~~~~~~~~~~~~
    if df_life2 is not None and not df_life2.empty:
        sns.boxplot(
            data=df_life2,
            x="lifetime_ns",
            y="species",
            order=full_order,
            orient='h',
            palette=color_dict,
            showmeans=True,
            meanprops={'marker':'o','markerfacecolor':'black','markeredgecolor':'black','markersize':8},
            showfliers=True,
            ax=ax_life
        )
        ax_life.set_title("(a) Lifetime (Box)")
        ax_life.set_xlabel("Lifetime (ns)")
        ax_life.set_ylabel("Species")
        ax_life.set_xscale("log")

        # Hide spacer row labels
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

    # ~~~~~~~~~~~~~ (b) Distance (Box) ~~~~~~~~~~~~~
    if df_dist2 is not None and not df_dist2.empty:
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
        ax_dist.set_title("(b) Distance (Box)")
        ax_dist.set_xlabel("Distance (nm)")
        ax_dist.set_ylabel("")
        ax_dist.set_yticklabels([])
    else:
        ax_dist.set_title("No Distance Data")
        ax_dist.set_xlabel("")
        ax_dist.set_ylabel("")
        ax_dist.set_yticklabels([])

    # ~~~~~~~~~~~~~ (d) H-bond Count (Violin) ~~~~~~~~~~~~~
    if df_hbond2 is not None and not df_hbond2.empty:
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
        ax_hbond.set_title("(d) H-bond Count (Violin)")
        ax_hbond.set_xlabel("H-bond Count")
        ax_hbond.set_ylabel("")
        ax_hbond.set_yticklabels([])

        # Overlay means
        df_means = df_hbond2.groupby('species')['hbond_count'].mean().reset_index()
        df_means = df_means[~df_means['species'].str.startswith(spacer_prefix)]
        ax_hbond.scatter(
            df_means['hbond_count'],
            df_means['species'],
            color='black',
            marker='o',
            s=80,
            label='Mean',
            zorder=5
        )
    else:
        ax_hbond.set_title("No H-bond Data")
        ax_hbond.set_xlabel("")
        ax_hbond.set_ylabel("")
        ax_hbond.set_yticklabels([])

    # ~~~~~~~~~~~~~ (c) Time/Atom (Stacked) ~~~~~~~~~~~~~
    ax_stack.set_title("(c) Time / H-Atom")
    ax_stack.set_xlabel("Engaged Time (ns)")
    ax_stack.set_ylabel("Species")

    # We'll map each species/spacer -> a y-position
    y_positions = {}
    real_index = 0
    for sp, val in grouped_species_order:
        y_positions[sp] = real_index
        real_index += 1

    # For each species, sum of engaged time in pivot_phos => stacked bars
    for sp in y_positions:
        if sp not in pivot_phos.index:
            continue
        species_color = color_dict.get(sp, "gray")  
        row_vals = pivot_phos.loc[sp]
        left_val = 0.0
        for pgroup in ['P1','P2','P3','P4','P5','P6']:
            seg_val = row_vals[pgroup]
            if seg_val>0:
                ax_stack.barh(
                    y=y_positions[sp],
                    width=seg_val,
                    left=left_val,
                    height=0.8,
                    color=species_color, 
                    edgecolor='k',
                    alpha=0.9
                )
                if seg_val>12.5:
                    ax_stack.annotate(
                        f"{seg_val:.2f}",
                        (left_val + seg_val/2, y_positions[sp]),
                        ha='center', va='center',
                        color='white', fontsize=10
                    )
                left_val += seg_val

    ax_stack.set_ylim(-0.5, real_index-0.5)
    ax_stack.set_yticks(np.arange(real_index))
    # Build final y labels
    y_lbls_stack = []
    for sp, val in grouped_species_order:
        if val is None:
            y_lbls_stack.append("")
        else:
            y_lbls_stack.append(sp)
    ax_stack.set_yticklabels(y_lbls_stack)
    ax_stack.invert_yaxis()

    # ~~~~~~~~~~~~~ (e) Horizontal Phosphate Violin ~~~~~~~~~~~~~
    # We want to inject spacers for the expanded data, so the y-axis matches top plots
    expanded_df = expand_phosphate_data_for_violins(df_phosphate)
    if expanded_df.empty:
        logger.warning("No phosphate-group donor data => skipping subplot (e).")
        ax_violin.set_title("(e) No Phosphate Data")
    else:
        # We'll treat 'pgroup_index' as measure column so we can inject spacers
        # though they won't have pgroup_index in the spacer row. It's still okay.
        # We'll create a new column measure for injection
        expanded_for_injection = expanded_df.copy()
        expanded_for_injection['pgroup_index'] = expanded_for_injection['pgroup_index'].astype(float)

        # Now we do the same injection approach
        # We'll rename measure_col to 'pgroup_index'
        def inject_violin_spacers(df_in):
            if df_in is None or df_in.empty:
                return pd.DataFrame(columns=["species", "pgroup_index"])
            new_df = df_in.copy()
            row_spacers = []
            for sp, val in grouped_species_order:
                if val is None:  
                    row_spacers.append({
                        "species": sp,
                        "pgroup_index": np.nan
                    })
            if row_spacers:
                dummy_df = pd.DataFrame(row_spacers)
                return pd.concat([new_df, dummy_df], ignore_index=True)
            return new_df

        expanded_with_spacers = inject_violin_spacers(expanded_for_injection)

        # Now we can do a horizontal violin
        sns.violinplot(
            data=expanded_with_spacers,
            x='pgroup_index',
            y='species',
            hue='species',
            order=full_order,
            palette=color_dict,
            orient='h',
            cut=0,
            inner=None,
            ax=ax_violin
        )
        ax_violin.set_title("(e) Phosphate Group Horizontal Violin")

        # x-axis ticks => 0..5 => P1..P6
        phosphate_labels = ['P1','P2','P3','P4','P5','P6']
        ax_violin.set_xticks(range(len(phosphate_labels)))
        ax_violin.set_xticklabels(phosphate_labels)
        ax_violin.set_xlabel("Phosphate Group")
        ax_violin.set_ylabel("Species")

        # We'll create the same y ticks as the other subplots => full_order
        ax_violin.set_yticks(np.arange(real_index))
        y_lbls_violin = []
        for sp, val in grouped_species_order:
            if val is None:
                y_lbls_violin.append("")
            else:
                y_lbls_violin.append(sp)
        ax_violin.set_yticklabels(y_lbls_violin)

        # If the legend is large, you may remove or place it differently
        # ax_violin.legend(loc='best')

    # ------------------------------------------------------
    # 5) Final
    # ------------------------------------------------------
    plt.tight_layout()
    out_pdf = f"{output_prefix}_3x2_with_violin.pdf"
    plt.savefig(out_pdf, dpi=300)
    plt.close()
    logger.info(f"Saved 3×2 figure with the 5th plot => '{out_pdf}'!")

# ------------------------------------------------------------------------------
# 2) Summary Plot (Occurrence bars + Existence Heatmap)
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

    # ----------------------------------------------------------------------
    # 1) Remove 'OW1' from each species' occurrence_data
    # ----------------------------------------------------------------------
    for sp in list(occurrence_data.keys()):
        ser = occurrence_data[sp]
        if "OW1" in ser.index:
            print(ser)
            print(ser.index)
            occurrence_data[sp] = ser.drop("OW1")

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
# 4) Textual Table: Donor→Acceptor usage
# ------------------------------------------------------------------------------

def build_donor_acceptor_summary(species_id, bond_df, hbonds_per_index, time_per_frame_ns, donor_map, acceptor_map):
    """
    For each bond => map donor->pgroup, acceptor->pgroup => accumulate times
    """
    result = {}
    hbonds_length = len(hbonds_per_index)  # Cache length for efficiency
    for i, row in bond_df.iterrows():
        b_idx = row['idx']
        d_atom = row['donor']
        a_atom = row['acceptor']
        
        # Safe access to hbonds_per_index
        if 0 <= b_idx < hbonds_length:
            frames = hbonds_per_index[b_idx]
        else:
            frames = 0  # Default to 0 if index is out of bounds
        
        engaged_ns = frames * time_per_frame_ns

        donor_pgroup = donor_map.get(d_atom, None)
        acceptor_pgroup = acceptor_map.get(a_atom, None)
        if donor_pgroup is None or acceptor_pgroup is None:
            continue

        if donor_pgroup not in result:
            result[donor_pgroup] = {"sum_donor_time": 0.0, "targets": {}}
        result[donor_pgroup]["sum_donor_time"] += engaged_ns

        key_t = (acceptor_pgroup, a_atom)
        if key_t not in result[donor_pgroup]["targets"]:
            result[donor_pgroup]["targets"][key_t] = 0.0
        result[donor_pgroup]["targets"][key_t] += engaged_ns
    return result

def print_donor_acceptor_table(species_id, data, logger, summary_format='arrow', draw_figure=True):
    """
    Generate a summary table for each donor pgroup, log it, and (optionally) draw a figure
    for the top 3 phosphorus-target relationships (1st=red, 2nd=black, 3rd=blue).
    Now we also produce separate textual summaries for first, second, and third P-targets.

    Parameters:
    - species_id (str): Identifier for the species (e.g. "101111").
    - data (dict): Summary data from build_donor_acceptor_summary.
    - logger (logging.Logger): Configured logger instance.
    - summary_format (str): Format of the summary line ('arrow', 'graphical', etc.).
    - draw_figure (bool): If True, generate a color-coded figure for each species.
    """
    table_lines = []
    table_lines.append(f"Species: {species_id}\n")

    # We'll store up to 3 P-targets for each donor in top_three_map => used by the figure
    # Also store them in rank_p_targets => used by textual summary
    top_three_map   = {}  # e.g. {1: {"first":2, "second":4, "third":5}, ...}
    rank_p_targets  = { "first": {}, "second": {}, "third": {} }  
    # e.g. rank_p_targets["first"] = { "P2":"P3", "P4":"P5", ... }

    donor_sorted = sorted(data.keys())  # e.g. ['P1','P2','P3','P4','P5','P6'] if they exist
    for dpgroup in donor_sorted:
        sum_time = data[dpgroup]["sum_donor_time"]
        table_lines.append(f"  Donor {dpgroup} => total time: {sum_time:.2f} ns")
        
        # Sort all acceptor targets by descending time
        t_map = data[dpgroup]["targets"]  # e.g. {(acc_pg, acc_atom): time_ns, ...}
        t_list = sorted(t_map.items(), key=lambda x: x[1], reverse=True)

        if not t_list:
            table_lines.append("    No acceptor targets found.")
            table_lines.append("    No Phosphorus acceptor targets found.\n")
            continue
        
        # Print each target with rank (1st, 2nd, 3rd, etc.) for logging
        for i, ((acc_pg, acc_atom), val_ns) in enumerate(t_list, start=1):
            # ordinal suffix
            if 10 <= i % 100 <= 20:
                suffix = 'th'
            else:
                suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(i % 10, 'th')
            target_label = f"{i}{suffix} target"
            percentage = (val_ns / sum_time * 100) if sum_time > 0 else 0.0
            table_lines.append(f"    {target_label:15s} => {acc_pg}({acc_atom}) = {percentage:.1f} %")

        # Aggregate phosphorus targets only
        p_group_summary = {}
        for ((acc_pg, acc_atom), val_ns) in t_list:
            if acc_pg.startswith('P'):
                p_group_summary[acc_pg] = p_group_summary.get(acc_pg, 0.0) + val_ns

        if p_group_summary:
            # Sort P-targets by time, descending
            p_targets_sorted = sorted(p_group_summary.items(), key=lambda x: x[1], reverse=True)
            table_lines.append("    Phosphorus Targets:")

            # Donor numeric index (P3 => 3)
            donor_num = int(dpgroup.replace('P',''))
            top_three_map[donor_num] = {}

            # We’ll store up to 3 in top_three_map: 'first','second','third'
            rank_to_key = {1: 'first', 2: 'second', 3: 'third'}

            for i, (acc_pg, total_ns) in enumerate(p_targets_sorted, start=1):
                if 10 <= i % 100 <= 20:
                    suffix = 'th'
                else:
                    suffix = {1: 'st', 2: 'nd', 3: 'rd'}.get(i % 10, 'th')
                p_target_label = f"{i}{suffix} P target"

                p_percentage = (total_ns / sum_time * 100) if sum_time > 0 else 0.0
                table_lines.append(f"      {p_target_label:15s} => {acc_pg} = {p_percentage:.1f} %")

                # If i <= 3 => store in top_three_map for the figure
                if i <= 3:
                    top_three_map[donor_num][rank_to_key[i]] = int(acc_pg.replace('P',''))

                # Also store in rank_p_targets => for the textual summary
                if i <= 3:
                    # e.g. rank_p_targets["first"][dpgroup] = acc_pg
                    rank_str = rank_to_key[i]
                    rank_p_targets[rank_str][dpgroup] = acc_pg
        else:
            table_lines.append("    No Phosphorus acceptor targets found.")
        
        table_lines.append("")  # spacing

    # Now we produce three separate summaries (first, second, third)
    # e.g. "Summary of First Phosphorus Targets:\n  P2 -> P3 | P4 -> P5 | ..."
    # If any donors have a first, second, third P-target.

    # Helper to produce lines like "P2 -> P3 | P4 -> P5 | ..."
    def build_summary_line(mapping):
        # mapping e.g. rank_p_targets["first"] => { "P2":"P3", "P4":"P5", ... }
        pairs = []
        for donor, first_p in sorted(mapping.items()):
            pairs.append(f"{donor} -> {first_p}")
        return " | ".join(pairs)

    # 1) First
    if rank_p_targets["first"]:
        line_str = build_summary_line(rank_p_targets["first"])
        table_lines.append("Summary of First Phosphorus Targets:")
        table_lines.append(f"  {line_str}\n")

    # 2) Second
    if rank_p_targets["second"]:
        line_str = build_summary_line(rank_p_targets["second"])
        table_lines.append("Summary of Second Phosphorus Targets:")
        table_lines.append(f"  {line_str}\n")

    # 3) Third
    if rank_p_targets["third"]:
        line_str = build_summary_line(rank_p_targets["third"])
        table_lines.append("Summary of Third Phosphorus Targets:")
        table_lines.append(f"  {line_str}\n")

    table_lines.append("-" * 70)
    table_lines.append("")

    # Combine lines + log
    table_str = "\n".join(table_lines)
    logger.info(table_str)

    # Finally, if desired, draw the figure
    if draw_figure:
        # Build a simple phosphate→phosphate weight map in ns
        edge_weights = {}
        for donor_pg, info in data.items():
            u = int(donor_pg.replace('P',''))
            for (acceptor_pg, _atom), t_ns in info['targets'].items():
                if acceptor_pg.startswith('P'):
                    v = int(acceptor_pg.replace('P',''))
                    edge_weights[(u,v)] = edge_weights.get((u,v), 0.0) + t_ns
        # Now call the figure, passing this map:
        draw_phosphorus_diagram(species_id, top_three_map, edge_weights)
        #draw_phosphorus_diagram(species_id, top_three_map)

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
def main(hbond_type):
    logger.info("Starting combined analysis with 3x2 figure, summary plot, phosphate violin, textual table...")

    # 1) Find project root & process directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    try:
        project_path = find_project_root(script_dir)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)
        
    process_dir = os.path.join(project_path, "data/results/final/process_IP6_11ps")
    if not os.path.isdir(process_dir):
        logger.error(f"Invalid process directory => {process_dir}")
        sys.exit(1)
    logger.info(f"Using process directory => {process_dir}")

    process_path = Path(process_dir)
    time_per_frame_ns = 0.01

    # 2) Lists/dicts for data accumulation
    folder_list = []
    all_species_ids = []
    hbonds_data = []      # for 2x2 hbond count
    lifetimes_data = []   # for 2x2 lifetime
    distance_data = []    # for 2x2 distance
    engaged_data = []     # for 2x2 time/atom
    df_phosphate_rows = []  # for phosphate group violin
    donor_atoms_by_species = {}
    occurrence_data = {}
    existence_data = {}
    donor_acceptor_summaries = {}

    # 3) Identify species subfolders
    for fname in sorted(os.listdir(process_path)):
        if not fname.startswith("IP_"):
            continue
        folder_path = process_path / fname
        if not folder_path.is_dir():
            continue
        species_id = fname.replace("IP_", "")
        n_ones = species_id.count('1')
        if n_ones not in [3,4,5]:
            continue

        folder_list.append((fname, species_id, n_ones))
        all_species_ids.append(species_id)

    # sort 5->4->3
    folder_list.sort(key=lambda x: (-x[2], x[1]))
    sorted_species_list = [x[1] for x in folder_list]

    # 4) Loop over species => parse data
    for folder_name, species_id, n_ones in folder_list:
        folder_path = process_path / folder_name

        # Required files
        hb_num_file = folder_path / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb_num.xvg"
        xpm_file    = folder_path / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb_matrix.xpm"
        dist_file   = folder_path / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb_dist.xvg"
        log_file    = folder_path / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb.log"

        # Check existence
        if not (hb_num_file.is_file() and xpm_file.is_file() and dist_file.is_file() and log_file.is_file()):
            logger.warning(f"Skipping species {species_id}: missing hbond files.")
            continue

        # (a) H-bond count vs time
        data_num = load_hb_num_xvg(str(hb_num_file))
        if data_num.size<2:
            logger.warning(f"No hbond count data => {species_id}")
            continue
        hbond_counts = data_num[:,1]
        for val in hbond_counts:
            hbonds_data.append((species_id, n_ones, val))

        # (b) XPM => analyze => lifetimes
        xpm_matrix, meta = parse_xpm(str(xpm_file))
        analysis_res = analyze_hydrogen_bonds(xpm_matrix, meta)
        lf_frames = [lf for bond_lf in analysis_res['lifetimes'] for lf in bond_lf]
        lf_ns = [f*time_per_frame_ns for f in lf_frames]
        for val_ns in lf_ns:
            lifetimes_data.append((species_id, n_ones, val_ns))

        # (c) Distances
        d_list = []
        with open(dist_file, 'r') as df_:
            for line in df_:
                if not line.strip() or line.startswith(('#','@')):
                    continue
                parts = line.split()
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
        for row_val in arr_d:
            distance_data.append((species_id,n_ones,row_val[0], row_val[1]))

        # (d) bond_df => parse log
        bond_df = parse_hbond_log_to_dataframe(str(log_file))
        if bond_df.empty:
            logger.warning(f"No hbond pairs => {species_id}")
            continue
        # store donor sets
        donor_atoms_by_species[species_id] = set(bond_df['donor'].unique())

        # from analysis => frames/bond
        hpi = analysis_res['hbonds_per_index']
        # donor-based sum => time/atom
        donor_counts = count_donor_occurrences_from_matrix(bond_df, hpi)
        total_donor_frames = donor_counts.sum()
        engaged_ns = (total_donor_frames * time_per_frame_ns)/n_ones
        engaged_data.append((species_id, n_ones, engaged_ns))

        # Build phosphate group data => for violin
        for donor_atom, frames in donor_counts.items():
            engaged_time_ns = frames*time_per_frame_ns
            pgroup = donor_to_phosphate.get(donor_atom, None)
            if pgroup:
                df_phosphate_rows.append({
                    'species': species_id,
                    'num_ones': n_ones,
                    'donor_atom': donor_atom,
                    'phosphate_group': pgroup,
                    'engaged_ns': engaged_time_ns
                })

        # Build donor→acceptor summary => textual table
        summ_data = build_donor_acceptor_summary(
            species_id, bond_df, hpi, time_per_frame_ns,
            donor_to_phosphate, acceptor_to_phosphate
        )
        donor_acceptor_summaries[species_id] = summ_data

    # 5) Build DataFrames for the 2×2 figure
    df_hbonds    = pd.DataFrame(hbonds_data, columns=["species","num_ones","hbond_count"])
    df_lifetimes = pd.DataFrame(lifetimes_data, columns=["species","num_ones","lifetime_ns"])
    df_distance  = pd.DataFrame(distance_data, columns=["species","num_ones","distance_nm","frequency"])
    df_engaged   = pd.DataFrame(engaged_data,  columns=["species","num_ones","engaged_time_ns"])

    # Build color_dict
    all_species_ids = sorted(set(all_species_ids))
    palette = sns.color_palette("tab10", len(all_species_ids))
    color_dict = {}
    for i, sp in enumerate(all_species_ids):
        color_dict[sp] = palette[i % len(palette)]

    # Build a DataFrame for the phosphate group violin
    df_phosphate = pd.DataFrame(df_phosphate_rows)
    # 6) Optionally create the 3×2 figure
    three_by_two_plots_box_violin(
        df_lifetimes=df_lifetimes,
        df_distance=df_distance,
        df_engaged=df_engaged,
        df_hbonds=df_hbonds,
        df_phosphate=df_phosphate,
        color_dict=color_dict,
        output_prefix=f"{hbond_type}_hbonds",
        figure_size=(16,18)
    )

    # 7) Create data for summary existence plot
    occurrence_data.clear()
    existence_data.clear()

    # We'll also collect all oxygens across these species
    all_oxys = set()
    for fname, species_id, n_ones in folder_list:
        log_path = process_path / fname / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb.log"
        if not log_path.is_file():
            continue
        bdf = parse_hbond_log_to_dataframe(str(log_path))
        if bdf.empty:
            continue
        all_oxys.update(bdf['donor'].unique())
        all_oxys.update(bdf['acceptor'].unique())

    for fname, species_id, n_ones in folder_list:
        fpath = process_path / fname
        xpm_file = fpath / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb_matrix.xpm"
        log_file = fpath / "analyze_final_sim" / "h_bonds" / f"{hbond_type}_hb.log"
        if not (xpm_file.is_file() and log_file.is_file()):
            continue
        bond_df = parse_hbond_log_to_dataframe(str(log_file))
        if bond_df.empty:
            continue
        mat, meta = parse_xpm(str(xpm_file))
        an = analyze_hydrogen_bonds(mat, meta)
        hpi = an['hbonds_per_index']

        # Sum donor+acceptor → occurrence_data
        from_both = count_oxygen_occurrences_from_matrix(bond_df, hpi)
        occurrence_data[species_id] = from_both

        # Build existence map
        bin_size_ns = 0.2
        frames_per_bin = int(bin_size_ns / time_per_frame_ns)
        if frames_per_bin < 1:
            continue
        n_frames = mat.shape[1]
        nbins = n_frames // frames_per_bin
        if n_frames % frames_per_bin != 0:
            nbins += 1

        binned = np.zeros((mat.shape[0], nbins))
        for i in range(nbins):
            st = i * frames_per_bin
            en = min((i + 1) * frames_per_bin, n_frames)
            chunk = mat[:, st:en]
            binned[:, i] = np.mean(chunk, axis=1)

        # Build a mapping from oxygen to bond indices
        oxy_map = {}
        for idx, row in bond_df.iterrows():
            d = row['donor']
            a = row['acceptor']
            oxy_map.setdefault(d, []).append(idx)
            oxy_map.setdefault(a, []).append(idx)

        # Filter out "OW1" from the global oxygen list
        sorted_oxy = sorted([o for o in all_oxys if o != "OW1"],
                            key=lambda x: int(re.findall(r'\d+', x)[0]) if re.findall(r'\d+', x) else 0)
        
        agg = np.zeros((len(sorted_oxy), nbins))
        for o_i, oxy in enumerate(sorted_oxy):
            b_list = oxy_map.get(oxy, [])
            if b_list:
                for bn in range(nbins):
                    agg[o_i, bn] = np.max(binned[b_list, bn])
        existence_data[species_id] = agg

    
    # 8) Once we have occurrence_data + existence_data:
    # This is time consuming...
    generate_summary_plot(
        occurrence_data=occurrence_data,
        existence_data=existence_data,
        donor_atoms_by_species=donor_atoms_by_species,
        color_dict=color_dict,
        time_per_frame_ns=time_per_frame_ns,
        output_file=f"{hbond_type}_oxygen_occurrences_summary.pdf",
        folder_order=sorted_species_list
    )

    # 10) Print textual table
    for sp, data in donor_acceptor_summaries.items():
        print_donor_acceptor_table(sp, data, logger, "arrow")

    logger.info("All analysis steps completed successfully!")


if __name__=="__main__":
    main("intra")