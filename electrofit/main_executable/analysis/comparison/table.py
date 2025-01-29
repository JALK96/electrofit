#!/usr/bin/env python3
"""
Combined Script for H-Bond Analysis with:
  1) A 2×2 Figure (Lifetime, Distance, Time/Atom, H-Bond Count)
  2) A Summary Plot (Occurrence bars + Existence Heatmap)
  3) A NEW Phosphate Group Violin Plot:
     - Aggregating donor occurrences by phosphate label (P1..P6)
     - Each species (micro protonation state) shown as separate distribution
  4) A textual table showing donor→acceptor phosphate usage.

We focus on IP6 microprotonation states. 
If a donor is in [O1,O7,O8,O9] => P1, [O2,O10,O11,O12] => P2, ...,
in line with the final corrected dictionary.

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

# Acceptors use the same dictionary, if they follow the same O1..O24 pattern
acceptor_to_phosphate = donor_to_phosphate


# ------------------------------------------------------------------------------
# Common Utilities
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


# ------------------------------------------------------------------------------
# Donor-based Summation (Donor Only)
# ------------------------------------------------------------------------------
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


# ------------------------------------------------------------------------------
# Summaries: Occurrence Bars + Existence Heatmap
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
# Build/Print Table: Donor Phosphate => top acceptor targets
# ------------------------------------------------------------------------------
def build_donor_acceptor_summary(species_id, bond_df, hbonds_per_index, time_per_frame_ns, donor_map, acceptor_map):
    """
    For each bond => map donor->pgroup, acceptor->pgroup => accumulate times
    """
    result={}
    for i, row in bond_df.iterrows():
        b_idx=row['idx']
        d_atom=row['donor']
        a_atom=row['acceptor']
        frames=hbonds_per_index[b_idx]
        engaged_ns=frames*time_per_frame_ns

        donor_pgroup=donor_map.get(d_atom,None)
        acceptor_pgroup=acceptor_map.get(a_atom,None)
        if donor_pgroup is None or acceptor_pgroup is None:
            continue

        if donor_pgroup not in result:
            result[donor_pgroup]={"sum_donor_time":0.0,"targets":{}}
        result[donor_pgroup]["sum_donor_time"]+=engaged_ns

        key_t=(acceptor_pgroup,a_atom)
        if key_t not in result[donor_pgroup]["targets"]:
            result[donor_pgroup]["targets"][key_t]=0.0
        result[donor_pgroup]["targets"][key_t]+=engaged_ns
    return result


def print_donor_acceptor_table(species_id, data):
    """
    Print lines for each donor pgroup => top acceptor targets
    """
    print(f"Species: {species_id}\n")

    donor_sorted=sorted(data.keys())  # e.g. P1..P6
    for dpgroup in donor_sorted:
        sum_time=data[dpgroup]["sum_donor_time"]
        print(f"  Donor {dpgroup} => total time: {sum_time:.2f} ns")
        # top 3 targets
        t_map=data[dpgroup]["targets"]
        t_list=sorted(t_map.items(), key=lambda x:x[1], reverse=True)
        top_names=["main target Pi","second target Pi","third target Pi"]
        for i, ((acc_pg,acc_atom),val_ns) in enumerate(t_list):
            if i<3:
                print(f"    {top_names[i]:20s} => {acc_pg}({acc_atom}) = {val_ns/sum_time*100:.0f} %")
        print()
    print("-"*70)
    print()


# ------------------------------------------------------------------------------
# 2×2 Figure Layout
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# The 2×2 Figure Layout
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# New Plot: Donor → Phosphate Group Distribution (Violin)
# ------------------------------------------------------------------------------
#def plot_phosphate_group_violin(...):
#    pass  # keep your code or adapt

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
def main():
    logger.info("Starting combined analysis script with textual table output...")

    script_dir=os.path.dirname(os.path.abspath(__file__))
    try:
        project_path=find_project_root(script_dir)
    except FileNotFoundError as e:
        logger.error(str(e))
        sys.exit(1)

    process_dir=os.path.join(project_path,"process.nobackup")
    if not os.path.isdir(process_dir):
        logger.error(f"Invalid process dir => {process_dir}")
        sys.exit(1)
    logger.info(f"Using process directory => {process_dir}")

    time_per_frame_ns=0.01
    process_path=Path(process_dir)

    # We'll do the same steps, but define bond_df=None at start of each loop
    donor_acceptor_summaries={}
    # gather everything else => ...
    folder_list=[]
    all_species_ids=[]

    for fname in sorted(os.listdir(process_path)):
        if not fname.startswith("IP_"):
            continue
        folder_path=process_path/fname
        if not folder_path.is_dir():
            continue

        species_id=fname.replace("IP_","")
        n_ones=species_id.count('1')
        if n_ones not in [3,4,5]:
            continue

        # add to folder_list
        folder_list.append((fname,species_id,n_ones))
        all_species_ids.append(species_id)

    # sort
    folder_list.sort(key=lambda x:(-x[2], x[1]))

    # define placeholders
    hbonds_data=[]
    lifetimes_data=[]
    distance_data=[]
    engaged_data=[]
    df_phosphate_rows=[]
    donor_atoms_by_species={}
    occurrence_data={}
    existence_data={}

    for folder_name,species_id,n_ones in folder_list:
        fpath=process_path/folder_name

        hb_num_file=fpath/"analyze_final_sim"/"h_bonds"/"intra_hb_num.xvg"
        xpm_file   =fpath/"analyze_final_sim"/"h_bonds"/"intra_hb_matrix.xpm"
        dist_file  =fpath/"analyze_final_sim"/"h_bonds"/"intra_hb_dist.xvg"
        log_file   =fpath/"analyze_final_sim"/"h_bonds"/"intra_hb.log"

        # define bond_df=None => so if we skip, it's still defined
        bond_df=None

        # check files
        if not (hb_num_file.is_file() and xpm_file.is_file() and dist_file.is_file() and log_file.is_file()):
            logger.warning(f"Skipping {species_id} => missing hbond files.")
            continue

        # load hbond counts => if empty => skip
        data_num=load_hb_num_xvg(str(hb_num_file))
        if data_num.size<2:
            logger.warning(f"No hbond count data => {species_id}")
            continue
        hbond_counts=data_num[:,1]
        for val in hbond_counts:
            hbonds_data.append((species_id,n_ones,val))

        # parse XPM => analyze => lifetimes
        xpm_mat,meta=parse_xpm(str(xpm_file))
        analysis_res=analyze_hydrogen_bonds(xpm_mat, meta)

        lf_frames=[lf for bond_lf in analysis_res['lifetimes'] for lf in bond_lf]
        lf_ns=[f*time_per_frame_ns for f in lf_frames]
        for val in lf_ns:
            lifetimes_data.append((species_id,n_ones,val))

        # distances
        d_list=[]
        with open(dist_file,'r') as df_:
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
        arr_d=np.array(d_list)
        if arr_d.size<2:
            logger.warning(f"No distance data => {species_id}")
            continue
        arr_d=arr_d[arr_d[:,0]>=0.2]
        for row in arr_d:
            distance_data.append((species_id,n_ones,row[0],row[1]))

        # parse .log => bond_df
        bond_df=parse_hbond_log_to_dataframe(str(log_file))
        if bond_df.empty:
            logger.warning(f"No hbond pairs => {species_id}")
            continue
        # now we have bond_df, proceed
        donor_atoms_by_species[species_id]=set(bond_df['donor'].unique())

        # frames per bond
        hpi=analysis_res['hbonds_per_index']
        donor_counts=count_donor_occurrences_from_matrix(bond_df,hpi)
        total_donor_frames=donor_counts.sum()
        engaged_ns=(total_donor_frames*time_per_frame_ns)/n_ones
        engaged_data.append((species_id,n_ones,engaged_ns))

        # build phosphate data
        for donor_atom,frames in donor_counts.items():
            engaged_time_ns=frames*time_per_frame_ns
            phosphate_label=donor_to_phosphate.get(donor_atom,None)
            if phosphate_label is not None:
                df_phosphate_rows.append({
                    'species':species_id,
                    'num_ones':n_ones,
                    'donor_atom':donor_atom,
                    'phosphate_group':phosphate_label,
                    'engaged_ns':engaged_time_ns
                })

        # build donor-acceptor summary
        summary_data=build_donor_acceptor_summary(
            species_id,bond_df,hpi,time_per_frame_ns,
            donor_to_phosphate,acceptor_to_phosphate
        )
        # store
        donor_acceptor_summaries[species_id]=summary_data

    # now build final DataFrames
    df_hbonds=pd.DataFrame(hbonds_data,columns=["species","num_ones","hbond_count"])
    df_lifetimes=pd.DataFrame(lifetimes_data,columns=["species","num_ones","lifetime_ns"])
    df_distance=pd.DataFrame(distance_data,columns=["species","num_ones","distance_nm","frequency"])
    df_engaged=pd.DataFrame(engaged_data,columns=["species","num_ones","engaged_time_ns"])

    # build color dict
    all_species_ids=sorted(set(all_species_ids))
    palette=sns.color_palette("tab10",len(all_species_ids))
    color_dict={}
    for i,sp in enumerate(all_species_ids):
        color_dict[sp]=palette[i%len(palette)]

    # optionally do 2x2 figure, summary plot, phosphate violin
    # ...
    # at the end, print textual table
    for sp, data in donor_acceptor_summaries.items():
        print_donor_acceptor_table(sp, data)

    logger.info("All analysis steps completed successfully!")

if __name__=="__main__":
    main()