import os
import re
import numpy as np
import matplotlib.pyplot as plt
import sys

import sys
import os
import re
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from itertools import combinations

sns.set_context("talk")


def load_hb_num_xvg(filename):
    """
    Load data from a .xvg file, skipping lines that start with '@' or '#'.

    Parameters
    ----------
    filename : str
        Path to the .xvg file.

    Returns
    -------
    np.ndarray
        A 2D array with the data parsed from the file.
    """
    data = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith(('@', '#')):
                # Skip header/comment lines
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    floats = [float(part) for part in parts]
                    data.append(floats)
                except ValueError:
                    # Skip lines that don't contain valid floats
                    continue
    return np.array(data)


def refine_atom_name(atom):
    """
    Refines the atom name based on specified rules:
      - If the atom is named 'O', change it to 'O1'.
      - If the atom is named 'O<number>', increment the number by 1 
        (e.g., 'O10' -> 'O11').

    Parameters
    ----------
    atom : str
        Original atom name.

    Returns
    -------
    str
        Refined atom name (with incremented numbering if applicable).
    """
    # Match atoms like 'O' (no number)
    if re.fullmatch(r'[A-Za-z]+', atom):
        return atom + '1'

    # Match atoms like 'O10', 'O2', etc.
    match = re.fullmatch(r'([A-Za-z]+)(\d+)', atom)
    if match:
        name = match.group(1)
        number = int(match.group(2)) + 1  # Increment the number by 1
        return f"{name}{number}"

    # If atom name doesn't match expected patterns, return as is
    return atom


def parse_xpm(file_path):
    """
    Parses an XPM file and converts it to a binary NumPy array.

    Assumes that '#FF0000' in the file indicates a "present bond" (1),
    and all other or unknown colors map to 0.

    Parameters
    ----------
    file_path : str
        Path to the XPM file.

    Returns
    -------
    tuple of (np.ndarray, dict)
        A (data_matrix, metadata) pair, where data_matrix is a 2D array 
        (binary, 0 or 1) representing hydrogen bonds, and metadata is a 
        dictionary with possible keys like 'title'.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize variables
    metadata = {}
    data_lines = []
    color_map = {}
    header_found = False
    data_started = False

    width = None
    height = None
    num_colors = None
    chars_per_pixel = None

    for line in lines:
        line = line.strip()

        # Extract metadata
        if line.startswith("/*") and not data_started:
            comment = line.strip("/* ").strip(" */")
            if ':' in comment:
                key, value = comment.split(":", 1)
                metadata[key.strip().lower()] = value.strip().strip('"')
            continue

        # Identify the start of data
        if line.startswith('static char'):
            continue  # Skip static declaration line

        if line.startswith('"') and not header_found:
            # This line should contain width, height, num_colors, chars_per_pixel
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
            # Read color definitions
            color_def = line.strip('",')
            # e.g. "   c #FFFFFF " or "o  c #FF0000 "
            match = re.match(
                r'(.{%d})\s+c\s+(\S+)' % chars_per_pixel,
                color_def
            )
            if match:
                symbol = match.group(1)
                color = match.group(2)
                color_map[symbol] = color
            if len(color_map) == num_colors:
                data_started = True
            continue

        if data_started:
            # Read data lines for the matrix
            if line.startswith('"'):
                data_line = line.strip('",')
                data_lines.append(data_line)

    data_matrix = np.zeros((height, width), dtype=int)

    for y, line in enumerate(data_lines):
        for x, char in enumerate(line):
            if char in color_map:
                if color_map[char] == '#FF0000':  # Present bond
                    data_matrix[y, x] = 1
                else:
                    data_matrix[y, x] = 0
            else:
                data_matrix[y, x] = 0  # Default to 0 if unknown symbol

    return data_matrix, metadata


def analyze_hydrogen_bonds(data_matrix, metadata):
    """
    Analyzes the hydrogen bond data to provide some summary metrics,
    such as hydrogen bonds over time and lifetimes.

    Parameters
    ----------
    data_matrix : np.ndarray
        A binary matrix (rows = hydrogen bonds, columns = frames).
    metadata : dict
        Additional info from the XPM file, e.g., title, etc.

    Returns
    -------
    dict
        A dictionary containing:
            - 'hbonds_over_time': 1D array with sum of bonds over frames.
            - 'hbonds_per_index': 1D array with sum of bonds for each bond index.
            - 'lifetimes': A list of lists, where each sublist contains the 
                           continuous "on" durations (in frames) for each bond.
    """
    analysis_results = {}

    # Total number of hydrogen bonds over time (summing rows for each column)
    hbonds_over_time = np.sum(data_matrix, axis=0)

    # Total occurrence of each hydrogen bond (summing columns for each row)
    hbonds_per_index = np.sum(data_matrix, axis=1)

    # Lifetime distribution requires tracking continuous sequences of 1s
    lifetimes = []
    for bond in data_matrix:
        current_lifetime = 0
        bond_lifetimes = []
        for state in bond:
            if state == 1:
                current_lifetime += 1
            else:
                if current_lifetime > 0:
                    bond_lifetimes.append(current_lifetime)
                    current_lifetime = 0
        # If the final frames ended with 1s:
        if current_lifetime > 0:
            bond_lifetimes.append(current_lifetime)
        lifetimes.append(bond_lifetimes)

    analysis_results['hbonds_over_time'] = hbonds_over_time
    analysis_results['hbonds_per_index'] = hbonds_per_index
    analysis_results['lifetimes'] = lifetimes

    return analysis_results


def parse_hbond_log_to_dataframe(file_path):
    """
    Parses a GROMACS .log file to extract donor-acceptor pairs with 
    refined atom names.

    Parameters
    ----------
    file_path : str
        Path to the .log file.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns: [idx, donor, acceptor].
        'idx' is a zero-based index for each H-bond pair.
    """
    hbond_pairs = []

    line_pattern = re.compile(r'^\s*(\S+)\s+-\s+(\S+)\s*$')
    atom_pattern = re.compile(r'^[A-Za-z]+\d+([A-Za-z]+\d*)$')

    with open(file_path, 'r') as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('"""') or line.startswith('*'):
                continue

            match = line_pattern.match(line)
            if match:
                donor_full = match.group(1)
                acceptor_full = match.group(2)

                donor_match = atom_pattern.match(donor_full)
                acceptor_match = atom_pattern.match(acceptor_full)

                if donor_match and acceptor_match:
                    donor_atom = donor_match.group(1)
                    acceptor_atom = acceptor_match.group(1)

                    refined_donor = refine_atom_name(donor_atom)
                    refined_acceptor = refine_atom_name(acceptor_atom)

                    hbond_pairs.append({
                        'donor': refined_donor,
                        'acceptor': refined_acceptor
                    })
                else:
                    print(f"Warning (Line {line_number}): Could not parse atoms in line: {line}")
            else:
                print(f"Warning (Line {line_number}): "
                      f"Line did not match expected format: {line}")

    df = pd.DataFrame(hbond_pairs)
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'idx'}, inplace=True)
    df['idx'] = df.index

    return df


def find_project_root(current_dir, project_name="electrofit"):
    """
    Find the root directory of the project by looking for the given 
    folder name (`project_name`) upward from `current_dir`.

    Parameters
    ----------
    current_dir : str
        Starting directory path.
    project_name : str
        Name of the project folder to look for.

    Returns
    -------
    str
        Path to the outermost matched project directory.

    Raises
    ------
    FileNotFoundError
        If the project_name directory is not found up the chain.
    """
    root = None
    while True:
        parent_dir = os.path.dirname(current_dir)
        if os.path.basename(current_dir) == project_name:
            root = current_dir
        if parent_dir == current_dir:
            # Reached filesystem root
            if root is None:
                raise FileNotFoundError(f"Project root directory '{project_name}' not found.")
            return root
        current_dir = parent_dir


# --------first naive approach 

def generate_intra_hbond_summary(
        process_dir, 
        output_file="intra_hbonds_summary.pdf", 
        time_per_frame_ns=0.01
    ):
    """
    Generates a single summary figure (3 columns x N rows) of intramolecular
    hydrogen-bond data for all species found under `process_dir`.
    
    Columns:
      (1) Violin plot (horizontal) of the distribution of # of H-bonds over time,
          with a dashed line showing the mean.
      (2) Lifetime histogram (in ns).
      (3) D–A distance distribution plot.
    
    Each row corresponds to one folder/species. The species identifier is
    derived by stripping 'IP_' from the folder name (e.g., 'IP_101010' -> '101010').

    Parameters
    ----------
    process_dir : str
        Path to the main directory containing subfolders for each species.
    output_file : str
        Filename (PDF) to save the final summary figure.
    time_per_frame_ns : float
        Duration of each frame in nanoseconds. Default is 0.01 ns.
    """

    sns.set_context("talk")  # or use your preferred style/theme

    # This list will store data from each species folder
    species_data_list = []

    # Loop over subdirectories looking for something like "IP_XXXXXX"
    for folder_name in sorted(os.listdir(process_dir)):
        folder_path = os.path.join(process_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue

        # Skip folders that do not start with "IP_"
        if not folder_name.startswith("IP_"):
            continue

        # Extract a simpler identifier by stripping "IP_"
        identifier = folder_name.replace("IP_", "")

        # Path to the directory where H-bond files are
        analyze_dir = os.path.join(folder_path, "analyze_final_sim", "h_bonds")

        # The intramolecular H-bond files needed
        intra_num_xvg = os.path.join(analyze_dir, "intra_hb_num.xvg")
        intra_matrix_xpm = os.path.join(analyze_dir, "intra_hb_matrix.xpm")
        intra_dist_xvg = os.path.join(analyze_dir, "intra_hb_dist.xvg")

        # Check file existence
        if (not os.path.isfile(intra_num_xvg) or
            not os.path.isfile(intra_matrix_xpm) or
            not os.path.isfile(intra_dist_xvg)):
            print(f"Warning: Missing required files for {folder_name}. Skipped.")
            continue

        # -------------------------------------------------------
        # (1) Load # of H-bonds vs. time from intra_hb_num.xvg
        # -------------------------------------------------------
        hbond_time_data = load_hb_num_xvg(intra_num_xvg)
        if hbond_time_data.size < 2:
            print(f"Warning: No valid data in {intra_num_xvg} for {folder_name}. Skipped.")
            continue

        # Typically:
        #   column 0: time in ps
        #   column 1: # of H-bonds
        #   column 2: # of pairs within 0.35 nm (sometimes)
        # We'll use column 1 for the distribution
        hbond_count_dist = hbond_time_data[:, 1]

        # -------------------------------------------------------
        # (2) Load matrix to compute lifetime from .xpm
        # -------------------------------------------------------
        data_matrix, meta = parse_xpm(intra_matrix_xpm)
        analysis = analyze_hydrogen_bonds(data_matrix, meta)

        # Flatten lifetime data in frames
        all_frame_lifetimes = [
            lifetime
            for bond_lifetimes in analysis['lifetimes']
            for lifetime in bond_lifetimes
        ]
        # Convert frames -> ns
        lifetime_ns = [lf * time_per_frame_ns for lf in all_frame_lifetimes]

        # -------------------------------------------------------
        # (3) Distance distribution from intra_hb_dist.xvg
        # -------------------------------------------------------
        dist_data = []
        with open(intra_dist_xvg, 'r') as f:
            for line in f:
                line = line.strip()
                if (not line) or line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) == 2:
                    try:
                        xval = float(parts[0])
                        yval = float(parts[1])
                        dist_data.append((xval, yval))
                    except ValueError:
                        pass
        dist_data = np.array(dist_data)
        if dist_data.size < 2:
            print(f"Warning: No distance data in {intra_dist_xvg} for {folder_name}. Skipped.")
            continue

        dist_x = dist_data[:, 0]
        dist_y = dist_data[:, 1]

        # Collect all data for this species
        species_data_list.append({
            "identifier": identifier,
            "hbond_count_dist": hbond_count_dist,
            "lifetime_ns": lifetime_ns,
            "dist_x": dist_x,
            "dist_y": dist_y
        })

    if not species_data_list:
        print("No data found for any species. Exiting without plotting.")
        return

    # Sort or keep them in the order discovered:
    # species_data_list is already sorted by folder_name due to sorted()

    # We'll create a figure with N rows, 3 columns
    n_species = len(species_data_list)
    fig, axes = plt.subplots(nrows=n_species, ncols=3, figsize=(18, 4 * n_species))
    if n_species == 1:
        # Ensure 2D indexing if there's only one species
        axes = np.array([axes])

    # Create each row of plots
    for i, sp_data in enumerate(species_data_list):
        identifier = sp_data['identifier']
        hbond_count_dist = sp_data['hbond_count_dist']
        lifetime_ns = sp_data['lifetime_ns']
        dist_x = sp_data['dist_x']
        dist_y = sp_data['dist_y']

        # ===== (1) Violin Plot of # of H-bonds (horizontal) =====
        ax_violin = axes[i, 0]
        # Seaborn usage: pass an array to x=... with orient='h'
        sns.violinplot(
            x=hbond_count_dist,
            orient='h',
            color='lightblue',
            cut=0,
            ax=ax_violin
        )
        # Mark the mean with a dashed vertical line
        mean_val = np.mean(hbond_count_dist)
        ax_violin.axvline(mean_val, color='black', linestyle='--', label=f"Mean={mean_val:.1f}")
        ax_violin.legend()

        # Aesthetics
        ax_violin.set_title(f"{identifier}\n# of H-bonds")
        ax_violin.set_xlabel("H-bond Count")
        ax_violin.set_ylabel("")  # removing y-axis label, since it's a single distribution

        # ===== (2) Lifetime Histogram =====
        ax_life = axes[i, 1]
        if len(lifetime_ns) > 0:
            ax_life.hist(lifetime_ns, bins=30, color='darkred', alpha=0.7)
        ax_life.set_title(f"{identifier}\nH-bond Lifetime")
        ax_life.set_xlabel("Lifetime (ns)")
        ax_life.set_ylabel("Frequency")

        # ===== (3) Distance Distribution =====
        ax_dist = axes[i, 2]
        ax_dist.plot(dist_x, dist_y, color='darkblue')
        ax_dist.set_title(f"{identifier}\nD–A Distance")
        ax_dist.set_xlabel("Distance (nm)")
        ax_dist.set_ylabel("Frequency")

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"[DONE] Summary figure saved as '{output_file}'")

# --------- three plots: arbitrary ordering

def create_three_plots_for_intra_hbonds(process_dir,
                                        output_prefix="intra_hbonds",
                                        time_per_frame_ns=0.01):
    """
    Creates three separate plots (each in its own figure) for all species
    under 'process_dir', focusing on intramolecular H-bonds:
    
    1) Violin plot of H-bond counts distribution (x=species, y=H-bond count).
    2) Lifetime distribution, overlaid or side-by-side, for all species.
    3) D–A distance distribution, overlaid for all species.

    Files expected in each species folder (e.g., 'IP_101010/analyze_final_sim/h_bonds'):
      - intra_hb_num.xvg
      - intra_hb_matrix.xpm
      - intra_hb_dist.xvg

    Parameters
    ----------
    process_dir : str
        Path to the directory containing species folders named like 'IP_XXXXXX'.
    output_prefix : str
        Prefix for the saved figures (e.g. 'intra_hbonds_violin.pdf', etc).
    time_per_frame_ns : float
        Duration of each frame in nanoseconds (used to convert lifetimes).
    """
    sns.set_context("talk")  # or style of your choice

    # Prepare empty lists to gather data from all species
    hbonds_data = []     # For violin plot: each row => (species, hbond_count)
    lifetimes_data = []  # For lifetime distribution => (species, lifetime_ns)
    distance_data = []   # For D-A distribution => (species, distance_nm, frequency)

    # ------------------------------------------------------------------
    # 1) Discover species folders
    # ------------------------------------------------------------------
    from pathlib import Path
    process_path = Path(process_dir)
    if not process_path.is_dir():
        raise ValueError(f"'{process_dir}' is not a valid directory.")


    # Loop through subdirectories
    for folder_name in sorted(os.listdir(process_path)):
        if not folder_name.startswith("IP_"):
            continue  # skip non-species folders
        folder_path = process_path / folder_name
        if not folder_path.is_dir():
            continue

        # Extract a simpler identifier (strip off 'IP_')
        species_id = folder_name.replace("IP_", "")

        # The intramolecular H-bond files:
        analyze_dir = folder_path / "analyze_final_sim" / "h_bonds"
        num_file = analyze_dir / "intra_hb_num.xvg"
        matrix_file = analyze_dir / "intra_hb_matrix.xpm"
        dist_file = analyze_dir / "intra_hb_dist.xvg"

        if not (num_file.is_file() and matrix_file.is_file() and dist_file.is_file()):
            print(f"Warning: Missing intramolecular H-bond files in {analyze_dir}.")
            continue

        # ------------------------------------------------------
        # Load # of H-bonds vs. time
        # ------------------------------------------------------
        data_num = load_hb_num_xvg(str(num_file))
        if data_num.size < 2:
            print(f"Warning: No valid data in {num_file}. Skipping species {species_id}.")
            continue
        
        # 2nd column => # of H-bonds at each frame
        hbond_counts = data_num[:, 1]

        # Append to the "hbonds_data" list
        for count_val in hbond_counts:
            hbonds_data.append((species_id, count_val))

        # ------------------------------------------------------
        # Lifetime data from XPM
        # ------------------------------------------------------
        xpm_matrix, meta = parse_xpm(str(matrix_file))
        analysis = analyze_hydrogen_bonds(xpm_matrix, meta)

        # Flatten all lifetimes (in frames)
        all_lifetimes_frames = [
            x for bond_lifetimes in analysis['lifetimes'] 
                for x in bond_lifetimes
        ]
        # Convert frames -> ns
        all_lifetimes_ns = [lf * time_per_frame_ns for lf in all_lifetimes_frames]

        # Store
        for lf_ns in all_lifetimes_ns:
            lifetimes_data.append((species_id, lf_ns))

        # ------------------------------------------------------
        # Distance distribution from .xvg
        # ------------------------------------------------------
        dist_list = []
        with open(dist_file, 'r') as fdist:
            for line in fdist:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) == 2:
                    try:
                        dx = float(parts[0])  # distance
                        dy = float(parts[1])  # frequency or count
                        dist_list.append((dx, dy))
                    except ValueError:
                        pass
        
        dist_list = np.array(dist_list)
        if dist_list.size < 2:
            print(f"Warning: Invalid distance data in {dist_file}. Skipping.")
            continue

        # Each row => (species, distance_nm, frequency)
        for row in dist_list:
            distance_data.append((species_id, row[0], row[1]))

    # ------------------------------------------------------------------
    # Convert data into Pandas DataFrames (long format)
    # ------------------------------------------------------------------
    df_hbonds = pd.DataFrame(hbonds_data, columns=["species", "hbond_count"])
    df_lifetimes = pd.DataFrame(lifetimes_data, columns=["species", "lifetime_ns"])
    df_distance = pd.DataFrame(distance_data, columns=["species", "distance_nm", "frequency"])

    # If no data was collected, we stop
    if df_hbonds.empty or df_lifetimes.empty or df_distance.empty:
        print("No data found for any species. No plots will be generated.")
        return

    # ------------------------------------------------------------------
    # PLOT #1: Violin plot (H-bond counts) for all species in one figure
    # ------------------------------------------------------------------

    plt.figure(figsize=(8, 6))
    sns.violinplot(
        data=df_hbonds,
        y="species",       # species on the Y-axis
        x="hbond_count",   # H-bond count on the X-axis
        palette="husl",
        cut=0,
        orient='h'         # ensure horizontal orientation (redundant with y=..., x=...)
    )
    # Optionally add a stripplot on top for more detail:
    sns.stripplot(
        data=df_hbonds,
        y="species",
        x="hbond_count",
        color='k',
        alpha=0.3,
        size=2
    )

    plt.title("Intra H-bond Count Distribution (Violin, Horizontal)")
    plt.xlabel("H-bond Count")
    plt.ylabel("Species")
    plt.tight_layout()
    violin_out = f"{output_prefix}_violin.pdf"
    plt.savefig(violin_out, dpi=300)
    plt.close()
    print(f"Saved violin plot to {violin_out}")


    # ------------------------------------------------------------------
    # PLOT #2: Lifetime distribution for all species (Overlaid)
    # ------------------------------------------------------------------
    plt.figure(figsize=(10, 6))
    # One approach: overlay with histplot (stacked or layered)
    # Another approach: KDE curves. Example with multiple='layer'
    # For clarity, let's do a KDE with a distinct color per species:
    unique_species = df_lifetimes["species"].unique()
    palette = sns.color_palette("hls", len(unique_species))

    for sp, color in zip(unique_species, palette):
        data_subset = df_lifetimes[df_lifetimes["species"] == sp]["lifetime_ns"]
        # Use a KDE curve to avoid too much overlap
        sns.kdeplot(data_subset, label=sp, color=color, fill=True, alpha=0.2)
    
    plt.title("Intra H-bond Lifetime Distribution (All Species)")
    plt.xlabel("Lifetime (ns)")
    plt.ylabel("Density")
    plt.legend(title="Species")
    plt.tight_layout()
    plt.xlim(left=0.1, right=0.2)
    plt.ylim(top=3)
    plt.savefig(f"{output_prefix}_lifetime.pdf", dpi=300)
    plt.close()
    print(f"Saved lifetime distribution plot to {output_prefix}_lifetime.pdf")

    # ------------------------------------------------------------------
    # PLOT #3: D–A Distance distribution, overlaid lines per species
    # ------------------------------------------------------------------
    plt.figure(figsize=(10, 6))
    # Each species presumably has "distance_nm" vs. "frequency"
    # We'll do a line plot with hue='species'
    # But we need to be sure each species has a curve. If it's binned data, 
    # each species might have ~some discrete points. We'll group by species.
    sns.lineplot(
        data=df_distance,
        x="distance_nm",
        y="frequency",
        hue="species",
        palette="tab10"
    )
    plt.title("Intra H-bond D–A Distance Distribution (All Species)")
    plt.xlabel("Distance (nm)")
    plt.ylabel("Frequency")
    plt.legend(title="Species")
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_distance.pdf", dpi=300)
    plt.close()
    print(f"Saved D–A distance plot to {output_prefix}_distance.pdf")

def create_three_plots_for_intra_hbonds_s(process_dir,
                                        output_prefix="intra_hbonds",
                                        time_per_frame_ns=0.01):
    """
    Creates three separate plots (each in its own figure) for all species
    under 'process_dir', focusing on intramolecular H-bonds:

    1) Inverted violin plot of H-bond counts distribution:
       - x = H-bond count, y = species

    2) Lifetime distribution (KDE), but only for lifetimes >= 0.2 ns.

    3) D–A distance distribution, overlaid line plots for each species,
       showing only distances >= 0.2 nm on the x-axis.

    Files expected in each species folder:
      - intra_hb_num.xvg
      - intra_hb_matrix.xpm
      - intra_hb_dist.xvg

    Parameters
    ----------
    process_dir : str
        Path to the directory containing species folders named like 'IP_XXXXXX'.
    output_prefix : str
        Prefix for the saved figures (PDF/PNG).
    time_per_frame_ns : float
        Duration of each frame in nanoseconds (used to convert lifetimes).
    """
    sns.set_context("talk")  # or style/theme of your choice

    # Prepare empty lists to gather data from all species
    hbonds_data = []     # For violin plot: (species, hbond_count)
    lifetimes_data = []  # For lifetime distribution: (species, lifetime_ns)
    distance_data = []   # For distance distribution: (species, distance_nm, frequency)

    from pathlib import Path
    process_path = Path(process_dir)
    if not process_path.is_dir():
        raise ValueError(f"'{process_dir}' is not a valid directory.")


    # Loop through subdirectories
    for folder_name in sorted(os.listdir(process_path)):
        if not folder_name.startswith("IP_"):
            continue  # skip non-species folders
        folder_path = process_path / folder_name
        if not folder_path.is_dir():
            continue

        # Extract a simpler identifier (strip off 'IP_')
        species_id = folder_name.replace("IP_", "")

        # The intramolecular H-bond files:
        analyze_dir = folder_path / "analyze_final_sim" / "h_bonds"
        num_file = analyze_dir / "intra_hb_num.xvg"
        matrix_file = analyze_dir / "intra_hb_matrix.xpm"
        dist_file = analyze_dir / "intra_hb_dist.xvg"

        if not (num_file.is_file() and matrix_file.is_file() and dist_file.is_file()):
            print(f"Warning: Missing intramolecular H-bond files in {analyze_dir}.")
            continue

        # ------------------------------------------------------
        # (1) Load # of H-bonds vs. time
        # ------------------------------------------------------
        data_num = load_hb_num_xvg(str(num_file))
        if data_num.size < 2:
            print(f"Warning: No valid data in {num_file}. Skipping {species_id}.")
            continue
        
        # 2nd column => # of H-bonds at each frame
        hbond_counts = data_num[:, 1]
        for count_val in hbond_counts:
            hbonds_data.append((species_id, count_val))

        # ------------------------------------------------------
        # (2) Lifetime data from XPM
        # ------------------------------------------------------
        xpm_matrix, meta = parse_xpm(str(matrix_file))
        analysis = analyze_hydrogen_bonds(xpm_matrix, meta)

        # Flatten lifetimes (in frames)
        lifetimes_frames = [
            lf for bond_lifetimes in analysis['lifetimes']
               for lf in bond_lifetimes
        ]
        # Convert frames -> ns
        lifetimes_ns = [f * time_per_frame_ns for f in lifetimes_frames]

        # We only keep lifetimes >= 0.2 ns
        lifetimes_ns = [lf for lf in lifetimes_ns if lf >= 0.2]

        for lf_ns in lifetimes_ns:
            lifetimes_data.append((species_id, lf_ns))

        # ------------------------------------------------------
        # (3) Distance distribution (distance, frequency)
        # ------------------------------------------------------
        dist_list = []
        with open(dist_file, 'r') as fdist:
            for line in fdist:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) == 2:
                    try:
                        dx = float(parts[0])  # distance (nm)
                        dy = float(parts[1])  # frequency or count
                        dist_list.append((dx, dy))
                    except ValueError:
                        pass
        
        dist_list = np.array(dist_list)
        if dist_list.size < 2:
            print(f"Warning: No valid distance data in {dist_file}. Skipping.")
            continue
        
        # Only consider distances >= 0.2 nm
        dist_list = dist_list[dist_list[:, 0] >= 0.2]

        for row in dist_list:
            distance_data.append((species_id, row[0], row[1]))

    # ------------------------------------------------------------------
    # Convert into DataFrames
    # ------------------------------------------------------------------
    df_hbonds = pd.DataFrame(hbonds_data, columns=["species", "hbond_count"])
    df_lifetimes = pd.DataFrame(lifetimes_data, columns=["species", "lifetime_ns"])
    df_distance = pd.DataFrame(distance_data, columns=["species", "distance_nm", "frequency"])

    if df_hbonds.empty and df_lifetimes.empty and df_distance.empty:
        print("No data found for any species. No plots will be generated.")
        return

    # ------------------------------------------------------------------
    # PLOT 1: Violin plot (H-bond counts), species on y-axis
    # ------------------------------------------------------------------
    plt.figure(figsize=(8, 6))
    sns.violinplot(
        data=df_hbonds,
        y="species",       # species on the Y-axis
        x="hbond_count",   # H-bond count on the X-axis
        palette="husl",
        cut=0,
        orient='h'         # ensure horizontal orientation (redundant with y=..., x=...)
    )
    # Optionally add a stripplot on top for more detail:
    sns.stripplot(
        data=df_hbonds,
        y="species",
        x="hbond_count",
        color='k',
        alpha=0.3,
        size=2
    )

    plt.title("Intra H-bond Count Distribution (Violin, Horizontal)")
    plt.xlabel("H-bond Count")
    plt.ylabel("Species")
    plt.tight_layout()
    violin_out = f"{output_prefix}_violin.pdf"
    plt.savefig(violin_out, dpi=300)
    plt.close()
    print(f"Saved violin plot to {violin_out}")

    # ------------------------------------------------------------------
    # PLOT 2: Lifetime distribution (KDE), only lifetimes >= 0.2 ns
    # ------------------------------------------------------------------
    if not df_lifetimes.empty:
        plt.figure(figsize=(8, 6))
        # Overlaid KDE with one curve per species
        species_unique = df_lifetimes["species"].unique()
        palette = sns.color_palette("hls", len(species_unique))

        for sp, col in zip(species_unique, palette):
            subset = df_lifetimes[df_lifetimes["species"] == sp]["lifetime_ns"]
            if len(subset) == 0:
                continue
            sns.kdeplot(subset, label=sp, color=col, fill=True, alpha=0.2)

        plt.title("Intra H-bond Lifetime Distribution (>= 0.2 ns)")
        plt.xlabel("Lifetime (ns)")
        plt.ylabel("Density")
        plt.legend(title="Species")
        plt.xlim(left=0.1, right=0.4)  # Force x-axis to start at 0.2 ns
        plt.tight_layout()
        lifetime_out = f"{output_prefix}_lifetime_kde.pdf"
        plt.savefig(lifetime_out, dpi=300)
        plt.close()
        print(f"Saved lifetime distribution to {lifetime_out}")
    else:
        print("No lifetime data (>=0.2 ns) found. Skipping lifetime plot.")

    # ------------------------------------------------------------------
    # PLOT 3: D–A distance distribution, overlaid lines
    #         Only distances >= 0.2 nm
    # ------------------------------------------------------------------
    if not df_distance.empty:
        plt.figure(figsize=(8, 6))
        sns.lineplot(
            data=df_distance,
            x="distance_nm",
            y="frequency",
            hue="species",
            palette="tab10"
        )
        plt.title("Intra H-bond D–A Distance Distribution (>= 0.2 nm)")
        plt.xlabel("Distance (nm)")
        plt.ylabel("Frequency")
        plt.xlim(left=0.2)  # Show only from 0.2 nm onwards
        plt.legend(title="Species")
        plt.tight_layout()
        dist_out = f"{output_prefix}_distance.pdf"
        plt.savefig(dist_out, dpi=300)
        plt.close()
        print(f"Saved D–A distance distribution to {dist_out}")
    else:
        print("No distance data (>=0.2 nm) found. Skipping distance plot.")

# ---------- three plots grouped 

def create_three_plots_for_intra_hbonds_grouped(
    process_dir,
    output_prefix="intra_hbonds_grouped",
    time_per_frame_ns=0.01
):
    """
    Creates three separate figures (violin, lifetime, distance), each grouped by 
    how many '1's a species ID contains (5, 4, or 3).

    FIGURE 1 (Violin plots of H-bond counts): 3 subplots side by side.
      Subplot 1 => species with 5 ones
      Subplot 2 => species with 4 ones
      Subplot 3 => species with 3 ones

    FIGURE 2 (Lifetime distributions, KDE): 3 subplots STACKED (vertical).
      Subplot 1 (top) => species with 5 ones
      Subplot 2 (middle) => species with 4 ones
      Subplot 3 (bottom) => species with 3 ones

    FIGURE 3 (Distance distributions): 3 subplots side by side.
      Subplot 1 => species with 5 ones
      Subplot 2 => species with 4 ones
      Subplot 3 => species with 3 ones
      Each line for a species in a given subplot is vertically offset by 2.5*i 
      to make the curves distinguishable.

    Only folders whose IDs have exactly 3, 4, or 5 ones are considered. Others 
    are skipped.

    Parameters
    ----------
    process_dir : str
        Path to the directory containing species folders named like 'IP_XXXXXX'.
    output_prefix : str
        Prefix for the saved figure files.
    time_per_frame_ns : float
        Duration of each frame in nanoseconds (used to convert lifetimes).
    """

    sns.set_context("talk")

    # Gather data from all species
    hbonds_data = []     # (species, num_ones, hbond_count)
    lifetimes_data = []  # (species, num_ones, lifetime_ns)
    distance_data = []   # (species, num_ones, distance_nm, frequency)

    from pathlib import Path
    process_path = Path(process_dir)
    if not process_path.is_dir():
        raise ValueError(f"'{process_dir}' is not a valid directory.")


    # Loop over subdirectories (species)
    for folder_name in sorted(os.listdir(process_path)):
        if not folder_name.startswith("IP_"):
            continue
        folder_path = process_path / folder_name
        if not folder_path.is_dir():
            continue

        # Extract ID by removing "IP_"
        species_id = folder_name.replace("IP_", "")
        # Count how many '1' bits
        num_ones = species_id.count('1')
        if num_ones not in [3, 4, 5]:
            continue  # skip if not 3,4,5

        # The intramolecular H-bond files:
        analyze_dir = folder_path / "analyze_final_sim" / "h_bonds"
        num_file = analyze_dir / "intra_hb_num.xvg"
        matrix_file = analyze_dir / "intra_hb_matrix.xpm"
        dist_file = analyze_dir / "intra_hb_dist.xvg"

        if not (num_file.is_file() and matrix_file.is_file() and dist_file.is_file()):
            print(f"Warning: Missing H-bond files for {species_id}. Skipping.")
            continue

        # -- (1) # of H-bonds vs time --
        data_num = load_hb_num_xvg(str(num_file))
        if data_num.size < 2:
            print(f"Warning: No valid data in {num_file} for {species_id}.")
            continue
        hbond_counts = data_num[:, 1]  # 2nd column
        for cval in hbond_counts:
            hbonds_data.append((species_id, num_ones, cval))

        # -- (2) Lifetime data --
        xpm_matrix, meta = parse_xpm(str(matrix_file))
        analysis = analyze_hydrogen_bonds(xpm_matrix, meta)
        lifetimes_frames = [
            lf for bond_lf in analysis['lifetimes']
               for lf in bond_lf
        ]
        lifetimes_ns = [f * time_per_frame_ns for f in lifetimes_frames]
        print(np.shape(lifetimes_ns))
        for lf_ns in lifetimes_ns:
            lifetimes_data.append((species_id, num_ones, lf_ns))

        # -- (3) Distance distribution (distance, frequency) --
        dist_list = []
        with open(dist_file, 'r') as fdist:
            for line in fdist:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) == 2:
                    try:
                        dx = float(parts[0])
                        dy = float(parts[1])
                        dist_list.append((dx, dy))
                    except ValueError:
                        pass

        dist_list = np.array(dist_list)
        if dist_list.size < 2:
            print(f"Warning: No valid distance data in {dist_file} for {species_id}.")
            continue


        for row in dist_list:
            distance_data.append((species_id, num_ones, row[0], row[1]))

    # Convert to DataFrames
    df_hbonds = pd.DataFrame(hbonds_data, columns=["species", "num_ones", "hbond_count"])
    df_lifetimes = pd.DataFrame(lifetimes_data, columns=["species", "num_ones", "lifetime_ns"])
    df_distance = pd.DataFrame(distance_data, columns=["species", "num_ones", "distance_nm", "frequency"])

    if df_hbonds.empty and df_lifetimes.empty and df_distance.empty:
        print("No data found for any species with 3,4,5 ones. Exiting.")
        return

    # Helper to pick subplot index based on # of ones
    # For the violin & distance (side by side), we want order 5 -> col=0, 4->1, 3->2
    def subplot_index_3col(n_ones):
        mapping = {5: 0, 4: 1, 3: 2}
        return mapping[n_ones]

    # Helper for lifetime subplots (stacked): 
    # 5 -> row=0, 4->1, 3->2
    def subplot_index_3row(n_ones):
        mapping = {5: 0, 4: 1, 3: 2}
        return mapping[n_ones]
    
    # Create global color mapping for all species
    all_species = []
    if not df_hbonds.empty:
        all_species.extend(df_hbonds['species'].unique())
    if not df_lifetimes.empty:
        all_species.extend(df_lifetimes['species'].unique())
    if not df_distance.empty:
        all_species.extend(df_distance['species'].unique())
    unique_species = sorted(list(set(all_species)))
    color_palette = sns.color_palette("husl", len(unique_species))
    color_dict = {species: color for species, color in zip(unique_species, color_palette)}

    # ========================================================
    # FIGURE 1: Violin plots (3 subplots side-by-side)
    # ========================================================
    fig1, axes1 = plt.subplots(1, 3, figsize=(20, 6), sharey=True)
    fig1.suptitle("Intra H-bond Count Distribution (Grouped by # of ones)")

    for n_ones in [5, 4, 3]:
        ax = axes1[subplot_index_3col(n_ones)]
        sub_df = df_hbonds[df_hbonds["num_ones"] == n_ones]

        if sub_df.empty:
            ax.set_title(f"No data for {n_ones} ones")
            ax.set_xlabel("H-bond Count")
            ax.set_ylabel("")
            continue

        # Modified violinplot section:
        species_order = sorted(sub_df['species'].unique())
        palette = [color_dict[sp] for sp in species_order]
        
        sns.violinplot(
            data=sub_df,
            y="species",
            x="hbond_count",
            palette=palette,
            order=species_order,
            cut=0,
            orient='h',
            ax=ax
        )
        # Optional strip plot
        sns.stripplot(
            data=sub_df,
            y="species",
            x="hbond_count",
            color='k',
            alpha=0.3,
            size=2,
            ax=ax
        )
        ax.set_title(f"{n_ones} ones\n({sub_df['species'].nunique()} species)")
        ax.set_xlabel("H-bond Count")
        ax.set_ylabel("Species")

    plt.tight_layout()
    violin_out = f"{output_prefix}_violin.pdf"
    plt.savefig(violin_out, dpi=300)
    plt.close()
    print(f"Saved grouped violin plot to {violin_out}")

    # ========================================================
    # FIGURE 2: Lifetime distribution (3 subplots stacked)
    # ========================================================
    if not df_lifetimes.empty:
        # Calculate global bins using all data
        all_lifetimes = df_lifetimes["lifetime_ns"]
        bins_global = np.histogram_bin_edges(all_lifetimes, bins=100)
        
        fig2, axes2 = plt.subplots(3, 1, figsize=(8, 18), sharex=True)
        fig2.suptitle("Intra H-bond Lifetime Distribution (Grouped, stacked)")

        for n_ones in [5, 4, 3]:
            row_idx = subplot_index_3row(n_ones)
            ax = axes2[row_idx]
            sub_df = df_lifetimes[df_lifetimes["num_ones"] == n_ones]
            
            if sub_df.empty:
                ax.set_title(f"No data for {n_ones} ones")
                ax.set_ylabel("Density")
                continue

            species_list = sub_df["species"].unique()
                
            # Modified lifetime plot section:
            for sp in species_list:
                data_subset = sub_df[sub_df["species"] == sp]["lifetime_ns"]
                if len(data_subset) == 0:
                    continue
                col = color_dict[sp]
                ax.hist(data_subset,
                        bins=bins_global,
                        density=True,
                        histtype='step',
                        linewidth=1.5,
                        alpha=0.7,
                        ec=col,
                        label=sp)


            ax.set_title(f"{n_ones} ones ({len(species_list)} species)")
            ax.set_xlabel("Lifetime (ns)")
            ax.set_ylabel("Density")
            ax.set_xlim(0, 1)
            ax.set_yscale("log")
            ax.legend(title="Species", fontsize='small', frameon=False)

        plt.tight_layout()
        lifetime_out = f"{output_prefix}_lifetime_steps.pdf"
        plt.savefig(lifetime_out, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved stacked lifetime distribution to {lifetime_out}")
    else:
        print("No lifetime data. Skipping lifetime plot.")

    # ========================================================
    # FIGURE 3: Distance distribution (3 subplots side-by-side),
    #           each line offset by 2.5*i
    # ========================================================
    if not df_distance.empty:
        # Single figure + single axes
        fig, ax = plt.subplots(figsize=(12, 8))
        fig.suptitle("Intra H-bond Distance Distribution (All in One Plot, Vertical Offsets)")
        
        # Order: 3 ones at the bottom, 4 in the middle, 5 at the top
        group_order = [3, 4, 5]
        
        # Base offsets to separate the groups well (adjust to suit your data)
        # e.g. group 3 starts at 0, group 4 starts at 15, group 5 starts at 30
        base_offset = {
            3: 0,
            4: 25,
            5: 50
        }
        
        for n_ones in group_order:
            # Filter for this group of species
            sub_df = df_distance[df_distance["num_ones"] == n_ones]
            if sub_df.empty:
                continue
            
            # Sort species so the offset is applied consistently
            species_list = sorted(sub_df["species"].unique())
            
            # Plot each species in this group with an additional offset (i * 2.5)
            for i, sp in enumerate(species_list):
                sp_data = sub_df[sub_df["species"] == sp].copy()
                if sp_data.empty:
                    continue
                
                # Sort data by distance
                sp_data = sp_data.sort_values(by="distance_nm")
                
                # Total offset = group offset + species index offset
                offset = base_offset[n_ones] + i * 5
                
                # Get the color for this species
                col = color_dict.get(sp, "black")  # fallback color if not in dict
                
                ax.plot(
                    sp_data["distance_nm"],
                    sp_data["frequency"] + offset,
                    color=col,
                    label=f"{sp} (offset={offset})"
                )
        
        ax.set_xlabel("Distance (nm)")
        ax.set_ylabel("Frequency + offset")
        ax.legend(fontsize='small')
        
        plt.tight_layout()
        dist_out = f"{output_prefix}_distance_offset_single_plot.pdf"
        plt.savefig(dist_out, dpi=300)
        plt.close()
        print(f"Saved distance distribution with offsets (single plot) to {dist_out}")
    else:
        print("No distance data available. Skipping distance plot.")


# -------------- 4 plots

def create_four_plots_for_intra_hbonds_grouped(
    process_dir,
    output_prefix="intra_hbonds_grouped",
    time_per_frame_ns=0.01
):
    """
    Creates four separate figures (violin, lifetime, distance with offsets, box plots),
    each grouped by how many '1's a species ID contains (5, 4, or 3).

    - FIGURE 1: Violin plots of H-bond counts
    - FIGURE 2: Lifetime distributions (KDE) stacked vertically
    - FIGURE 3: Distance distributions on one plot with vertical offsets
    - FIGURE 4: Box plots of D–A distances

    Only folders whose IDs have exactly 3, 4, or 5 ones are considered. Others are skipped.

    Parameters
    ----------
    process_dir : str
        Path to the directory containing species folders named like 'IP_XXXXXX'.
    output_prefix : str
        Prefix for the saved figure files.
    time_per_frame_ns : float
        Duration of each frame in nanoseconds (used to convert lifetimes).
    """

    sns.set_context("talk")
    
    # Prepare empty lists to gather data from all species
    hbonds_data = []     # (species, num_ones, hbond_count)
    lifetimes_data = []  # (species, num_ones, lifetime_ns)
    distance_data = []   # (species, num_ones, distance_nm, frequency)

    from pathlib import Path
    process_path = Path(process_dir)
    if not process_path.is_dir():
        raise ValueError(f"'{process_dir}' is not a valid directory.")

    # -- Define or import your helper functions --
    # Ensure these functions are defined in your script or imported from modules:
    #   load_hb_num_xvg(file_path)
    #   parse_xpm(file_path)
    #   analyze_hydrogen_bonds(data_matrix, metadata)

    # Example color dictionary (ensure this is defined or imported appropriately)
    # You might have a predefined color dictionary; otherwise, define one here.
    # For demonstration, we'll generate a color palette dynamically.
    # You can replace this with your own `color_dict` if you have specific colors per species.
    species_colors = sns.color_palette("tab10", 11)  # Assuming up to 11 species
    color_dict = {}  # {species_id: color}
    
    # Collect unique species identifiers first to assign colors
    unique_species = []
    for folder_name in sorted(os.listdir(process_path)):
        if not folder_name.startswith("IP_"):
            continue
        folder_path = process_path / folder_name
        if not folder_path.is_dir():
            continue

        species_id = folder_name.replace("IP_", "")
        num_ones = species_id.count('1')
        if num_ones not in [3, 4, 5]:
            continue
        unique_species.append(species_id)
    
    # Assign colors
    for i, sp in enumerate(sorted(unique_species)):
        color_dict[sp] = species_colors[i % len(species_colors)]

    # Loop through subdirectories (species)
    for folder_name in sorted(os.listdir(process_path)):
        if not folder_name.startswith("IP_"):
            continue
        folder_path = process_path / folder_name
        if not folder_path.is_dir():
            continue

        # Extract ID by removing "IP_"
        species_id = folder_name.replace("IP_", "")
        # Count how many '1' bits
        num_ones = species_id.count('1')
        if num_ones not in [3, 4, 5]:
            continue  # skip if not 3,4,5

        # The intramolecular H-bond files:
        analyze_dir = folder_path / "analyze_final_sim" / "h_bonds"
        num_file = analyze_dir / "intra_hb_num.xvg"
        matrix_file = analyze_dir / "intra_hb_matrix.xpm"
        dist_file = analyze_dir / "intra_hb_dist.xvg"

        if not (num_file.is_file() and matrix_file.is_file() and dist_file.is_file()):
            print(f"Warning: Missing H-bond files for {species_id}. Skipping.")
            continue

        # -- (1) Load # of H-bonds vs time --
        data_num = load_hb_num_xvg(str(num_file))
        if data_num.size < 2:
            print(f"Warning: No valid data in {num_file} for {species_id}.")
            continue
        hbond_counts = data_num[:, 1]  # 2nd column
        for cval in hbond_counts:
            hbonds_data.append((species_id, num_ones, cval))

        # -- (2) Lifetime data --
        xpm_matrix, meta = parse_xpm(str(matrix_file))
        analysis = analyze_hydrogen_bonds(xpm_matrix, meta)
        lifetimes_frames = [
            lf for bond_lf in analysis['lifetimes']
               for lf in bond_lf
        ]
        lifetimes_ns = [f * time_per_frame_ns for f in lifetimes_frames]
        for lf_ns in lifetimes_ns:
            lifetimes_data.append((species_id, num_ones, lf_ns))

        # -- (3) Distance distribution (distance, frequency) --
        dist_list = []
        with open(dist_file, 'r') as fdist:
            for line in fdist:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) == 2:
                    try:
                        dx = float(parts[0])
                        dy = float(parts[1])
                        dist_list.append((dx, dy))
                    except ValueError:
                        pass

        dist_list = np.array(dist_list)
        if dist_list.size < 2:
            print(f"Warning: No valid distance data in {dist_file} for {species_id}.")
            continue
        # Only consider distances >= 0.2 nm
        dist_list = dist_list[dist_list[:, 0] >= 0.2]

        for row in dist_list:
            distance_data.append((species_id, num_ones, row[0], row[1]))

    # ------------------------------------------------------------------
    # Convert into DataFrames
    # ------------------------------------------------------------------
    df_hbonds = pd.DataFrame(hbonds_data, columns=["species", "num_ones", "hbond_count"])
    df_lifetimes = pd.DataFrame(lifetimes_data, columns=["species", "num_ones", "lifetime_ns"])
    df_distance = pd.DataFrame(distance_data, columns=["species", "num_ones", "distance_nm", "frequency"])

    if df_hbonds.empty and df_lifetimes.empty and df_distance.empty:
        print("No data found for any species with 3,4,5 ones. Exiting.")
        return

    # Helper to pick subplot index based on # of ones
    # For the violin & distance (side by side), we want order 5 -> col=0, 4->1, 3->2
    def subplot_index_3col(n_ones):
        mapping = {5: 0, 4: 1, 3: 2}
        return mapping[n_ones]

    # Helper for lifetime subplots (stacked): 
    # 5 -> row=0, 4->1, 3->2
    def subplot_index_3row(n_ones):
        mapping = {5: 0, 4: 1, 3: 2}
        return mapping[n_ones]

    # ========================================================
    # FIGURE 1: Violin plots (3 subplots side-by-side)
    # ========================================================
    fig1, axes1 = plt.subplots(1, 3, figsize=(20, 6), sharey=True)
    fig1.suptitle("Intra H-bond Count Distribution (Grouped by # of ones)")

    for n_ones in [5, 4, 3]:
        ax = axes1[subplot_index_3col(n_ones)]
        sub_df = df_hbonds[df_hbonds["num_ones"] == n_ones]

        if sub_df.empty:
            ax.set_title(f"No data for {n_ones} ones")
            ax.set_xlabel("H-bond Count")
            ax.set_ylabel("")
            continue

        sns.violinplot(
            data=sub_df,
            y="species",         # species on the Y-axis
            x="hbond_count",     # H-bond count on X-axis
            palette="husl",
            cut=0,
            orient='h',
            ax=ax
        )
        # Optional strip plot
        sns.stripplot(
            data=sub_df,
            y="species",
            x="hbond_count",
            color='k',
            alpha=0.3,
            size=2,
            ax=ax
        )
        ax.set_title(f"{n_ones} ones\n({sub_df['species'].nunique()} species)")
        ax.set_xlabel("H-bond Count")
        ax.set_ylabel("Species")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust to accommodate suptitle
    violin_out = f"{output_prefix}_violin.pdf"
    plt.savefig(violin_out, dpi=300)
    plt.close()
    print(f"Saved grouped violin plot to {violin_out}")

    # ========================================================
    # FIGURE 2: Lifetime distribution (3 subplots stacked)
    # ========================================================
    if not df_lifetimes.empty:
        fig2, axes2 = plt.subplots(3, 1, figsize=(10, 18), sharex=True)
        fig2.suptitle("Intra H-bond Lifetime Distribution (Grouped by # of ones)")

        for n_ones in [5, 4, 3]:
            row_idx = subplot_index_3row(n_ones)
            ax = axes2[row_idx]
            sub_df = df_lifetimes[df_lifetimes["num_ones"] == n_ones]
            if sub_df.empty:
                ax.set_title(f"No data for {n_ones} ones")
                ax.set_xlabel("Lifetime (ns)")
                ax.set_ylabel("Density")
                continue

            species_list = sorted(sub_df["species"].unique())
            palette = sns.color_palette("hls", len(species_list))
            for sp, col in zip(species_list, palette):
                data_subset = sub_df[sub_df["species"] == sp]["lifetime_ns"]
                if len(data_subset) == 0:
                    continue
                sns.kdeplot(data_subset, label=sp, color=col, fill=True, alpha=0.2, ax=ax)

            ax.set_title(f"{n_ones} ones\n({len(species_list)} species)")
            ax.set_xlabel("Lifetime (ns)")
            ax.set_ylabel("Density")
            ax.legend(title="Species", fontsize='small')

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        lifetime_out = f"{output_prefix}_lifetime_kde.pdf"
        plt.savefig(lifetime_out, dpi=300)
        plt.close()
        print(f"Saved stacked lifetime distribution to {lifetime_out}")
    else:
        print("No lifetime data. Skipping lifetime plot.")

    # ========================================================
    # FIGURE 3: Distance distribution (single plot with offsets)
    # ========================================================
    if not df_distance.empty:
        # Single figure + single axes
        fig3, ax3 = plt.subplots(figsize=(12, 8))
        fig3.suptitle("Intra H-bond Distance Distribution (All in One Plot, Vertical Offsets)")

        # Order: 3 ones at the bottom, 4 in the middle, 5 at the top
        group_order = [3, 4, 5]

        # Base offsets to separate the groups well (adjust to suit your data)
        # e.g. group 3 starts at 0, group 4 starts at 15, group 5 starts at 30
        base_offset = {
            3: 0,
            4: 15,
            5: 30
        }

        for n_ones in group_order:
            # Filter for this group of species
            sub_df = df_distance[df_distance["num_ones"] == n_ones]
            if sub_df.empty:
                continue

            # Sort species so the offset is applied consistently
            species_list = sorted(sub_df["species"].unique())

            # Plot each species in this group with an additional offset (i * 2.5)
            for i, sp in enumerate(species_list):
                sp_data = sub_df[sub_df["species"] == sp].copy()
                if sp_data.empty:
                    continue

                # Sort data by distance
                sp_data = sp_data.sort_values(by="distance_nm")

                # Total offset = group offset + species index offset
                offset = base_offset[n_ones] + i * 2.5

                # Get the color for this species
                col = color_dict.get(sp, "black")  # fallback color if not in dict

                ax3.plot(
                    sp_data["distance_nm"],
                    sp_data["frequency"] + offset,
                    color=col,
                    label=f"{sp} (offset={offset})"
                )

        ax3.set_xlabel("Distance (nm)")
        ax3.set_ylabel("Frequency + offset")
        ax3.legend(fontsize='small', bbox_to_anchor=(1.05, 1), loc='upper left')  # Place legend outside
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to make space for legend
        dist_out = f"{output_prefix}_distance_offset_single_plot.pdf"
        plt.savefig(dist_out, dpi=300)
        plt.close()
        print(f"Saved distance distribution with offsets (single plot) to {dist_out}")
    else:
        print("No distance data available. Skipping distance plot.")

    # ========================================================
    # FIGURE 4: Box plots of D–A Distance Distributions
    # ========================================================
    if not df_distance.empty:
        # Reconstruct raw distance data by repeating distance_nm according to frequency
        # WARNING: This can be memory-intensive if 'frequency' is large
        # Consider downsampling if necessary

        # To handle large data, you might sample instead of repeating
        # Here, we'll proceed with repeating, but be cautious

        print("Reconstructing raw distance data for box plots...")
        # Create a new DataFrame where each distance_nm is repeated 'frequency' times
        # First, ensure 'frequency' is integer
        df_distance['frequency'] = df_distance['frequency'].astype(int)
        # Expand the DataFrame
        df_distance_expanded = df_distance.loc[df_distance.index.repeat(df_distance['frequency'])].reset_index(drop=True)

        # Now, create the box plot
        plt.figure(figsize=(10, 12))  # Tall figure to accommodate all species
        sns.boxplot(
            data=df_distance_expanded,
            y="species",
            x="distance_nm",
            orient='h',
            palette=color_dict,
            showfliers=False  # Optionally hide outliers for clarity
        )
        # Optional: Add swarm plot for individual data points
        sns.swarmplot(
            data=df_distance_expanded,
            y="species",
            x="distance_nm",
            color='k',
            alpha=0.3,
            size=2
        )

        plt.title("Intra H-bond D–A Distance Distribution (Box Plots)")
        plt.xlabel("Distance (nm)")
        plt.ylabel("Species")
        plt.tight_layout()
        boxplot_out = f"{output_prefix}_distance_boxplot.pdf"
        plt.savefig(boxplot_out, dpi=300)
        plt.close()
        print(f"Saved distance distribution box plots to {boxplot_out}")
    else:
        print("No distance data available. Skipping box plot.")


def create_four_plots_for_intra_hbonds_grouped2(
    process_dir,
    output_prefix="intra_hbonds_grouped",
    time_per_frame_ns=0.01
):
    """
    Creates four separate figures (violin, lifetime, distance with offsets, box plots),
    each grouped by how many '1's a species ID contains (5, 4, or 3).

    - FIGURE 1: Violin plots of H-bond counts
    - FIGURE 2: Lifetime distributions (KDE) stacked vertically
    - FIGURE 3: Distance distributions on one plot with vertical offsets
    - FIGURE 4: Box plots of D–A distances

    Only folders whose IDs have exactly 3, 4, or 5 ones are considered. Others are skipped.

    Parameters
    ----------
    process_dir : str
        Path to the directory containing species folders named like 'IP_XXXXXX'.
    output_prefix : str
        Prefix for the saved figure files.
    time_per_frame_ns : float
        Duration of each frame in nanoseconds (used to convert lifetimes).
    """

    sns.set_context("talk")
    
    # Prepare empty lists to gather data from all species
    hbonds_data = []     # (species, num_ones, hbond_count)
    lifetimes_data = []  # (species, num_ones, lifetime_ns)
    distance_data = []   # (species, num_ones, distance_nm, frequency)

    from pathlib import Path
    process_path = Path(process_dir)
    if not process_path.is_dir():
        raise ValueError(f"'{process_dir}' is not a valid directory.")

    # Example color dictionary (ensure this is defined or imported appropriately)
    # You might have a predefined color dictionary; otherwise, define one here.
    # For demonstration, we'll generate a color palette dynamically.
    # You can replace this with your own `color_dict` if you have specific colors per species.
    species_colors = sns.color_palette("tab10", 11)  # Assuming up to 11 species
    color_dict = {}  # {species_id: color}
    
    # Collect unique species identifiers first to assign colors
    unique_species = []
    for folder_name in sorted(os.listdir(process_path)):
        if not folder_name.startswith("IP_"):
            continue
        folder_path = process_path / folder_name
        if not folder_path.is_dir():
            continue

        species_id = folder_name.replace("IP_", "")
        num_ones = species_id.count('1')
        if num_ones not in [3, 4, 5]:
            continue
        unique_species.append(species_id)
    
    # Assign colors
    for i, sp in enumerate(sorted(unique_species)):
        color_dict[sp] = species_colors[i % len(species_colors)]

    # Loop through subdirectories (species)
    for folder_name in sorted(os.listdir(process_path)):
        if not folder_name.startswith("IP_"):
            continue
        folder_path = process_path / folder_name
        if not folder_path.is_dir():
            continue

        # Extract ID by removing "IP_"
        species_id = folder_name.replace("IP_", "")
        # Count how many '1' bits
        num_ones = species_id.count('1')
        if num_ones not in [3, 4, 5]:
            continue  # skip if not 3,4,5

        # The intramolecular H-bond files:
        analyze_dir = folder_path / "analyze_final_sim" / "h_bonds"
        num_file = analyze_dir / "intra_hb_num.xvg"
        matrix_file = analyze_dir / "intra_hb_matrix.xpm"
        dist_file = analyze_dir / "intra_hb_dist.xvg"

        if not (num_file.is_file() and matrix_file.is_file() and dist_file.is_file()):
            print(f"Warning: Missing H-bond files for {species_id}. Skipping.")
            continue

        # -- (1) Load # of H-bonds vs time --
        data_num = load_hb_num_xvg(str(num_file))
        if data_num.size < 2:
            print(f"Warning: No valid data in {num_file} for {species_id}.")
            continue
        hbond_counts = data_num[:, 1]  # 2nd column
        for cval in hbond_counts:
            hbonds_data.append((species_id, num_ones, cval))

        # -- (2) Lifetime data --
        xpm_matrix, meta = parse_xpm(str(matrix_file))
        analysis = analyze_hydrogen_bonds(xpm_matrix, meta)
        lifetimes_frames = [
            lf for bond_lf in analysis['lifetimes']
               for lf in bond_lf
        ]
        lifetimes_ns = [f * time_per_frame_ns for f in lifetimes_frames]
        for lf_ns in lifetimes_ns:
            lifetimes_data.append((species_id, num_ones, lf_ns))

        # -- (3) Distance distribution (distance, frequency) --
        dist_list = []
        with open(dist_file, 'r') as fdist:
            for line in fdist:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('@'):
                    continue
                parts = line.split()
                if len(parts) == 2:
                    try:
                        dx = float(parts[0])
                        dy = float(parts[1])
                        dist_list.append((dx, dy))
                    except ValueError:
                        pass

        dist_list = np.array(dist_list)
        if dist_list.size < 2:
            print(f"Warning: No valid distance data in {dist_file} for {species_id}.")
            continue
        # Only consider distances >= 0.2 nm
        dist_list = dist_list[dist_list[:, 0] >= 0.2]

        for row in dist_list:
            distance_data.append((species_id, num_ones, row[0], row[1]))

    # ------------------------------------------------------------------
    # Convert into DataFrames
    # ------------------------------------------------------------------
    df_hbonds = pd.DataFrame(hbonds_data, columns=["species", "num_ones", "hbond_count"])
    df_lifetimes = pd.DataFrame(lifetimes_data, columns=["species", "num_ones", "lifetime_ns"])
    df_distance = pd.DataFrame(distance_data, columns=["species", "num_ones", "distance_nm", "frequency"])

    if df_hbonds.empty and df_lifetimes.empty and df_distance.empty:
        print("No data found for any species with 3,4,5 ones. Exiting.")
        return

    # Helper to pick subplot index based on # of ones
    # For the violin & distance (side by side), we want order 5 -> col=0, 4->1, 3->2
    def subplot_index_3col(n_ones):
        mapping = {5: 0, 4: 1, 3: 2}
        return mapping[n_ones]

    # Helper for lifetime subplots (stacked): 
    # 5 -> row=0, 4->1, 3->2
    def subplot_index_3row(n_ones):
        mapping = {5: 0, 4: 1, 3: 2}
        return mapping[n_ones]

    # ========================================================
    # FIGURE 1: Violin plots (3 subplots side-by-side)
    # ========================================================
    fig1, axes1 = plt.subplots(1, 3, figsize=(20, 6), sharey=True)
    fig1.suptitle("Intra H-bond Count Distribution (Grouped by # of ones)")

    for n_ones in [5, 4, 3]:
        ax = axes1[subplot_index_3col(n_ones)]
        sub_df = df_hbonds[df_hbonds["num_ones"] == n_ones]

        if sub_df.empty:
            ax.set_title(f"No data for {n_ones} ones")
            ax.set_xlabel("H-bond Count")
            ax.set_ylabel("")
            continue

        sns.violinplot(
            data=sub_df,
            y="species",         # species on the Y-axis
            x="hbond_count",     # H-bond count on X-axis
            palette="husl",
            cut=0,
            orient='h',
            ax=ax
        )
        # Optional strip plot
        sns.stripplot(
            data=sub_df,
            y="species",
            x="hbond_count",
            color='k',
            alpha=0.3,
            size=2,
            ax=ax
        )
        ax.set_title(f"{n_ones} ones\n({sub_df['species'].nunique()} species)")
        ax.set_xlabel("H-bond Count")
        ax.set_ylabel("Species")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust to accommodate suptitle
    violin_out = f"{output_prefix}_violin.pdf"
    plt.savefig(violin_out, dpi=300)
    plt.close()
    print(f"Saved grouped violin plot to {violin_out}")

    # ========================================================
    # FIGURE 2: Lifetime distribution (3 subplots stacked)
    # ========================================================
    if not df_lifetimes.empty:
        fig2, axes2 = plt.subplots(3, 1, figsize=(10, 18), sharex=True)
        fig2.suptitle("Intra H-bond Lifetime Distribution (Grouped by # of ones)")

        for n_ones in [5, 4, 3]:
            row_idx = subplot_index_3row(n_ones)
            ax = axes2[row_idx]
            sub_df = df_lifetimes[df_lifetimes["num_ones"] == n_ones]
            if sub_df.empty:
                ax.set_title(f"No data for {n_ones} ones")
                ax.set_xlabel("Lifetime (ns)")
                ax.set_ylabel("Density")
                continue

            species_list = sorted(sub_df["species"].unique())
            palette = sns.color_palette("hls", len(species_list))
            for sp, col in zip(species_list, palette):
                data_subset = sub_df[sub_df["species"] == sp]["lifetime_ns"]
                if len(data_subset) == 0:
                    continue
                sns.kdeplot(data_subset, label=sp, color=col, fill=True, alpha=0.2, ax=ax)

            ax.set_title(f"{n_ones} ones\n({len(species_list)} species)")
            ax.set_xlabel("Lifetime (ns)")
            ax.set_ylabel("Density")
            ax.legend(title="Species", fontsize='small')

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        lifetime_out = f"{output_prefix}_lifetime_kde.pdf"
        plt.savefig(lifetime_out, dpi=300)
        plt.close()
        print(f"Saved stacked lifetime distribution to {lifetime_out}")
    else:
        print("No lifetime data. Skipping lifetime plot.")

    # ========================================================
    # FIGURE 3: Distance distribution (single plot with offsets)
    # ========================================================
    if not df_distance.empty:
        # Single figure + single axes
        fig3, ax3 = plt.subplots(figsize=(12, 8))
        fig3.suptitle("Intra H-bond Distance Distribution (All in One Plot, Vertical Offsets)")

        # Order: 3 ones at the bottom, 4 in the middle, 5 at the top
        group_order = [3, 4, 5]

        # Base offsets to separate the groups well (adjust to suit your data)
        # e.g. group 3 starts at 0, group 4 starts at 15, group 5 starts at 30
        base_offset = {
            3: 0,
            4: 15,
            5: 30
        }

        for n_ones in group_order:
            # Filter for this group of species
            sub_df = df_distance[df_distance["num_ones"] == n_ones]
            if sub_df.empty:
                continue

            # Sort species so the offset is applied consistently
            species_list = sorted(sub_df["species"].unique())

            # Plot each species in this group with an additional offset (i * 2.5)
            for i, sp in enumerate(species_list):
                sp_data = sub_df[sub_df["species"] == sp].copy()
                if sp_data.empty:
                    continue

                # Sort data by distance
                sp_data = sp_data.sort_values(by="distance_nm")

                # Total offset = group offset + species index offset
                offset = base_offset[n_ones] + i * 2.5

                # Get the color for this species
                col = color_dict.get(sp, "black")  # fallback color if not in dict

                ax3.plot(
                    sp_data["distance_nm"],
                    sp_data["frequency"] + offset,
                    color=col,
                    label=f"{sp} (offset={offset})"
                )

        ax3.set_xlabel("Distance (nm)")
        ax3.set_ylabel("Frequency + offset")
        ax3.legend(fontsize='small', bbox_to_anchor=(1.05, 1), loc='upper left')  # Place legend outside
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to make space for legend
        dist_out = f"{output_prefix}_distance_offset_single_plot.pdf"
        plt.savefig(dist_out, dpi=300)
        plt.close()
        print(f"Saved distance distribution with offsets (single plot) to {dist_out}")
    else:
        print("No distance data available. Skipping distance plot.")

    # ========================================================
    # FIGURE 4: Box plots of D–A Distance Distributions
    # ========================================================
    if not df_distance.empty:
        # Reconstruct raw distance data by repeating distance_nm according to frequency
        # WARNING: This can be memory-intensive if 'frequency' is large
        # Consider downsampling if necessary

        print("Reconstructing raw distance data for box plots...")
        # Create a new DataFrame where each distance_nm is repeated 'frequency' times
        # First, ensure 'frequency' is integer
        print(df_distance['frequency'])
        df_distance['frequency'] = df_distance['frequency'].astype(int)
        # Expand the DataFrame
        df_distance_expanded = df_distance.loc[df_distance.index.repeat(df_distance['frequency'])].reset_index(drop=True)

        # Define species order: 5 ones first, then 4, then 3
        # Within each group, sort species alphabetically
        species_order = []
        for n_ones in [5, 4, 3]:
            species_in_group = sorted(df_distance[df_distance["num_ones"] == n_ones]["species"].unique())
            species_order.extend(species_in_group)

        # Now, create the box plot with the specified order
        plt.figure(figsize=(10, 12))  # Tall figure to accommodate all species
        sns.boxplot(
            data=df_distance_expanded,
            y="species",
            x="distance_nm",
            orient='h',
            palette=color_dict,
            order=species_order,
            showfliers=False  # Optionally hide outliers for clarity
        )
        # Optional: Add swarm plot for individual data points
        sns.swarmplot(
            data=df_distance_expanded,
            y="species",
            x="distance_nm",
            color='k',
            alpha=0.3,
            size=2,
            order=species_order
        )

        plt.title("Intra H-bond D–A Distance Distribution (Box Plots)")
        plt.xlabel("Distance (nm)")
        plt.ylabel("Species")
        plt.tight_layout()
        boxplot_out = f"{output_prefix}_distance_boxplot.pdf"
        plt.savefig(boxplot_out, dpi=300)
        plt.close()
        print(f"Saved distance distribution box plots to {boxplot_out}")
    else:
        print("No distance data available. Skipping box plot.")

#--------------------
#--------------------
script_dir = os.path.dirname(os.path.abspath(__file__))
project_path = find_project_root(current_dir=script_dir)

sys.path.append(project_path)

from electrofit.helper.file_manipulation import find_file_with_extension
from electrofit.commands.run_commands import run_command

#--------------------
#--------------------

if __name__ == "__main__":
    process_directory = os.path.join(project_path, "dummy_process")
    create_four_plots_for_intra_hbonds_grouped2(process_dir=process_directory)