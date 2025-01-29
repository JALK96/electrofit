import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# Set Seaborn context
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

    width = height = num_colors = chars_per_pixel = None

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
                width, height, num_colors, chars_per_pixel = map(int, tokens[:4])
                header_found = True
            continue

        if header_found and not data_started:
            # Read color definitions
            color_def = line.strip('",')
            # e.g. "   c #FFFFFF " or "o  c #FF0000 "
            match = re.match(
                rf'(.{{{chars_per_pixel}}})\s+c\s+(\S+)',
                color_def
            )
            if match:
                symbol, color = match.groups()
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
            color = color_map.get(char, None)
            if color == '#FF0000':  # Present bond
                data_matrix[y, x] = 1
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


def parse_gromacs_log(file_path):
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
            if not line or line.startswith(('#', '"""', '*')):
                continue

            match = line_pattern.match(line)
            if match:
                donor_full, acceptor_full = match.groups()

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
                print(f"Warning (Line {line_number}): Line did not match expected format: {line}")

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


def one_row_three_columns_box_violin(
    df_lifetimes,    # columns: [species, num_ones, lifetime_ns]
    df_distance,     # columns: [species, num_ones, distance_nm, frequency]
    df_hbonds,       # columns: [species, num_ones, hbond_count]
    color_dict=None, # e.g. { "species_id": "#RRGGBB" }
    output_prefix="intra_hbonds",
    figure_size=(20, 8),
    spacer_prefix="~space_"   # prefix for dummy categories
):
    """
    Single figure with 1 row × 3 columns:
      (1) Lifetime (Box Plot)
      (2) Distance (Box Plot expanded by frequency)
      (3) H-bond Count (Violin Plot)
    We insert dummy rows to create spacing between groups of species
    (5 -> 4 -> 3). Only the left column shows y-tick labels. The 2nd and 3rd
    columns hide y-labels. We add curly braces on the left column to mark
    each group.

    Naive approach:
      - Convert 'frequency' to int, replicate rows for distance.
      - Insert dummy categories like "~space_5to4~" between groups.
    Then use curlyBrace(...) to draw brackets along the left side.

    Parameters
    ----------
    df_lifetimes : DataFrame
        Must have: [species, num_ones, lifetime_ns]
    df_distance : DataFrame
        Must have: [species, num_ones, distance_nm, frequency]
    df_hbonds : DataFrame
        Must have: [species, num_ones, hbond_count]
    color_dict : dict, optional
        Maps species -> color. If None, a default palette is used.
    output_prefix : str
        Filename prefix for the saved figure.
    figure_size : tuple
        Size of the final figure.
    spacer_prefix : str
        Prefix for dummy category labels inserted between group transitions.
    """

    # 1) Check if all DF are empty
    if (df_lifetimes is None or df_lifetimes.empty) and \
       (df_distance is None or df_distance.empty) and \
       (df_hbonds is None or df_hbonds.empty):
        print("All data frames empty. Nothing to plot.")
        return

    # 2) Collect all species from the 3 DFs
    species_in_any = set()
    for df in [df_lifetimes, df_distance, df_hbonds]:
        if df is not None and not df.empty:
            species_in_any.update(df["species"].unique())

    # 3) Determine num_ones for each species, keep only {3,4,5}
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
        if val in [3, 4, 5]:
            all_with_num.append((sp, val))
    if not all_with_num:
        print("No species with num_ones in [3,4,5]. Exiting.")
        return

    # 4) Sort so 5-> top, 4-> middle, 3-> bottom; alpha tiebreak
    all_with_num.sort(key=lambda x: (-x[1], x[0]))

    # Insert dummy categories between group transitions
    grouped_species_order = []
    prev_n_ones = None
    for i, (sp, n_ones) in enumerate(all_with_num):
        if i > 0 and n_ones != prev_n_ones:
            # Insert a dummy
            dummy_label = f"{spacer_prefix}{prev_n_ones}to{n_ones}"
            grouped_species_order.append((dummy_label, None))
        grouped_species_order.append((sp, n_ones))
        prev_n_ones = n_ones

    # Full y-axis order includes real species and dummy spacers
    full_order = [t[0] for t in grouped_species_order]

    # 5) Expand distance if needed
    if df_distance is not None and not df_distance.empty:
        df_distance = df_distance.copy()
        df_distance['frequency'] = df_distance['frequency'].astype(int)
        df_distance_expanded = df_distance.loc[
            df_distance.index.repeat(df_distance['frequency'])
        ].reset_index(drop=True)
    else:
        df_distance_expanded = pd.DataFrame()

    # 6) Add dummy rows to each DF for the spacers so Seaborn leaves them blank
    def inject_spacers(df, measure_col):
        """
        Insert rows with species= dummy label, measure=NaN for each spacer.
        Returns a new DF that includes them.
        """
        if df is None or df.empty:
            return pd.DataFrame(columns=["species", measure_col, "num_ones"])
        new_df = df.copy()
        # For each spacer in grouped_species_order with n_ones=None
        spacers = []
        for sp, val in grouped_species_order:
            if val is None:  # indicates dummy
                spacers.append({
                    "species": sp,
                    measure_col: np.nan,
                    "num_ones": -1
                })
        if not spacers:
            return new_df
        dummy_df = pd.DataFrame(spacers)
        combined = pd.concat([new_df, dummy_df], ignore_index=True)
        return combined

    df_life2 = inject_spacers(df_lifetimes, "lifetime_ns")
    df_dist2 = inject_spacers(df_distance_expanded, "distance_nm")
    df_hbond2 = inject_spacers(df_hbonds, "hbond_count")

    # Extend the color_dict to include spacer labels with a default color
    if color_dict is None:
        color_dict = {}
    else:
        color_dict = color_dict.copy()  # Avoid mutating the original dict

    # Define a default color for spacers
    spacer_color = "#D3D3D3"  # Light grey

    for sp, n in grouped_species_order:
        if n is None and sp not in color_dict:
            color_dict[sp] = spacer_color

    # 7) Create figure with 3 columns
    fig, axes = plt.subplots(1, 3, figsize=figure_size)
    ax_life, ax_dist, ax_hbond = axes

    # ----- Column 1: Lifetime Box Plot -----
    if df_life2 is not None and not df_life2.empty:
        sns.boxplot(
            data=df_life2,
            x="lifetime_ns",
            y="species",
            order=full_order,
            showmeans=True,
            meanprops={
                'marker': 'o',
                'markerfacecolor': 'black',
                'markeredgecolor': 'black',
                'markersize': 8
            },
            orient='h',
            palette=color_dict,
            showfliers=False,
            ax=ax_life
        )
        ax_life.set_title("Lifetime (Box)")
        ax_life.set_xlabel("Lifetime (ns)")
        ax_life.set_ylabel("Species")
        # Hide the NaN dummy labels
        y_labels = ax_life.get_yticklabels()
        new_ylab = []
        for lbl in y_labels:
            txt = lbl.get_text()
            if txt.startswith(spacer_prefix):
                new_ylab.append("")
            else:
                new_ylab.append(txt)
        ax_life.set_yticklabels(new_ylab)
    else:
        ax_life.set_title("No Lifetime Data")
        ax_life.set_xlabel("")
        ax_life.set_ylabel("")

    # ----- Column 2: Distance Box Plot (expanded) -----
    if not df_dist2.empty:
        sns.boxplot(
            data=df_dist2,
            x="distance_nm",
            y="species",
            showmeans=True,
            meanprops={
                'marker': 'o',
                'markerfacecolor': 'black',
                'markeredgecolor': 'black',
                'markersize': 8
            },
            order=full_order,
            orient='h',
            palette=color_dict,
            showfliers=False,
            ax=ax_dist
        )
        ax_dist.set_title("Distance (Box)")
        ax_dist.set_xlabel("Distance (nm)")
        ax_dist.set_ylabel("")
        # Hide ytick labels
        ax_dist.set_yticklabels([])
    else:
        ax_dist.set_title("No Distance Data")
        ax_dist.set_xlabel("")
        ax_dist.set_ylabel("")
        ax_dist.set_yticklabels([])

    # ----- Column 3: H-bond Count Violin -----
    if df_hbond2 is not None and not df_hbond2.empty:
        sns.violinplot(
            data=df_hbond2,
            x="hbond_count",
            y="species",
            order=full_order,
            orient='h',
            inner=None,
            palette=color_dict,
            cut=0,
            ax=ax_hbond
        )
        ax_hbond.set_title("H-bond Count (Violin)")
        ax_hbond.set_xlabel("H-bond Count")
        ax_hbond.set_ylabel("")
        ax_hbond.set_yticklabels([])

        # ----- Overlay Mean Values -----
        # Calculate mean hbond_count per species
        df_means = df_hbond2.groupby('species')['hbond_count'].mean().reset_index()

        # Remove spacer rows from df_means
        df_means = df_means[~df_means['species'].str.startswith(spacer_prefix)]

        # Plot mean values as diamonds
        ax_hbond.scatter(
            df_means['hbond_count'],
            df_means['species'],
            color='black',
            marker='o',  
            s=80,
            label='Mean',
            zorder=5  # Ensure it appears above the violins
        )

        # Add a legend for the mean
        ax_hbond.legend(loc='lower right')
    else:
        ax_hbond.set_title("No H-bond Data")
        ax_hbond.set_xlabel("")
        ax_hbond.set_ylabel("")
        ax_hbond.set_yticklabels([])

    plt.tight_layout()

    out_path = f"{output_prefix}_3col_naive_space.pdf"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved single 3-column figure with bracket grouping to '{out_path}'!")

def main():
    # Determine the script directory and project path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_path = find_project_root(current_dir=script_dir)

    # Add project path to sys.path
    sys.path.append(project_path)

    # Define the base process directory
    process_dir = os.path.join(project_path, "process.nobackup")

    time_per_frame_ns = 0.01

    # Prepare empty lists to gather data from all species
    hbonds_data = []     # (species, num_ones, hbond_count)
    lifetimes_data = []  # (species, num_ones, lifetime_ns)
    distance_data = []   # (species, num_ones, distance_nm, frequency)

    process_path = Path(process_dir)
    if not process_path.is_dir():
        raise ValueError(f"'{process_dir}' is not a valid directory.")

    # Example color dictionary (ensure this is defined or imported appropriately)
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
                if not line or line.startswith(('#', '@')):
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
    print(df_hbonds)
    df_lifetimes = pd.DataFrame(lifetimes_data, columns=["species", "num_ones", "lifetime_ns"])
    print(df_lifetimes)
    df_distance = pd.DataFrame(distance_data, columns=["species", "num_ones", "distance_nm", "frequency"])
    print(df_distance)

    # Generate the plots
    one_row_three_columns_box_violin(
        df_lifetimes=df_lifetimes,
        df_distance=df_distance,
        df_hbonds=df_hbonds,
        color_dict=color_dict if color_dict else None
    )

if __name__ == "__main__":
    main()