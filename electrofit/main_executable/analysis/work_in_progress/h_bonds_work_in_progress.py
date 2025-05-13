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


def visualize_data_donor_accpetor_pair(xpm_file='intra_hb_matrix.xpm',
                                       hbond_df=None,
                                       output_prefix='intra',
                                       bin_size_ns=0.1,
                                       time_per_frame_ns=0.01):
    """
    Visualizes hydrogen bond data with time binning, using donor-acceptor pairs 
    as y-axis labels.

    Parameters
    ----------
    xpm_file : str
        Path to the XPM file.
    hbond_df : pd.DataFrame
        DataFrame containing donor and acceptor labels. 
        Must have columns 'donor' and 'acceptor'.
    output_prefix : str
        Prefix for the output file names.
    bin_size_ns : float
        Size of each time bin in nanoseconds.
    time_per_frame_ns : float
        Time duration of each frame in nanoseconds.

    Returns
    -------
    None
        Saves multiple plots as files (PDF) in the current directory.
    """
    data_matrix, metadata = parse_xpm(xpm_file)
    analysis_results = analyze_hydrogen_bonds(data_matrix, metadata)

    matrix_shape = np.shape(data_matrix)
    print(f"Data matrix shape: {matrix_shape}")

    # Extract donor-acceptor labels
    if hbond_df is not None:
        hbond_labels = hbond_df['donor'] + '-' + hbond_df['acceptor']
    else:
        hbond_labels = [str(i) for i in range(matrix_shape[0])]

    # Verify alignment
    if len(hbond_labels) != data_matrix.shape[0]:
        print(f"Number of labels: {len(hbond_labels)}, "
              f"Number of rows in data_matrix: {data_matrix.shape[0]}")
        raise ValueError("Mismatch between number of hydrogen bonds in "
                         "hbond_df and data_matrix.")

    # Calculate number of frames per bin
    frames_per_bin = int(bin_size_ns / time_per_frame_ns)
    if frames_per_bin <= 0:
        raise ValueError("Bin size must be larger than the time per frame.")

    # Number of bins (account for partial bin if it doesn't divide evenly)
    num_bins = data_matrix.shape[1] // frames_per_bin
    if data_matrix.shape[1] % frames_per_bin != 0:
        num_bins += 1

    # Initialize binned data matrix
    binned_matrix = np.zeros((data_matrix.shape[0], num_bins))

    # Aggregate data into bins
    for i in range(num_bins):
        start_idx = i * frames_per_bin
        end_idx = min((i + 1) * frames_per_bin, data_matrix.shape[1])
        bin_data = data_matrix[:, start_idx:end_idx]
        binned_matrix[:, i] = np.mean(bin_data, axis=1)

    # Adjust figure size based on the number of hydrogen bonds
    fig_height = max(6, 0.3 * len(hbond_labels))
    plt.figure(figsize=(8, fig_height))

    # Heatmap of hydrogen bonds
    plt.imshow(binned_matrix, aspect='auto', cmap="Reds",
               origin='lower', vmin=0, vmax=1)
    plt.title(f"{metadata.get('title', 'Hydrogen Bond Existence Map')} (Binned)")
    plt.xlabel('Time (ns)')
    plt.ylabel('Donor-Acceptor Pairs')

    # Adjust time axis labels
    bin_times_ns = np.arange(num_bins) * bin_size_ns
    num_ticks = 5
    tick_positions = np.linspace(0, num_bins - 1, num_ticks, dtype=int)
    tick_labels = [f"{bin_times_ns[pos]:.1f}" for pos in tick_positions]
    plt.xticks(tick_positions, labels=tick_labels)

    plt.yticks(np.arange(len(hbond_labels)), hbond_labels)

    plt.tick_params(axis='y', which='major', labelsize=8)

    # Color bar for fraction of bond presence
    cbar = plt.colorbar(label='Fraction of Bond Presence')
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar.set_ticklabels(['0%', '25%', '50%', '75%', '100%'])

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_hb_existence_map_binned.pdf")
    plt.close()

    print("Hydrogen bonds per index:", analysis_results['hbonds_per_index'])

    # Histogram: Hydrogen Bonds per Index (Horizontal Bar Plot)
    plt.figure(figsize=(8, fig_height))
    plt.barh(range(len(analysis_results['hbonds_per_index'])),
             analysis_results['hbonds_per_index'],
             color="darkred")
    plt.title('Total Occurrence of Each Hydrogen Bond')
    plt.xlabel('Occurrences')
    plt.ylabel('Donor-Acceptor Pairs')
    plt.yticks(np.arange(len(hbond_labels)), hbond_labels)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_hb_occurrences.pdf")
    plt.close()

    # Lifetime Distribution
    # Each frame is 0.01 ns (10 ps)
    time_per_frame = 0.01

    # Flatten all lifetimes
    all_lifetimes_frames = [
        lifetime
        for bond_lifetimes in analysis_results['lifetimes']
        for lifetime in bond_lifetimes
    ]

    # Convert frames to ns
    all_lifetimes_ns = [l * time_per_frame for l in all_lifetimes_frames]

    max_lifetime = max(all_lifetimes_ns) if all_lifetimes_ns else 0
    bin_width = 0.01  # 0.01 ns bin width
    if max_lifetime > 0:
        num_bins_lifetime = int(max_lifetime / bin_width) + 2
    else:
        num_bins_lifetime = 1
    bins = np.arange(0, num_bins_lifetime * bin_width, bin_width)

    plt.figure(figsize=(8, 6))
    plt.hist(all_lifetimes_ns, bins=bins, edgecolor=None, color="darkred")

    plt.title('Hydrogen Bond Lifetime Distribution')
    plt.xlabel('Lifetime (ns)')
    plt.ylabel('Frequency')
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"{output_prefix}_hb_lifetime.pdf")
    plt.close()


def visualize_data_donor_acceptor(xpm_file='inter_hb_matrix.xpm',
                                  hbond_df=None,
                                  output_prefix='inter',
                                  bin_size_ns=0.1,
                                  time_per_frame_ns=0.01):
    """
    Visualizes hydrogen bond data and analysis with time binning, focusing on
    oxygen atoms in the molecule (both donor & acceptor). Aggregates data to
    show if any H-bond is present for a given oxygen, across time bins.

    Parameters
    ----------
    xpm_file : str
        Path to the XPM file.
    hbond_df : pd.DataFrame
        DataFrame containing donor and acceptor labels (must include 'donor'
        and 'acceptor' columns).
    output_prefix : str
        Prefix for output file names.
    bin_size_ns : float
        Size of each time bin in nanoseconds.
    time_per_frame_ns : float
        Duration of each frame in nanoseconds.

    Returns
    -------
    None
        Saves several plots (heatmaps, histograms) to the current directory.
    """
    data_matrix, metadata = parse_xpm(xpm_file)
    analysis_results = analyze_hydrogen_bonds(data_matrix, metadata)

    matrix_shape = np.shape(data_matrix)
    print(f"Data matrix shape: {matrix_shape}")

    # Ensure hbond_df is provided
    if hbond_df is None:
        raise ValueError("hbond_df must be provided.")

    hbond_df = hbond_df.reset_index(drop=True)
    hbond_df['idx'] = hbond_df.index

    # Identify oxygen atoms in the molecule.
    # Here, we define a function to check if an atom is from the molecule:
    def is_molecule_atom(atom_name):
        return bool(re.match(r'^O\d+$', atom_name))

    molecule_donors = hbond_df['donor'][hbond_df['donor'].apply(is_molecule_atom)]
    molecule_acceptors = hbond_df['acceptor'][hbond_df['acceptor'].apply(is_molecule_atom)]
    molecule_oxygens = set(molecule_donors).union(set(molecule_acceptors))

    # Create mapping from oxygen to indices
    oxygen_to_indices = {}
    for oxygen in molecule_oxygens:
        indices = hbond_df[
            (hbond_df['donor'] == oxygen) | (hbond_df['acceptor'] == oxygen)
        ]['idx'].tolist()
        oxygen_to_indices[oxygen] = indices

    frames_per_bin = int(bin_size_ns / time_per_frame_ns)
    if frames_per_bin <= 0:
        raise ValueError("Bin size must be larger than the time per frame.")

    num_bins = data_matrix.shape[1] // frames_per_bin
    if data_matrix.shape[1] % frames_per_bin != 0:
        num_bins += 1

    binned_matrix = np.zeros((data_matrix.shape[0], num_bins))
    for i in range(num_bins):
        start_idx = i * frames_per_bin
        end_idx = min((i + 1) * frames_per_bin, data_matrix.shape[1])
        bin_data = data_matrix[:, start_idx:end_idx]
        binned_matrix[:, i] = np.mean(bin_data, axis=1)

    # Aggregate data per oxygen
    num_oxygens = len(oxygen_to_indices)
    aggregated_binned_matrix = np.zeros((num_oxygens, num_bins))
    oxygen_list = list(oxygen_to_indices.keys())

    for i, oxygen in enumerate(oxygen_list):
        indices = oxygen_to_indices[oxygen]
        oxygen_data = binned_matrix[indices, :]  # shape: (N_bonds_for_oxygen, num_bins)
        # Use max to indicate presence in a given bin
        aggregated_data = np.max(oxygen_data, axis=0)
        aggregated_binned_matrix[i, :] = aggregated_data

    # Prepare roles for coloring y-tick labels
    oxygen_roles = {}
    for oxygen in oxygen_list:
        is_donor = any(hbond_df['donor'] == oxygen)
        is_acceptor = any(hbond_df['acceptor'] == oxygen)
        if is_donor and is_acceptor:
            role = 'both'
        elif is_donor:
            role = 'donor'
        elif is_acceptor:
            role = 'acceptor'
        else:
            role = 'unknown'
        oxygen_roles[oxygen] = role

    fig_height = max(6, 0.3 * len(oxygen_list))
    plt.figure(figsize=(8, fig_height))

    plt.imshow(aggregated_binned_matrix, aspect='auto', cmap="Reds",
               origin='lower', vmin=0, vmax=1)
    plt.title(f"{metadata.get('title', 'Hydrogen Bond Existence Map')} (Binned)")
    plt.xlabel('Time (ns)')
    plt.ylabel('Oxygen Atoms')

    bin_times_ns = np.arange(num_bins) * bin_size_ns
    num_ticks = 5
    tick_positions = np.linspace(0, num_bins - 1, num_ticks, dtype=int)
    tick_labels = [f"{bin_times_ns[pos]:.1f}" for pos in tick_positions]
    plt.xticks(tick_positions, labels=tick_labels)

    plt.yticks(np.arange(len(oxygen_list)), oxygen_list)
    plt.tick_params(axis='y', which='major', labelsize=8)

    # Color y-tick labels
    ax = plt.gca()
    yticks = ax.get_yticklabels()
    for tick_label in yticks:
        text = tick_label.get_text()
        role = oxygen_roles.get(text, 'unknown')
        if role == 'donor':
            tick_label.set_color('blue')
        elif role == 'acceptor':
            tick_label.set_color('green')
        elif role == 'both':
            tick_label.set_color('purple')

    cbar = plt.colorbar(label='Fraction of Bond Presence')
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cbar.set_ticklabels(['0%', '25%', '50%', '75%', '100%'])

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_hb_existence_map_binned.pdf")
    plt.close()

    # Histogram of occurrences per oxygen atom
    oxygen_occurrences = np.sum(aggregated_binned_matrix * frames_per_bin, axis=1)
    plt.figure(figsize=(8, fig_height))
    plt.barh(range(len(oxygen_occurrences)), oxygen_occurrences, color="darkred")
    plt.title('Total Occurrence of Hydrogen Bonds per Oxygen Atom')
    plt.xlabel('Occurrences')
    plt.ylabel('Oxygen Atoms')
    plt.yticks(np.arange(len(oxygen_list)), oxygen_list)
    plt.grid(False)

    ax = plt.gca()
    yticks = ax.get_yticklabels()
    for tick_label in yticks:
        text = tick_label.get_text()
        role = oxygen_roles.get(text, 'unknown')
        if role == 'donor':
            tick_label.set_color('blue')
        elif role == 'acceptor':
            tick_label.set_color('green')
        elif role == 'both':
            tick_label.set_color('purple')

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_hb_occurrences.pdf")
    plt.close()


def plot_hb_dist_xvg(file_path, plot_file_name="hb_distribution.pdf"):
    """
    Reads an .xvg file and plots the data using matplotlib.

    Parameters
    ----------
    file_path : str
        Path to the XVG file.
    plot_file_name : str
        The output filename (PDF), default "hb_distribution.pdf".

    Returns
    -------
    None
        Saves the plot as a PDF.
    """
    title = "XVG Plot"
    x_label = "X-axis"
    y_label = "Y-axis"

    data = []
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()

            # Skip empty lines
            if not line:
                continue

            # Handle metadata lines
            if line.startswith('#'):
                continue
            elif line.startswith('@'):
                tokens = line.split()
                if len(tokens) >= 3:
                    if tokens[1] == "title":
                        title = " ".join(tokens[2:]).strip('"')
                    elif tokens[1] == "xaxis" and tokens[2] == "label":
                        x_label = " ".join(tokens[3:]).strip('"')
                    elif tokens[1] == "yaxis" and tokens[2] == "label":
                        y_label = " ".join(tokens[3:]).strip('"')
            else:
                # Data lines
                try:
                    x, y = map(float, line.split())
                    data.append((x, y))
                except ValueError:
                    # Handle lines that do not have exactly two float numbers
                    continue

    if not data:
        print(f"No data found in the XVG file ({file_path}).")
        return

    data = np.array(data)
    x = data[:, 0]
    y = data[:, 1]

    plt.figure(figsize=(8, 6))
    plt.plot(x, y, marker='o', linestyle='-', color='darkblue')
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid(False)
    plt.tight_layout()
    plt.savefig(f"{plot_file_name}.pdf", format="pdf")
    plt.close()


def add_dashed_hbonds(svg_content, hbond_pairs, atom_coords,
                      dasharray='5,5', stroke_width='2',
                      angle_threshold=5, offset_step=5):
    """
    Adds dashed lines between specified atom pairs in the SVG content,
    adjusting overlapping lines to prevent them from overlaying.

    The lines have a gradient from black (donor) to light grey (acceptor).

    Parameters
    ----------
    svg_content : str
        The original SVG content.
    hbond_pairs : list of tuples
        List of atom index pairs representing H-bonds (donor_idx, acceptor_idx).
    atom_coords : dict
        Mapping atom indices to (x, y) drawing coordinates in the SVG.
    dasharray : str
        Dash pattern for the dashed lines (e.g., '5,5').
    stroke_width : str
        Width of the dashed lines (as a string, e.g., '2').
    angle_threshold : float
        Threshold angle in degrees to consider lines "overlapping".
    offset_step : float
        Step size (in drawing units) for offsetting overlapping lines.

    Returns
    -------
    str
        Modified SVG content with dashed H-bond lines inserted after the
        </rect> element(s).
    """
    import uuid

    lines_data = []
    gradient_defs = []

    for idx, (donor_idx, acceptor_idx) in enumerate(hbond_pairs):
        if donor_idx not in atom_coords or acceptor_idx not in atom_coords:
            print(f"Warning: Atom indices {donor_idx} or {acceptor_idx} "
                  f"not found in molecule.")
            continue

        x1, y1 = atom_coords[donor_idx]
        x2, y2 = atom_coords[acceptor_idx]

        dx = x2 - x1
        dy = y2 - y1
        length = math.hypot(dx, dy)
        if length == 0:
            print(f"Warning: Zero length line between atom {donor_idx} "
                  f"and {acceptor_idx}. Skipping line.")
            continue

        dir_vector = (dx / length, dy / length)

        line_data = {
            'donor_idx': donor_idx,
            'acceptor_idx': acceptor_idx,
            'x1': x1,
            'y1': y1,
            'x2': x2,
            'y2': y2,
            'dir_vector': dir_vector,
            'offset': 0,
            'gradient_id': f'grad_{uuid.uuid4().hex}'
        }
        lines_data.append(line_data)

    # Detect overlapping lines and apply offsets
    for (line1, line2) in combinations(lines_data, 2):
        common_atoms = {
            line1['donor_idx'], line1['acceptor_idx']
        }.intersection({
            line2['donor_idx'], line2['acceptor_idx']
        })
        if common_atoms:
            dot_product = (line1['dir_vector'][0] * line2['dir_vector'][0] +
                           line1['dir_vector'][1] * line2['dir_vector'][1])
            # Clamp dot_product between -1 and 1 for acos
            dot_product = max(-1.0, min(1.0, dot_product))
            angle_rad = math.acos(dot_product)
            angle_deg = math.degrees(angle_rad)

            # If lines are nearly colinear
            if abs(angle_deg) < angle_threshold or abs(angle_deg - 180) < angle_threshold:
                line1['offset'] += offset_step
                line2['offset'] -= offset_step

    dashed_lines = []
    for line_data in lines_data:
        x1, y1 = line_data['x1'], line_data['y1']
        x2, y2 = line_data['x2'], line_data['y2']
        dx = x2 - x1
        dy = y2 - y1
        length = math.hypot(dx, dy)
        nx = -dy / length
        ny = dx / length

        offset_distance = line_data['offset']
        x1_adj = x1 + nx * offset_distance
        y1_adj = y1 + ny * offset_distance
        x2_adj = x2 + nx * offset_distance
        y2_adj = y2 + ny * offset_distance

        donor_color_hex = '#000000'  # Black
        acceptor_color_hex = '#D3D3D3'  # Light grey

        gradient_def = f'''
        <linearGradient id="{line_data['gradient_id']}"
                        gradientUnits="userSpaceOnUse"
                        x1="{x1_adj}" y1="{y1_adj}"
                        x2="{x2_adj}" y2="{y2_adj}">
            <stop offset="0%" stop-color="{donor_color_hex}" />
            <stop offset="100%" stop-color="{acceptor_color_hex}" />
        </linearGradient>
        '''
        gradient_defs.append(gradient_def)

        line_svg = (f'<line x1="{x1_adj}" y1="{y1_adj}" '
                    f'x2="{x2_adj}" y2="{y2_adj}" '
                    f'stroke="url(#{line_data["gradient_id"]})" '
                    f'stroke-width="{stroke_width}" '
                    f'stroke-dasharray="{dasharray}" />')
        dashed_lines.append(line_svg)

    if dashed_lines:
        gradients_svg = '\n'.join(gradient_defs)

        # Insert gradient definitions into the SVG
        header_end_pos = svg_content.find('<!-- END OF HEADER -->')
        if header_end_pos != -1:
            svg_content = (svg_content[:header_end_pos] +
                           f'\n<defs>\n{gradients_svg}\n</defs>\n' +
                           svg_content[header_end_pos:])
        else:
            insert_pos = svg_content.find('>') + 1
            svg_content = (svg_content[:insert_pos] +
                           f'\n<defs>\n{gradients_svg}\n</defs>\n' +
                           svg_content[insert_pos:])
            print("Warning: <!-- END OF HEADER --> comment not found. "
                  "Inserting <defs> after <svg> tag.")

        dashed_svg_content = '\n'.join(dashed_lines)

        insert_pos = svg_content.find('</rect>') + len('</rect>')
        if insert_pos == -1 + len('</rect>'):
            # </rect> not found, default to inserting after </defs>
            insert_pos = svg_content.find('</defs>') + len('</defs>')
            if insert_pos == -1 + len('</defs>'):
                insert_pos = svg_content.find('>') + 1
                print("Warning: Neither </rect> nor </defs> found. "
                      "Inserting dashed lines after <svg> tag.")
        svg_content = (svg_content[:insert_pos] +
                       f'\n{dashed_svg_content}\n' +
                       svg_content[insert_pos:])
        print(f"Added {len(dashed_lines)} dashed H-bond line(s) with gradients to SVG.")
    else:
        print("No dashed lines were added to the SVG content.")

    return svg_content


def create_legend(x, y, legend_entries, font_size=12):
    """
    Creates SVG elements for a legend.

    Parameters
    ----------
    x : float
        X-coordinate of the top-left corner of the legend.
    y : float
        Y-coordinate of the top-left corner of the legend.
    legend_entries : list of dict
        Legend entries, each with keys like:
            - 'label' (string)
            - 'color' (string or hex color code)
            - 'dasharray' (optional, for dashed line entries)
    font_size : int
        Font size for legend text.

    Returns
    -------
    str
        SVG snippet representing the legend items.
    """
    legend_svg = []
    entry_height = font_size + 10
    for i, entry in enumerate(legend_entries):
        entry_y = y + i * entry_height

        if 'dasharray' in entry:
            line = (f'<line x1="{x}" y1="{entry_y + font_size / 2}" '
                    f'x2="{x + 20}" y2="{entry_y + font_size / 2}" '
                    f'stroke="{entry["color"]}" stroke-width="2" '
                    f'stroke-dasharray="{entry["dasharray"]}" />')
            legend_svg.append(line)
        else:
            rect = (f'<rect x="{x}" y="{entry_y}" width="20" '
                    f'height="{font_size}" fill="{entry["color"]}" '
                    f'stroke="black" />')
            legend_svg.append(rect)

        text_x = x + 25
        text_y = entry_y + font_size
        text = (f'<text x="{text_x}" y="{text_y}" '
                f'font-size="{font_size}" font-family="Arial">'
                f'{entry["label"]}</text>')
        legend_svg.append(text)

    return '\n'.join(legend_svg)


def add_legend_to_svg(svg_content, legend_svg_content):
    """
    Inserts the legend SVG content into the main SVG content before </svg>.

    Parameters
    ----------
    svg_content : str
        The original SVG content.
    legend_svg_content : str
        SVG snippet for the legend.

    Returns
    -------
    str
        Modified SVG with the legend inserted.
    """
    insert_pos = svg_content.rfind('</svg>')
    if insert_pos == -1:
        print("Warning: </svg> tag not found. Legend will not be added.")
        return svg_content

    svg_content = (svg_content[:insert_pos] +
                   f'\n{legend_svg_content}\n' +
                   svg_content[insert_pos:])
    print("Legend added to SVG.")
    return svg_content


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

def generate_intra_hbond_summary(process_dir, output_file="intra_hbonds_summary.pdf"):
    """
    Generates a single summary figure (3 columns x N rows) of intramolecular
    hydrogen-bond data for all species found under `process_dir`.
    
    Each row corresponds to one species/folder. The columns are:
      1) Box plot (distribution of # of H-bonds over time)
      2) Lifetime histogram
      3) Distance distribution
    
    The species identifier is derived from the folder name by stripping 'IP_'.
    
    Parameters
    ----------
    process_dir : str
        Path to the main directory containing subfolders for each species.
    output_file : str
        Filename (PDF) to save the final summary figure.
    """
    
    # We will collect data for each species in this list of dicts
    species_data_list = []
    
    # Loop over all subdirectories in `process_dir` to find 'IP_XXXXXX' folders
    for folder_name in sorted(os.listdir(process_dir)):
        folder_path = os.path.join(process_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue
        
        # We assume species folders are named like "IP_101010", etc.
        if not folder_name.startswith("IP_"):
            continue
        
        # Extract the binary code (identifier) by stripping "IP_"
        identifier = folder_name.replace("IP_", "")
        
        # Path to 'analyze_final_sim/h_bonds' 
        analyze_dir = os.path.join(folder_path, "analyze_final_sim", "h_bonds")
        
        # The intramolecular H-bond files we need
        intra_num_xvg = os.path.join(analyze_dir, "intra_hb_num.xvg")
        intra_matrix_xpm = os.path.join(analyze_dir, "intra_hb_matrix.xpm")
        intra_dist_xvg = os.path.join(analyze_dir, "intra_hb_dist.xvg")
        
        # Make sure these files exist (otherwise skip or warn)
        if (not os.path.isfile(intra_num_xvg) or 
            not os.path.isfile(intra_matrix_xpm) or 
            not os.path.isfile(intra_dist_xvg)):
            print(f"Warning: Missing intra-hbond files for folder {folder_name}. Skipped.")
            continue
        
        # -------------------------------------------------------
        # (1) Parse number of H-bonds vs. time
        # -------------------------------------------------------
        data_num = load_hb_num_xvg(intra_num_xvg)
        # Example: column 0 is time in ps, column 1 is # H-bonds
        # If there's a third column, column 2 might be pairs within 0.35 nm
        if data_num.size < 2:
            print(f"Warning: No valid data in {intra_num_xvg}. Skipped.")
            continue
        
        # The second column is the distribution of # of intramolecular H-bonds
        # over time. We'll store that for the box plot.
        hbonds_distribution = data_num[:, 1]
        
        # -------------------------------------------------------
        # (2) Compute lifetimes from the .xpm
        # -------------------------------------------------------
        data_matrix, meta = parse_xpm(intra_matrix_xpm)
        analysis = analyze_hydrogen_bonds(data_matrix, meta)
        
        # Flatten out all lifetimes (in frames)
        # analysis['lifetimes'] is a list of lists, one sublist per H-bond
        all_frames_lifetimes = [
            lf for bond_lifetimes in analysis['lifetimes']
               for lf in bond_lifetimes
        ]
        
        # Convert frames to ns if you know the time per frame
        # In your existing script, each frame is 0.01 ns.
        time_per_frame_ns = 0.01
        lifetimes_ns = [x * time_per_frame_ns for x in all_frames_lifetimes]
        
        # -------------------------------------------------------
        # (3) Distance distribution from .xvg
        # -------------------------------------------------------
        # Typically, "intra_hb_dist.xvg" has 2 columns: distance and frequency
        dist_data = []
        with open(intra_dist_xvg, 'r') as fdist:
            for line in fdist:
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
            print(f"Warning: No valid distance data in {intra_dist_xvg}. Skipped.")
            continue
        
        # X is the distance, Y is the frequency (or count)
        dist_x = dist_data[:, 0]
        dist_y = dist_data[:, 1]
        
        # Add everything to our species_data_list
        species_data_list.append({
            'identifier': identifier,
            'hbonds_distribution': hbonds_distribution,
            'lifetimes_ns': lifetimes_ns,
            'dist_x': dist_x,
            'dist_y': dist_y
        })
    
    # If no data found, just return
    if not species_data_list:
        print("No species data found. Exiting without plotting.")
        return
    
    # We expect 11 species, but let's handle however many we found
    n_species = len(species_data_list)
    
    # Create the figure with 3 columns and n_species rows
    fig, axes = plt.subplots(nrows=n_species, ncols=3, figsize=(18, 4*n_species))
    
    # In case n_species=1, axes is not 2D. Ensure 2D indexing:
    if n_species == 1:
        # Expand dimensions so we can always do axes[i, j]
        axes = np.array([axes])
    
    for i, data_dict in enumerate(species_data_list):
        identifier = data_dict['identifier']
        hbonds_dist = data_dict['hbonds_distribution']
        lifetimes_ns = data_dict['lifetimes_ns']
        dist_x = data_dict['dist_x']
        dist_y = data_dict['dist_y']
        
        # 1) Box plot of # of H-bonds distribution
        ax_box = axes[i, 0]
        ax_box.boxplot(hbonds_dist, vert=True)
        ax_box.set_title(f"H-Bond Count (Box)\n{identifier}")
        ax_box.set_ylabel("H-bond count")
        # X-axis label is trivial; just show species code or nothing
        ax_box.set_xticks([])
        
        # 2) Lifetime histogram
        ax_life = axes[i, 1]
        if len(lifetimes_ns) > 0:
            ax_life.hist(lifetimes_ns, bins=30, color='darkred', alpha=0.7)
        ax_life.set_title(f"Lifetime Dist\n{identifier}")
        ax_life.set_xlabel("Lifetime (ns)")
        ax_life.set_ylabel("Frequency")
        
        # 3) Distance distribution
        ax_dist = axes[i, 2]
        ax_dist.plot(dist_x, dist_y, color='darkblue')
        ax_dist.set_title(f"D–A Distance\n{identifier}")
        ax_dist.set_xlabel("Distance (nm)")
        ax_dist.set_ylabel("Frequency")
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Summary figure saved to {output_file}!")

script_dir = os.path.dirname(os.path.abspath(__file__))
project_path = find_project_root(current_dir=script_dir)

sys.path.append(project_path)


from electrofit.helper.file_manipulation import find_file_with_extension
from electrofit.commands.run_commands import run_command

# Define the base process directory
process_dir = os.path.join(project_path, "process")

# Loop through each subdirectory in the process directory
for folder_name in os.listdir(process_dir):

    folder_path = os.path.join(process_dir, folder_name)

    # Check if it's a directory
    if os.path.isdir(folder_path):
        # Define the 'run_final_gmx_simulation' directory within this folder
        run_final_sim_dir = os.path.join(folder_path, 'run_final_gmx_simulation')
        run_gau_create_gmx_in_dir = os.path.join(folder_path, 'run_gau_create_gmx_in')

        # Check if 'run_final_gmx_simulation' exists
        if os.path.isdir(run_final_sim_dir):
            # Define the destination directory 'analyze_final_sim/h_bonds'
            dest_dir = os.path.join(folder_path, 'analyze_final_sim', 'h_bonds')
            os.makedirs(dest_dir, exist_ok=True)
            os.chdir(dest_dir)

            gro_file_path = os.path.join(run_final_sim_dir, 'md.tpr')
            xtc_file_path = os.path.join(run_final_sim_dir, 'md_center.xtc')

            # GROMACS commands to generate intra/inter H-bond data
            run_command(
                f'echo "2\n2\n" | gmx hbond -s {gro_file_path} -f {xtc_file_path} '
                f'-hbn intra_hb_idx.ndx -num intra_hb_num.xvg -dist intra_hb_dist.xvg '
                f'-g intra_hb.log -hbm intra_hb_matrix.xpm',
                cwd=dest_dir
            )
            run_command(
                f'echo "2\n5\n" | gmx hbond -s {gro_file_path} -f {xtc_file_path} '
                f'-hbn inter_hb_idx.ndx -num inter_hb_num.xvg -dist inter_hb_dist.xvg '
                f'-g inter_hb.log -hbm inter_hb_matrix.xpm',
                cwd=dest_dir
            )

            # ------------ Plot inter ------------
            xvg_filename = 'inter_hb_num.xvg'
            data = load_hb_num_xvg(xvg_filename)

            if data.size == 0:
                raise ValueError(f"No data found in {xvg_filename}.")

            time = data[:, 0]  # Time in ps
            time = time / 10**3  # Convert to ns
            num_hbonds = data[:, 1]
            pairs_within_0_35_nm = data[:, 2] if data.shape[1] > 2 else None

            # Plot
            plt.figure(figsize=(8, 6))
            plt.plot(time, num_hbonds, label='Hydrogen Bonds', color='darkblue', linewidth=2)
            if pairs_within_0_35_nm is not None:
                plt.plot(time, pairs_within_0_35_nm, label='Pairs within 0.35 nm',
                         color='darkred', linewidth=2)

            plt.title('Hydrogen Bonds Over Time')
            plt.xlabel('Time (ns)')
            plt.ylabel('Number')
            plt.legend()
            plt.grid(False)

            output_svg = 'inter_hb_num_over_time.pdf'
            plt.savefig(output_svg, format='pdf')
            plt.close()
            print(f"Plot saved as {output_svg}")

            plot_hb_dist_xvg('inter_hb_dist.xvg', plot_file_name='inter_hb_distriution')

            inter_log_file = 'inter_hb.log'
            inter_hbond_df = parse_hbond_log_to_dataframe(inter_log_file)
            print(inter_hbond_df)

            visualize_data_donor_acceptor(xpm_file="inter_hb_matrix.xpm",
                                          hbond_df=inter_hbond_df,
                                          output_prefix='inter')

            # ------------ Plot intra -------------
            log_file = 'intra_hb.log'
            output_csv = 'intra_hbond_pairs.csv'
            hbond_df = parse_hbond_log_to_dataframe(log_file)

            print("Extracted Donor-Acceptor Pairs DataFrame:")
            print(hbond_df)
            hbond_df.to_csv(output_csv, index=False)
            print(f"\nDonor-Acceptor pairs have been saved to '{output_csv}'.")

            visualize_data_donor_acceptor(xpm_file="intra_hb_matrix.xpm",
                                          hbond_df=hbond_df,
                                          output_prefix='intra2')

            data_matrix, meta_data = parse_xpm("intra_hb_matrix.xpm")
            analysis_results = analyze_hydrogen_bonds(data_matrix=data_matrix,
                                                      metadata=meta_data)
            hbonds_per_index = analysis_results['hbonds_per_index']
            hbond_df['count'] = hbonds_per_index

            # For each donor, get the two most abundant hydrogen bonds
            top_hbonds_df = hbond_df.groupby('donor').apply(
                lambda x: x.nlargest(2, 'count')
            ).reset_index(drop=True)

            print("Top hydrogen bonds:")
            print(top_hbonds_df)

            # Load the molecule from MOL2
            os.chdir(run_gau_create_gmx_in_dir)
            mol2_file_name = find_file_with_extension('mol2')
            os.chdir(dest_dir)
            mol2_file = os.path.join(run_gau_create_gmx_in_dir, mol2_file_name)

            molecule = Chem.MolFromMol2File(mol2_file, removeHs=False)
            if molecule is None:
                raise ValueError(f"Failed to load molecule from {mol2_file}")
            print("Molecule loaded successfully!")

            # Generate 2D coordinates
            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(molecule)

            # Create custom labels for the atoms
            atom_counters = {}
            atom_labels = {}
            for atom in molecule.GetAtoms():
                symbol = atom.GetSymbol()
                idx = atom.GetIdx()
                atom_counters[symbol] = atom_counters.get(symbol, 0) + 1
                atom_labels[idx] = f"{symbol}{atom_counters[symbol]}"

            label_to_atom_idx = {label: idx for idx, label in atom_labels.items()}

            donor_atom_indices = []
            acceptor_atom_indices = []

            for _, row in hbond_df.iterrows():
                donor_label = row['donor']
                acceptor_label = row['acceptor']
                donor_idx = label_to_atom_idx.get(donor_label)
                acceptor_idx = label_to_atom_idx.get(acceptor_label)
                if donor_idx is not None:
                    donor_atom_indices.append(donor_idx)
                if acceptor_idx is not None:
                    acceptor_atom_indices.append(acceptor_idx)

            donor_atom_indices = list(set(donor_atom_indices))
            acceptor_atom_indices = list(set(acceptor_atom_indices))
            both_donor_acceptor_indices = set(donor_atom_indices) & set(acceptor_atom_indices)

            donor_only_indices = set(donor_atom_indices) - both_donor_acceptor_indices
            acceptor_only_indices = set(acceptor_atom_indices) - both_donor_acceptor_indices

            print("Donor only atom indices:", donor_only_indices)
            print("Acceptor only atom indices:", acceptor_only_indices)
            print("Both donor and acceptor atom indices:", both_donor_acceptor_indices)

            # Create hbond_pairs from top_hbonds_df
            hbond_pairs = []
            for _, row in top_hbonds_df.iterrows():
                donor_label = row['donor']
                acceptor_label = row['acceptor']
                donor_idx = label_to_atom_idx.get(donor_label)
                acceptor_idx = label_to_atom_idx.get(acceptor_label)
                if donor_idx is not None and acceptor_idx is not None:
                    hbond_pairs.append((donor_idx, acceptor_idx))
                else:
                    print(f"Warning: Atom label {donor_label} or {acceptor_label} not found in molecule.")

            print("Hydrogen bond pairs for dashed lines:", hbond_pairs)

            highlight_atoms = donor_only_indices.union(
                acceptor_only_indices).union(
                both_donor_acceptor_indices
            )

            highlightAtomColors = {}
            for idx in donor_only_indices:
                highlightAtomColors[idx] = (0.6, 0.8, 1)  # Light blue
            for idx in acceptor_only_indices:
                highlightAtomColors[idx] = (0.6, 1, 0.6)  # Light green
            for idx in both_donor_acceptor_indices:
                highlightAtomColors[idx] = (0.6, 0.9, 0.8)  # Teal

            svg_size = 500
            drawer = rdMolDraw2D.MolDraw2DSVG(svg_size, svg_size)
            opts = drawer.drawOptions()
            opts.addAtomIndices = False
            opts.addBondIndices = False
            opts.baseFontSize = 0.3

            for idx, label in atom_labels.items():
                opts.atomLabels[idx] = label

            rdMolDraw2D.PrepareMolForDrawing(molecule)

            drawer.DrawMolecule(molecule,
                                highlightAtoms=list(highlight_atoms),
                                highlightAtomColors=highlightAtomColors)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()

            # Extract drawing coordinates for each atom
            draw_coords = {}
            for atom_idx in range(molecule.GetNumAtoms()):
                point = drawer.GetDrawCoords(atom_idx)
                draw_coords[atom_idx] = (point.x, point.y)

            # Add dashed lines
            modified_svg = add_dashed_hbonds(svg, hbond_pairs, draw_coords)

            with open('molecule_dashed_hbonds.svg', 'w') as f:
                f.write(modified_svg)
            print("Modified SVG with dashed H-bonds saved as 'molecule_dashed_hbonds.svg'")

            # Create a legend
            def rgb_to_hex(rgb_tuple):
                return '#' + ''.join(f'{int(255 * x):02X}' for x in rgb_tuple)

            legend_entries = []
            if donor_only_indices:
                legend_entries.append({'label': 'Donor (D)',
                                       'color': rgb_to_hex((0.6, 0.8, 1))})
            if acceptor_only_indices:
                legend_entries.append({'label': 'Acceptor (A)',
                                       'color': rgb_to_hex((0.6, 1, 0.6))})
            if both_donor_acceptor_indices:
                legend_entries.append({'label': 'D/A',
                                       'color': rgb_to_hex((0.6, 0.9, 0.8))})
            if hbond_pairs:
                legend_entries.append({'label': 'H-Bond',
                                       'color': 'lightgray',
                                       'dasharray': '5,5'})

            legend_x = 20
            legend_y = 20
            legend_svg_content = create_legend(legend_x, legend_y, legend_entries, font_size=14)
            modified_svg_with_legend = add_legend_to_svg(modified_svg, legend_svg_content)

            with open('molecule_dashed_hbonds_with_legend.svg', 'w') as f:
                f.write(modified_svg_with_legend)
            print("Modified SVG with dashed H-bonds and legend saved as 'molecule_dashed_hbonds_with_legend.svg'")

            # ---- Plot intra using donor-acceptor pair view ----
            visualize_data_donor_accpetor_pair(xpm_file="intra_hb_matrix.xpm",
                                               hbond_df=hbond_df)

            # Now plot the raw "intra_hb_num.xvg" time series
            xvg_filename = 'intra_hb_num.xvg'
            data = load_hb_num_xvg(xvg_filename)

            if data.size == 0:
                raise ValueError(f"No data found in {xvg_filename}.")

            time = data[:, 0]  # Time in ps
            num_hbonds = data[:, 1]
            pairs_within_0_35_nm = data[:, 2] if data.shape[1] > 2 else None

            plt.figure(figsize=(8, 6))
            plt.plot(time, num_hbonds, label='Hydrogen Bonds',
                     color='darkblue', linewidth=2)
            if pairs_within_0_35_nm is not None:
                plt.plot(time, pairs_within_0_35_nm,
                         label='Pairs within 0.35 nm',
                         color='darkred', linewidth=2)

            plt.title('Hydrogen Bonds Over Time')
            plt.xlabel('Time (ps)')
            plt.ylabel('Number')
            plt.legend()
            plt.grid(False)

            output_svg = 'intra_hb_num_over_time.pdf'
            plt.savefig(output_svg, format='pdf')
            plt.close()
            print(f"Plot saved as {output_svg}")

            plot_hb_dist_xvg('intra_hb_dist.xvg', plot_file_name='intra_hb_distriution')