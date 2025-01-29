import os
import sys
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import gridspec

sns.set_context("talk")

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

    Parameters:
    - file_path: str, path to the XPM file.

    Returns:
    - data_matrix: np.ndarray, binary matrix representing hydrogen bonds.
    - metadata: dict, contains title, labels, etc.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Initialize variables
    metadata = {}
    data_lines = []
    color_map = {}
    header_found = False
    data_started = False

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
            continue  # Skip the static declaration line

        if line.startswith('"') and not header_found:
            # Header line containing dimensions and color info
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
            # Example: "   c #FFFFFF " or "o  c #FF0000 "
            match = re.match(r'(.{%d})\s+c\s+(\S+)' % chars_per_pixel, color_def)
            if match:
                symbol = match.group(1)
                color = match.group(2)
                color_map[symbol] = color
            if len(color_map) == num_colors:
                data_started = True
            continue

        if data_started:
            # Read data lines
            if line.startswith('"'):
                data_line = line.strip('",')
                data_lines.append(data_line)

    # Convert data lines to binary matrix
    data_matrix = np.zeros((height, width), dtype=int)

    for y, line in enumerate(data_lines):
        for x, char in enumerate(line):
            if char in color_map:
                if color_map[char] == '#FF0000':  # Present bond
                    data_matrix[y, x] = 1
                else:
                    data_matrix[y, x] = 0
            else:
                data_matrix[y, x] = 0  # Default to 0 if unknown

    return data_matrix, metadata

def analyze_hydrogen_bonds(data_matrix, metadata):
    """
    Analyzes the hydrogen bond data.

    Parameters:
    - data_matrix: np.ndarray, binary matrix representing hydrogen bonds.
    - metadata: dict, contains title, labels, etc.

    Returns:
    - analysis_results: dict, contains various analysis metrics.
    """
    analysis_results = {}

    # Total number of hydrogen bonds over time
    hbonds_over_time = np.sum(data_matrix, axis=0)

    # Total occurrence of each hydrogen bond
    hbonds_per_index = np.sum(data_matrix, axis=1)

    # Lifetime distribution (how long each bond persists)
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
        if current_lifetime > 0:
            bond_lifetimes.append(current_lifetime)
        lifetimes.append(bond_lifetimes)
    analysis_results['hbonds_over_time'] = hbonds_over_time
    analysis_results['hbonds_per_index'] = hbonds_per_index
    analysis_results['lifetimes'] = lifetimes

    return analysis_results

def parse_hbond_log_to_dataframe(file_path):
    """
    Parses a GROMACS .log file to extract donor-acceptor pairs with indices and refined atom names.

    Parameters:
    - file_path (str): Path to the .log file.

    Returns:
    - pd.DataFrame: DataFrame containing idx, donor, and acceptor columns.
    """
    hbond_pairs = []

    # Regular expression patterns
    # Pattern to match lines with donor and acceptor information
    # Example line: I211O10              -      I211O2  
    line_pattern = re.compile(r'^\s*(\S+)\s+-\s+(\S+)\s*$')

    # Pattern to extract atom name and index from a string like 'I211O10'
    # Assuming residue identifier is 'I211' and atom name is 'O10'
    atom_pattern = re.compile(r'^[A-Za-z]+\d+([A-Za-z]+\d*)$')

    with open(file_path, 'r') as file:
        for line_number, line in enumerate(file, start=1):
            line = line.strip()

            # Skip empty lines and header lines
            if not line or line.startswith('#') or line.startswith('"""') or line.startswith('*'):
                continue

            # Match the line with donor and acceptor
            match = line_pattern.match(line)
            if match:
                donor_full = match.group(1)     # e.g., 'I211O10'
                acceptor_full = match.group(2)  # e.g., 'I211O2'

                # Extract atom names using regex
                donor_match = atom_pattern.match(donor_full)
                acceptor_match = atom_pattern.match(acceptor_full)

                if donor_match and acceptor_match:
                    donor_atom = donor_match.group(1)       # e.g., 'O10'
                    acceptor_atom = acceptor_match.group(1) # e.g., 'O2'

                    # Refine atom names
                    refined_donor = refine_atom_name(donor_atom)
                    refined_acceptor = refine_atom_name(acceptor_atom)

                    # Combine donor and acceptor with colon
                    pair = {
                        'donor': refined_donor,
                        'acceptor': refined_acceptor
                    }
                    hbond_pairs.append(pair)
                else:
                    print(f"Warning (Line {line_number}): Couldn't parse atoms in line: {line}")
            else:
                print(f"Warning (Line {line_number}): Line didn't match expected format and was skipped: {line}")

    # Create DataFrame with index
    df = pd.DataFrame(hbond_pairs)
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'idx'}, inplace=True)

    # Adjust idx to start from 0
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
            root = current_dir  # Set root to the current project_name directory
        if parent_dir == current_dir:
            # We've reached the filesystem root
            if root is None:
                raise FileNotFoundError(f"Project root directory '{project_name}' not found.")
            return root  # Return the outermost match found
        current_dir = parent_dir

def count_oxygen_occurrences_from_matrix(hbond_df, hbonds_per_index):
    """
    Counts the number of occurrences of each oxygen atom based on the matrix data.

    Parameters
    ----------
    hbond_df : pd.DataFrame
        DataFrame containing 'donor', 'acceptor', and 'idx' columns.
    hbonds_per_index : np.ndarray
        Array containing the count of occurrences for each hydrogen bond pair.

    Returns
    -------
    pd.Series
        Series with oxygen atom names as index and their summed counts as values.
    """
    # Assign the counts to the DataFrame
    hbond_df['count'] = hbonds_per_index

    # Melt the DataFrame to have one row per donor or acceptor with associated count
    melted = hbond_df.melt(id_vars=['idx', 'count'], value_vars=['donor', 'acceptor'], value_name='oxygen')

    # Drop any NaN values
    melted = melted.dropna(subset=['oxygen'])

    # Group by oxygen atom and sum the counts
    counts = melted.groupby('oxygen')['count'].sum()

    return counts

# --- new and better (correct ordering)

import matplotlib.pyplot as plt
import seaborn as sns
import re

def generate_summary_plot(occurrence_data, existence_data, time_per_frame_ns=0.01, output_file='oxygen_occurrences_summary.pdf', folder_order=None):
    """
    Generates a summary plot with an existence map and multiple horizontal bar subplots showing 
    oxygen atom occurrences per species. Each species's row contains an existence map and an occurrence plot.

    Parameters
    ----------
    occurrence_data : dict
        Dictionary where keys are species IDs and values are pd.Series of oxygen counts.
    existence_data : dict
        Dictionary where keys are species IDs and values are np.ndarray of aggregated binned data.
    time_per_frame_ns : float
        Time duration per frame in nanoseconds.
    output_file : str
        Filename for the saved summary plot.
    folder_order : list, optional
        List of species IDs in the desired order. Overrides internal sorting if provided.

    Returns
    -------
    None
    """
    num_species = len(occurrence_data)
    
    if num_species == 0:
        print("No data to plot for summary.")
        return
    
    # Determine all unique oxygen atoms across all species for consistent y-axis
    all_oxygens = set()
    for counts in occurrence_data.values():
        all_oxygens.update(counts.index.tolist())
    all_oxygens = sorted(all_oxygens, key=lambda x: int(re.findall(r'\d+', x)[0]) if re.findall(r'\d+', x) else 0)
    
    # Debugging: Print the number and list of oxygen atoms
    print(f"Total unique oxygen atoms: {len(all_oxygens)}")
    print(f"Oxygen atoms: {all_oxygens}")
    
    # Sort species IDs based on provided folder_order or default to sorted order
    if folder_order:
        sorted_species_ids = [sp for sp in folder_order if sp in occurrence_data]
    else:
        sorted_species_ids = sorted(occurrence_data.keys())
    
    # Determine the maximum occurrence to set consistent x-axis limits
    max_occurrence = max([counts.max() for counts in occurrence_data.values()]) if occurrence_data else 0
    x_limit = max_occurrence * time_per_frame_ns * 1.1  # Add 10% padding

    # Define the figure with a GridSpec layout
    fig_height = max(6, 5.5 * num_species)  # Adjust height based on number of species
    fig, axes = plt.subplots(
        nrows=num_species, 
        ncols=2, 
        figsize=(16, fig_height), 
        constrained_layout=True, 
        gridspec_kw={'width_ratios': [1, 1]}
    )
    
    # If there's only one species, axes is not a list of lists, so make it a list
    if num_species == 1:
        axes = [axes]
    
    # Define a color palette for occurrences
    palette = sns.color_palette("viridis", num_species)
    
    for i, species_id in enumerate(sorted_species_ids):
        # --- Subplot 1: Occurrence Bar Plot (Now on the Left) ---
        ax_occ = axes[i][0]
        counts = occurrence_data[species_id].reindex(all_oxygens, fill_value=0)
        time = time_per_frame_ns * counts
        
        sns.barplot(
            x=time.values, 
            y=all_oxygens, 
            ax=ax_occ, 
            color=palette[i]  # Use 'color' instead of 'palette' to set a single color per species
        )
        
        ax_occ.set_title(f"H-Bond Occurrence Time - {species_id}", fontweight='bold')
        ax_occ.set_xlabel('Time (ns)')
        
        # Set x-axis limits and invert it
        ax_occ.set_xlim(x_limit, 0)  # Reverse the axis

        # Add time labels at the end of each bar
        for p in ax_occ.patches:
            width = p.get_width()
            if width > 0:
                ax_occ.annotate(
                    f"{width:.2f}",
                    (width, p.get_y() + p.get_height() / 2),
                    ha='right',  # Adjusted for inverted axis
                    va='center',
                    xytext=(-5, 0),  # Shift left
                    textcoords='offset points',
                    fontsize=13
                )
        
        # Enhance grid
        ax_occ.grid(True, linestyle='--', alpha=0.5, axis='x')
        
        # Hide y-tick labels for the occurrence plot to reduce clutter
        ax_occ.set_yticklabels([])
        ax_occ.yaxis.set_ticks_position('right')

        # --- Subplot 2: Existence Map (Now on the Right) ---
        ax_exist = axes[i][1]
        aggregated_binned_matrix = existence_data[species_id]
        
        # Debugging: Print the shape of the matrix
        print(f"Species ID: {species_id}, Aggregated Binned Matrix Shape: {aggregated_binned_matrix.shape}")
        
        if aggregated_binned_matrix.shape[0] != len(all_oxygens):
            print(f"Error: Shape mismatch for species '{species_id}'. Expected {len(all_oxygens)} rows, got {aggregated_binned_matrix.shape[0]}.")
            # Handle the mismatch, e.g., skip plotting this species or adjust the matrix
            # For now, let's skip plotting this species
            ax_exist.set_visible(False)
            continue
        
        sns.heatmap(
            aggregated_binned_matrix, 
            cmap="Reds", 
            cbar=True, 
            ax=ax_exist, 
            linewidths=0,         # Set linewidths to create grid lines between cells
            linecolor='white',      # Set linecolor for the grid lines
            cbar_kws={"label": "Fraction of Bond Presence"}, 
            vmin=0, 
            vmax=1
        )
        ax_exist.set_title(f"H-Bond Existence Map - {species_id}", fontweight='bold')
        ax_exist.set_xlabel('Time (ns)')
        
        # Adjust x-ticks based on the number of time bins
        num_bins = aggregated_binned_matrix.shape[1]
        bin_size_ns = 0.2  # Example bin size; adjust as needed
        time_bins = np.arange(num_bins) * bin_size_ns
        num_ticks = min(6, num_bins)  # Limit the number of ticks for readability
        tick_positions = np.linspace(0, num_bins - 1, num_ticks, dtype=int)
        tick_labels = [f"{int(time_bins[pos])}" for pos in tick_positions]  # Convert to integers
        
        ax_exist.set_xticks(tick_positions + 0.5)
        ax_exist.set_xticklabels(tick_labels, rotation=0, ha='right')
        
        # Set y-ticks with oxygen labels and ensure horizontal orientation
        ax_exist.set_yticks(np.arange(len(all_oxygens)) + 0.5)
        ax_exist.set_yticklabels(all_oxygens, rotation=0)  # Set rotation=0 for horizontal labels

        for spine in ax_exist.spines.values():
            spine.set_visible(True)
    

    
    # Save the plot
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Summary plot with existence maps and occurrence plots saved as '{output_file}'")

def main():
    # Determine the script directory and project path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    try:
        project_path = find_project_root(current_dir=script_dir)
    except FileNotFoundError as e:
        print(e)
        sys.exit(1)

    # Define the base process directory
    process_dir = os.path.join(project_path, "process.nobackup")

    # Verify that process_dir exists
    if not os.path.isdir(process_dir):
        print(f"Error: '{process_dir}' is not a valid directory.")
        sys.exit(1)

    # Initialize dictionaries to hold oxygen occurrence and existence data
    occurrence_data = {}
    existence_data = {}

    # Collect all unique oxygen atoms across all folders
    all_oxygens = set()

    # First pass: Collect all unique oxygen atoms and sort folders
    folder_list = []
    for folder_name in sorted(os.listdir(process_dir)):
        if not folder_name.startswith("IP_"):
            continue  # Skip folders that do not start with 'IP_'

        folder_path = os.path.join(process_dir, folder_name)

        # Check if it's a directory
        if not os.path.isdir(folder_path):
            continue

        # Extract species ID and count the number of '1's
        species_id = folder_name.replace("IP_", "")
        num_ones = species_id.count('1')

        # Only consider folders with num_ones in [3, 4, 5]
        if num_ones not in [3, 4, 5]:
            continue

        folder_list.append((folder_name, species_id, num_ones))

    # Sort folders first by num_ones (descending), then alphabetically by species_id
    folder_list.sort(key=lambda x: (-x[2], x[1]))

    # Extract sorted folder names and species IDs
    sorted_folders = [x[0] for x in folder_list]
    sorted_species_ids = [x[1] for x in folder_list]

    # First pass: Collect all unique oxygen atoms
    for folder_name, species_id, num_ones in folder_list:
        folder_path = os.path.join(process_dir, folder_name)

        # Paths to required files
        intra_log_file = os.path.join(folder_path, 'analyze_final_sim', 'h_bonds', 'intra_hb.log')

        # Check if the necessary files exist
        if not os.path.isfile(intra_log_file):
            print(f"Warning: '{intra_log_file}' does not exist. Skipping folder '{folder_name}'.")
            continue

        # Parse the hydrogen bond log to get donor-acceptor pairs
        hbond_df = parse_hbond_log_to_dataframe(intra_log_file)
        if hbond_df.empty:
            print(f"Warning: No valid hydrogen bond data found in '{intra_log_file}'. Skipping folder '{folder_name}'.")
            continue

        # Collect oxygen atoms from donors and acceptors
        all_oxygens.update(hbond_df['donor'].unique())
        all_oxygens.update(hbond_df['acceptor'].unique())

    # Sort all unique oxygen atoms numerically
    all_oxygens = sorted(all_oxygens, key=lambda x: int(re.findall(r'\d+', x)[0]) if re.findall(r'\d+', x) else 0)

    # Debugging: Print the number and list of oxygen atoms
    print(f"Total unique oxygen atoms collected: {len(all_oxygens)}")
    print(f"Oxygen atoms: {all_oxygens}")

    # Second pass: Process each folder in the sorted order
    for folder_name, species_id, num_ones in folder_list:
        folder_path = os.path.join(process_dir, folder_name)

        # Paths to required files
        intra_matrix_xpm = os.path.join(folder_path, 'analyze_final_sim', 'h_bonds', 'intra_hb_matrix.xpm')
        intra_log_file = os.path.join(folder_path, 'analyze_final_sim', 'h_bonds', 'intra_hb.log')

        # Check if the necessary files exist
        if not os.path.isfile(intra_matrix_xpm):
            print(f"Warning: '{intra_matrix_xpm}' does not exist. Skipping folder '{folder_name}'.")
            continue
        if not os.path.isfile(intra_log_file):
            print(f"Warning: '{intra_log_file}' does not exist. Skipping folder '{folder_name}'.")
            continue

        # Parse the hydrogen bond log to get donor-acceptor pairs
        hbond_df = parse_hbond_log_to_dataframe(intra_log_file)
        if hbond_df.empty:
            print(f"Warning: No valid hydrogen bond data found in '{intra_log_file}'. Skipping folder '{folder_name}'.")
            continue

        # Parse the XPM file to get the data matrix and metadata
        data_matrix, metadata = parse_xpm(intra_matrix_xpm)

        # Analyze the hydrogen bonds to get hbonds_per_index
        analysis_results = analyze_hydrogen_bonds(data_matrix, metadata)
        hbonds_per_index = analysis_results['hbonds_per_index']
        print(f"Folder '{folder_name}' hbonds_per_index:")
        print(hbonds_per_index)

        # Count oxygen occurrences based on the matrix data
        counts = count_oxygen_occurrences_from_matrix(hbond_df, hbonds_per_index)
        occurrence_data[species_id] = counts

        # Determine existence data (aggregated_binned_matrix)
        # Define bin_size_ns and time_per_frame_ns
        bin_size_ns = 0.2  # Adjust as needed
        time_per_frame_ns = 0.01  # As per function parameters

        frames_per_bin = int(bin_size_ns / time_per_frame_ns)
        if frames_per_bin <= 0:
            print(f"Warning: Invalid bin size for folder '{folder_name}'. Skipping existence map.")
            continue

        # Number of bins
        num_bins = data_matrix.shape[1] // frames_per_bin
        if data_matrix.shape[1] % frames_per_bin != 0:
            num_bins += 1  # Include partial bin

        # Initialize binned data matrix
        binned_matrix = np.zeros((data_matrix.shape[0], num_bins))

        # Aggregate data into bins
        for i in range(num_bins):
            start_idx = i * frames_per_bin
            end_idx = min((i + 1) * frames_per_bin, data_matrix.shape[1])
            bin_data = data_matrix[:, start_idx:end_idx]
            # Calculate the fraction of time the bond is present in the bin
            binned_matrix[:, i] = np.mean(bin_data, axis=1)

        # Aggregate data per oxygen atom for the existence map
        # Create mapping from oxygen to indices
        oxygen_to_indices = {}
        for idx, row in hbond_df.iterrows():
            donor = row['donor']
            acceptor = row['acceptor']
            if donor not in oxygen_to_indices:
                oxygen_to_indices[donor] = []
            if acceptor not in oxygen_to_indices:
                oxygen_to_indices[acceptor] = []
            oxygen_to_indices[donor].append(idx)
            oxygen_to_indices[acceptor].append(idx)

        # Remove duplicate indices per oxygen atom
        for oxygen in oxygen_to_indices:
            oxygen_to_indices[oxygen] = list(set(oxygen_to_indices[oxygen]))

        # Create aggregated_binned_matrix: rows=all unique oxygen atoms, columns=bins
        aggregated_binned_matrix = np.zeros((len(all_oxygens), num_bins))
        for o_idx, oxygen in enumerate(all_oxygens):
            bond_indices = oxygen_to_indices.get(oxygen, [])
            if bond_indices:
                # For each bin, check if any of the bonds for this oxygen are present
                for bin_idx in range(num_bins):
                    # Check if any bond involving this oxygen is present in this bin
                    bin_values = binned_matrix[bond_indices, bin_idx]
                    aggregated_binned_matrix[o_idx, bin_idx] = np.max(bin_values)
            else:
                # If no bonds involve this oxygen, keep it as 0
                pass

        # Verify that the matrix has the correct shape
        if aggregated_binned_matrix.shape[0] != len(all_oxygens):
            print(f"Error: Aggregated binned matrix shape mismatch for species '{species_id}'.")
            print(f"Expected rows: {len(all_oxygens)}, Got: {aggregated_binned_matrix.shape[0]}")
            continue

        # Add to existence_data
        existence_data[species_id] = aggregated_binned_matrix

        print(f"Processed folder '{folder_name}' (Species ID: {species_id}):")
        print(f"Occurrences: {counts.to_dict()}")

    # Check if any data was collected
    if not occurrence_data:
        print("No oxygen occurrence data collected. Exiting.")
        sys.exit(1)

    # Generate the summary plot with existence maps and occurrence plots
    generate_summary_plot(
        occurrence_data, 
        existence_data, 
        time_per_frame_ns=time_per_frame_ns, 
        output_file='oxygen_occurrences_summary.pdf',
        folder_order=sorted_species_ids  # Pass the sorted order here
    )


if __name__ == "__main__":
    main()