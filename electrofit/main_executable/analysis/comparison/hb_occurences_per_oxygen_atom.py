import os
import sys
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set Seaborn context for better aesthetics
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
    # This is more complex and would require tracking continuous sequences of 1s
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

def generate_summary_plot(occurrence_data, output_file='oxygen_occurrences_summary.pdf'):
    """
    Generates a summary plot with multiple horizontal bar subplots showing oxygen atom occurrences per species.
    Each subplot corresponds to an oxygen atom, containing horizontal bars for each species with consistent coloring.
    Species tick labels (y-axis labels) are displayed only for subplots in the first column.
    
    Parameters
    ----------
    occurrence_data : dict
        Dictionary where keys are folder names (species) and values are pd.Series of oxygen counts.
    output_file : str
        Filename for the saved summary plot.
    
    Returns
    -------
    None
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np
    import math

    # Step 1: Transform data to have oxygen atoms as keys and species as values
    oxygen_species_data = {}
    species_list = sorted(occurrence_data.keys())  # Sorted for consistent order

    # Collect all unique oxygen atoms
    all_oxygens = set()
    for counts in occurrence_data.values():
        all_oxygens.update(counts.index.tolist())
    all_oxygens = sorted(all_oxygens, key=lambda x: int(re.findall(r'\d+', x)[0]) if re.findall(r'\d+', x) else 0)

    # Initialize dictionary
    for oxygen in all_oxygens:
        oxygen_species_data[oxygen] = pd.Series(index=species_list, dtype=int).fillna(0)

    # Populate the dictionary
    for species, counts in occurrence_data.items():
        for oxygen in all_oxygens:
            oxygen_species_data[oxygen][species] = counts.get(oxygen, 0)

    # Step 2: Assign colors to species
    num_species = len(species_list)
    # Use a palette with enough distinct colors
    if num_species <= 10:
        palette = sns.color_palette("tab10", num_species)
    elif num_species <= 20:
        palette = sns.color_palette("tab20", num_species)
    else:
        # Generate a palette with distinct colors
        palette = sns.color_palette("hsv", num_species)
    species_colors = dict(zip(species_list, palette))

    # Step 3: Determine subplot grid size (e.g., 3x4 for 12 oxygen atoms)
    num_oxygens = len(all_oxygens)
    cols = 3  # Number of columns in the subplot grid
    rows = math.ceil(num_oxygens / cols)

    # Step 4: Create subplots
    fig, axes = plt.subplots(nrows=rows, ncols=cols, figsize=(cols * 5, rows * 4), sharex=True)
    axes = axes.flatten()  # Flatten in case of multiple rows

    for idx, (oxygen, species_counts) in enumerate(oxygen_species_data.items()):
        ax = axes[idx]
        # Plot horizontal bars
        bars = ax.barh(species_list, species_counts.values, color=[species_colors[species] for species in species_list])
        
        # Set title
        ax.set_title(f'Oxygen Atom: {oxygen}', fontsize=12, fontweight='bold')
        
        # Set x-axis label
        ax.set_xlabel('Occurrences', fontsize=10)
        
        # **Conditional Y-axis Labels: Show only for the first column**
        if idx % cols == 0:
            ax.set_ylabel('Species', fontsize=10)
        else:
            ax.set_ylabel('')  # Remove y-axis label text
            ax.set_yticklabels([])  # Hide y-axis tick labels

        # Add occurrence labels next to bars
        for bar in bars:
            width = bar.get_width()
            ax.annotate(f'{int(width)}',
                        xy=(width, bar.get_y() + bar.get_height() / 2),
                        xytext=(5, 0),  # 5 points horizontal offset
                        textcoords='offset points',
                        ha='left', va='center',
                        fontsize=8)
        
        # Enhance grid
        ax.grid(True, linestyle='--', alpha=0.5, axis='x')

    # Step 5: Remove any unused subplots
    for j in range(idx + 1, len(axes)):
        fig.delaxes(axes[j])

    # Step 6: Create a single legend
    handles = [plt.Rectangle((0,0),1,1, color=species_colors[species]) for species in species_list]
    fig.legend(handles, species_list, title='Species', loc='upper right', bbox_to_anchor=(0.95, 0.95))

    # Step 7: Adjust layout to accommodate the legend
    plt.tight_layout(rect=[0, 0, 0.9, 1])

    # Save the plot
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Summary plot saved as '{output_file}'")

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

    # Initialize a dictionary to hold oxygen occurrence data
    occurrence_data = {}

    # Iterate through each subdirectory in the process directory
    for folder_name in sorted(os.listdir(process_dir)):
        if not folder_name.startswith("IP_"):
            continue  # Skip folders that do not start with 'IP_'

        folder_path = os.path.join(process_dir, folder_name)

        # Check if it's a directory
        if not os.path.isdir(folder_path):
            continue

        # Paths to required files
        intra_matrix_xpm = os.path.join(folder_path, 'analyze_final_sim', 'h_bonds', 'intra_hb_matrix.xpm')
        intra_log_file = os.path.join(folder_path, 'analyze_final_sim', 'h_bonds', 'intra_hb.log')
        intra_csv_output = os.path.join(folder_path, 'analyze_final_sim', 'h_bonds', 'intra_hbond_pairs.csv')

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
        print(hbonds_per_index)

        # Count oxygen occurrences based on the matrix data
        counts = count_oxygen_occurrences_from_matrix(hbond_df, hbonds_per_index)
        occurrence_data[folder_name] = counts

        print(f"Processed folder '{folder_name}':")
        print(counts.to_dict())

    # Check if any data was collected
    if not occurrence_data:
        print("No oxygen occurrence data collected. Exiting.")
        sys.exit(1)

    # Generate the summary plot
    generate_summary_plot(occurrence_data, output_file='oxygen_occurrences_summary.pdf')

if __name__ == "__main__":
    main()