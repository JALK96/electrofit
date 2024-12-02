import MDAnalysis as mda
from MDAnalysis.lib.distances import calc_dihedrals
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Function to load .xvg files, skipping header lines
def load_hb_num_xvg(filename):
    """
    Load data from an .xvg file, skipping lines that start with @ or #.

    Parameters:
        filename (str): Path to the .xvg file.

    Returns:
        np.ndarray: 2D array with the data.
    """
    data = []
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith(('@', '#')):
                continue  # Skip header lines
            parts = line.strip().split()
            if len(parts) >= 2:  # Ensure there are at least two columns
                try:
                    # Convert all parts to float
                    floats = [float(part) for part in parts]
                    data.append(floats)
                except ValueError:
                    # Skip lines that don't contain valid floats
                    continue
    return np.array(data)

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import re

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

def visualize_data(xpm_file='intra_hb_matrix.xpm'):
    """
    Visualizes the hydrogen bond data and analysis.
    
    Parameters:
    - file_path: str, path to the XPM file.

    Returns:
    - None
    """

    data_matrix, metadata = parse_xpm(xpm_file)
    analysis_results = analyze_hydrogen_bonds(data_matrix, metadata)

    matrix_shape = np.shape(data_matrix)
    print(np.shape(data_matrix))
    plt.figure(figsize=(8, 6))

    # Heatmap of hydrogen bonds
    plt.imshow(data_matrix, aspect='auto', cmap="Reds", origin='lower')
    plt.title(metadata.get('title', 'Hydrogen Bond Existence Map'))
    plt.xlabel(metadata.get('x-label', 'Time'))
    plt.ylabel(metadata.get('y-label', 'Hydrogen Bond Index'))
    plt.colorbar(label='Bond Presence')
    plt.yticks(np.arange(matrix_shape[0]))
    plt.tight_layout()
    plt.savefig("intra_hb_existence_map.pdf")


    # Additional Analysis: Example - Histogram of Hydrogen Bonds per Index
    plt.figure(figsize=(8, 6))
    plt.bar(range(len(analysis_results['hbonds_per_index'])), analysis_results['hbonds_per_index'], color="darkred")
    plt.title('Total Occurrence of Each Hydrogen Bond')
    plt.xlabel('Hydrogen Bond Index')
    plt.ylabel('Occurrences')
    plt.grid(False)
    plt.xticks(np.arange(matrix_shape[0]))
    plt.savefig("intra_hb_occurences.pdf")

    # Example: Lifetime Distribution
    # Flatten all lifetimes
    all_lifetimes = [lifetime for bond_lifetimes in analysis_results['lifetimes'] for lifetime in bond_lifetimes]
    plt.figure(figsize=(8, 6))
    plt.hist(all_lifetimes, bins=range(1, max(all_lifetimes)+2), edgecolor='black', color="darkred")
    plt.title('Hydrogen Bond Lifetime Distribution')
    plt.xlabel('Lifetime (number of frames)')
    plt.ylabel('Frequency')
    plt.grid(False)
    plt.savefig("intra_hb_lifetime.pdf")



def plot_hb_dist_xvg(file_path, plot_file_name="hb_distribution.pdf"):
    """
    Reads an XVG file and plots the data using matplotlib.
    
    Parameters:
    - file_path: str, path to the XVG file.
    
    Returns:
    - None
    """
    title = "XVG Plot"
    x_label = "X-axis"
    y_label = "Y-axis"
    plot_type = "line"  # default plot type
    data = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            
            # Skip empty lines
            if not line:
                continue
            
            # Handle metadata lines
            if line.startswith('#'):
                continue  # Ignore comments
            elif line.startswith('@'):
                tokens = line.split()
                if len(tokens) >= 3:
                    if tokens[1] == "title":
                        title = " ".join(tokens[2:]).strip('"')
                    elif tokens[1] == "xaxis":
                        if tokens[2] == "label":
                            x_label = " ".join(tokens[3:]).strip('"')
                    elif tokens[1] == "yaxis":
                        if tokens[2] == "label":
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
        print("No data found in the XVG file.")
        return

    # Convert data to NumPy arrays for easier handling
    data = np.array(data)
    x = data[:, 0]
    y = data[:, 1]

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, marker='o', linestyle='-', color='darkblue')
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"{plot_file_name}.pdf", format="pdf")

def find_project_root(current_dir, project_name="electrofit"):
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
script_dir = os.path.dirname(os.path.abspath(__file__))
project_path = find_project_root(current_dir=script_dir)

sys.path.append(project_path)

from electrofit.commands.run_commands import run_command

# Define the base process directory
process_dir = os.path.join(project_path, "process.nobackup")


# Loop through each subdirectory in the process directory
for folder_name in os.listdir(process_dir):
    folder_path = os.path.join(process_dir, folder_name)
    
    # Check if it's a directory
    if os.path.isdir(folder_path):
        # Define the 'run_final_gmx_simulation' directory within this folder
        run_final_sim_dir = os.path.join(folder_path, 'run_final_gmx_simulation')
        
        # Check if 'run_final_gmx_simulation' exists
        if os.path.isdir(run_final_sim_dir):
            # Define the destination directory 'analyze_final_sim' within the same working directory
            dest_dir = os.path.join(folder_path, 'analyze_final_sim')
            os.makedirs(dest_dir, exist_ok=True)
            os.chdir(dest_dir)

            gro_file_path = os.path.join(run_final_sim_dir, 'md.tpr')
            xtc_file_path = os.path.join(run_final_sim_dir, 'md_center.xtc')


            run_command(f'echo "2\n2\n" | gmx hbond -s {gro_file_path} -f {xtc_file_path} -hbn intra_hb_idx.ndx -num intra_hb_num.xvg -dist intra_hb_dist.xvg -g intra_hb.log -hbm intra_hb_matrix.xpm', cwd=dest_dir)
            run_command(f'echo "2\n5\n" | gmx hbond -s {gro_file_path} -f {xtc_file_path} -hbn inter_hb_idx.ndx -num inter_hb_num.xvg -dist inter_hb_dist.xvg -g inter_hb.log -hbm inter_hb_matrix.xpm', cwd=dest_dir)
            
            # ---- Plot intra ----
            visualize_data(xpm_file="intra_hb_matrix.xpm")
            # Load the hb_num.xvg file
            xvg_filename = 'intra_hb_num.xvg'  # Replace with your actual .xvg file name
            data = load_hb_num_xvg(xvg_filename)

            # Check if data was loaded correctly
            if data.size == 0:
                raise ValueError(f"No data found in {xvg_filename}.")

            # Extract columns
            time = data[:, 0]  # Time in ps
            num_hbonds = data[:, 1]  # Number of Hydrogen Bonds (s0)
            pairs_within_0_35_nm = data[:, 2] if data.shape[1] > 2 else None  # Pairs within 0.35 nm (s1)

            # Plotting
            plt.figure(figsize=(8, 6))
            plt.plot(time, num_hbonds, label='Hydrogen Bonds', color='darkblue', linewidth=2)
            if pairs_within_0_35_nm is not None:
                plt.plot(time, pairs_within_0_35_nm, label='Pairs within 0.35 nm', color='darkred', linewidth=2)

            plt.title('Hydrogen Bonds Over Time')
            plt.xlabel('Time (ps)')
            plt.ylabel('Number')
            plt.legend()
            plt.grid(True)

            # Save the plot as SVG
            output_svg = 'intra_hb_num_over_time.pdf'
            plt.savefig(output_svg, format='pdf')
            plt.close()

            print(f"Plot saved as {output_svg}")

            # Plot and save hb donor acceptor distance distribution:
            plot_hb_dist_xvg('intra_hb_dist.xvg', plot_file_name='intra_hb_distriution.pdf')

            # ---- Plot inter ----
            # Load the hb_num.xvg file
            xvg_filename = 'inter_hb_num.xvg'  # Replace with your actual .xvg file name
            data = load_hb_num_xvg(xvg_filename)

            # Check if data was loaded correctly
            if data.size == 0:
                raise ValueError(f"No data found in {xvg_filename}.")

            # Extract columns
            time = data[:, 0]  # Time in ps
            num_hbonds = data[:, 1]  # Number of Hydrogen Bonds (s0)
            pairs_within_0_35_nm = data[:, 2] if data.shape[1] > 2 else None  # Pairs within 0.35 nm (s1)

            # Plotting
            plt.figure(figsize=(8, 6))
            plt.plot(time, num_hbonds, label='Hydrogen Bonds', color='darkblue', linewidth=2)
            if pairs_within_0_35_nm is not None:
                plt.plot(time, pairs_within_0_35_nm, label='Pairs within 0.35 nm', color='darkred', linewidth=2)

            plt.title('Hydrogen Bonds Over Time')
            plt.xlabel('Time (ps)')
            plt.ylabel('Number')
            plt.legend()
            plt.grid(True)

            # Save the plot as SVG
            output_svg = 'inter_hb_num_over_time.pdf'
            plt.savefig(output_svg, format='pdf')
            plt.close()

            print(f"Plot saved as {output_svg}")

            # Plot and save hb donor acceptor distance distribution:
            plot_hb_dist_xvg('inter_hb_dist.xvg', plot_file_name='inter_hb_distriution.pdf')
            
            