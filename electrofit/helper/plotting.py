import numpy as np
import matplotlib.pyplot as plt
import re
import os

def plot_charges_by_atom(atoms_dict1, initial_charges_dict, base_dir, atoms_dict2=None):
    """
    Plot the charge distributions of atoms along with their average and initial charges,
    coloring the violin plots based on atom types. If a second charges dictionary is provided,
    compare the two datasets in the same plot.

    Parameters:
    - atoms_dict1: Dictionary of atoms and their charges (first dataset).
    - initial_charges_dict: Dictionary of initial charges for each atom.
    - base_dir: Directory to save the plot and charges.
    - atoms_dict2: (Optional) Dictionary of atoms and their charges (second dataset for comparison).
    """
    atom_names = list(atoms_dict1.keys())
    charges_data1 = [atoms_dict1[atom]['charges'] for atom in atom_names]
    avg_charges1 = [atoms_dict1[atom]['average_charge'] for atom in atom_names]
    init_charges = [initial_charges_dict[atom]['charges'] for atom in atom_names]
    avg_sum1 = sum(avg_charges1)

    if atoms_dict2 is not None:
        charges_data2 = [atoms_dict2[atom]['charges'] for atom in atom_names]
        avg_charges2 = [atoms_dict2[atom]['average_charge'] for atom in atom_names]
        avg_sum2 = sum(avg_charges2)
    else:
        charges_data2 = None
        avg_charges2 = None
        avg_sum2 = None

    # Define a color mapping for each atom type
    color_map = {
        'H': 'royalblue',     # Hydrogen
        'O': 'darkred',      # Oxygen
        'P': 'darkorange',   # Phosphorus
        'C': 'darkblue',     # Carbon
    }

    # Plotting
    fig, ax = plt.subplots(figsize=(12, 6))

    positions = np.arange(len(atom_names))
    width = 0.35  # Width of each violin

    positions1 = positions - width/2
    positions2 = positions + width/2

    # Plot violin plots for each atom's charge distribution (dataset 1)
    vp1 = ax.violinplot(
        charges_data1,
        positions=positions1,
        widths=width,
        showmeans=True,
        showmedians=False,
        showextrema=True,
    )

    # Plot violin plots for dataset 2 if provided
    if charges_data2 is not None:
        vp2 = ax.violinplot(
            charges_data2,
            positions=positions2,
            widths=width,
            showmeans=True,
            showmedians=False,
            showextrema=True,
        )

    # Function to set colors based on atom type
    def set_violin_colors(vp, alpha_value):
        for i, body in enumerate(vp['bodies']):
            atom_name = atom_names[i]
            # Extract the element symbol from the atom name
            match = re.match(r'^([A-Z][a-z]?)(\d*)', atom_name)
            if match:
                element = match.group(1)
            else:
                element = 'Unknown'  # Default value if no match is found
            # Get the color for this element from the color map
            color = color_map.get(element, 'black')  # Default to black if element not in color_map
            # Set the face and edge color of the violin body
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(alpha_value)

        # Customize the mean lines and other parts
        for partname in ('cmeans', 'cmins', 'cmaxes', 'cbars'):
            if partname in vp:
                items = vp[partname]
                if isinstance(items, list):
                    # For older versions of Matplotlib where items are lists
                    for i, item in enumerate(items):
                        atom_name = atom_names[i]
                        match = re.match(r'^([A-Z][a-z]?)(\d*)', atom_name)
                        if match:
                            element = match.group(1)
                        else:
                            element = 'Unknown'
                        color = color_map.get(element, 'black')
                        item.set_edgecolor(color)
                        item.set_linewidth(1.5)
                        item.set_alpha(alpha_value)
                else:
                    # For newer versions where items is a LineCollection
                    line_colors = []
                    for i in range(len(atom_names)):
                        atom_name = atom_names[i]
                        match = re.match(r'^([A-Z][a-z]?)(\d*)', atom_name)
                        if match:
                            element = match.group(1)
                        else:
                            element = 'Unknown'
                        color = color_map.get(element, 'black')
                        line_colors.append(color)
                    items.set_color(line_colors)
                    items.set_linewidth(1.5)
                    items.set_alpha(alpha_value)

    # Set colors for dataset 1
    set_violin_colors(vp1, alpha_value=0.5)

    # Set colors for dataset 2
    if charges_data2 is not None:
        set_violin_colors(vp2, alpha_value=0.2)

    # Add scatter points for initial charges
    ax.scatter(positions1, init_charges, color="black", marker='o', label="Initial Charge", zorder=5, s=5)

    # Labeling
    ax.set_xticks(positions)
    ax.set_xticklabels(atom_names, rotation=90)
    ax.set_xlabel('Atom Names')
    ax.set_ylabel('Atomic Partial Charge (e)')
    if avg_sum2 is not None:
        ax.set_title(f'Charge Distribution for Each Atom\n(Average Sum Dataset 1: {round(avg_sum1, 4)}, Dataset 2: {round(avg_sum2, 4)})')
    else:
        ax.set_title(f'Charge Distribution for Each Atom (Average Sum: {round(avg_sum1, 4)})')

    # Create custom legend
    handles = [
        plt.Line2D([], [], color='royalblue', marker='s', linestyle='None', markersize=10, label='Hydrogen'),
        plt.Line2D([], [], color='darkred', marker='s', linestyle='None', markersize=10, label='Oxygen'),
        plt.Line2D([], [], color='darkorange', marker='s', linestyle='None', markersize=10, label='Phosphorus'),
        plt.Line2D([], [], color='darkblue', marker='s', linestyle='None', markersize=10, label='Carbon'),
        plt.Line2D([], [], color='black', marker='o', linestyle='None', markersize=5, label='Initial Charge'),
    ]

    # Add dataset labels to legend if comparing
    if charges_data2 is not None:
        handles.extend([
            plt.Line2D([], [], color='grey', alpha=0.7, marker='s', linestyle='None', markersize=10, label='Dataset 1'),
            plt.Line2D([], [], color='grey', alpha=0.4, marker='s', linestyle='None', markersize=10, label='Dataset 2'),
        ])

    ax.legend(handles=handles, title='Legend', frameon=False)

    # Save and close the plot
    plt.tight_layout()
    figure_name = "charges.pdf" if atoms_dict2 is None else "charges_comparison.pdf"
    figure_path = os.path.join(base_dir, figure_name)
    plt.savefig(figure_path)
    plt.close(fig)  # Close the figure to free up memory

    # Save average charges to a file
    charges_name = "average_charges.chg" if atoms_dict2 is None else "average_charges_comparison.chg"
    charges_path = os.path.join(base_dir, charges_name)
    with open(charges_path, "w") as output:
        for i in avg_charges1:
            output.write(str(round(i, 4)) + '\n')
        if avg_charges2 is not None:
            output.write('\n')
            for i in avg_charges2:
                output.write(str(round(i, 4)) + '\n')

def create_atom_color_mapping(atom_names, symmetry_groups):
    # List of colors to use for the groups
    group_colors_list = ['darkred', 'darkgreen', 'darkorange', 'purple', 'royalblue', 'lightcoral', 'deepskyblue', 'mediumvioletred', 'orange', 'olive', 'teal', 'dodgerblue', 'darkkhaki', 'salmon', 'firebrick', 'olivedrab', 'palevioletred']
    group_colors = {}

    # Map each group to a color
    for i, (group_representative, group_atoms) in enumerate(symmetry_groups.items()):
        color = group_colors_list[i % len(group_colors_list)]  # Cycle through colors if needed
        # Include the representative atom in the group
        group = [group_representative] + group_atoms
        for atom in group:
            group_colors[atom] = color

    # Assign 'darkblue' to atoms not in any group
    atom_to_color = {}
    for atom in atom_names:
        color = group_colors.get(atom, 'darkblue')
        atom_to_color[atom] = color

    return atom_to_color


def plot_charges_by_symmetry(atoms_dict1, initial_charges_dict, base_dir, symmetry_groups, atoms_dict2=None):
    """
    Plot the charge distributions of atoms, coloring the violin plots based on symmetry groups.
    If a second charges dictionary is provided, compare the two datasets in the same plot.

    Parameters:
    - atoms_dict1: Dictionary of atoms and their charges (first dataset).
    - initial_charges_dict: Dictionary of initial charges for each atom.
    - base_dir: Directory to save the plot and charges.
    - symmetry_groups: Dictionary mapping representative atoms to lists of equivalent atoms.
    - atoms_dict2: (Optional) Dictionary of atoms and their charges (second dataset for comparison).
    """
    atom_names = list(atoms_dict1.keys())
    charges_data1 = [atoms_dict1[atom]['charges'] for atom in atom_names]
    avg_charges1 = [atoms_dict1[atom]['average_charge'] for atom in atom_names]
    init_charges = [initial_charges_dict[atom]['charges'] for atom in atom_names]
    avg_sum1 = sum(avg_charges1)

    if atoms_dict2 is not None:
        charges_data2 = [atoms_dict2[atom]['charges'] for atom in atom_names]
        avg_charges2 = [atoms_dict2[atom]['average_charge'] for atom in atom_names]
        avg_sum2 = sum(avg_charges2)
    else:
        charges_data2 = None
        avg_charges2 = None
        avg_sum2 = None

    # Create atom-to-color mapping
    atom_to_color = create_atom_color_mapping(atom_names, symmetry_groups)

    # Plotting
    fig, ax = plt.subplots(figsize=(12, 6))

    positions = np.arange(len(atom_names))
    width = 0.35  # Width of each violin

    positions1 = positions - width/2
    positions2 = positions + width/2

    # Plot violin plots for dataset 1
    vp1 = ax.violinplot(
        charges_data1,
        positions=positions1,
        widths=width,
        showmeans=False,
        showmedians=False,
        showextrema=False,
    )

    # Plot violin plots for dataset 2 if provided
    if charges_data2 is not None:
        vp2 = ax.violinplot(
            charges_data2,
            positions=positions2,
            widths=width,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )

    # Function to set colors based on symmetry groups
    def set_violin_colors(vp, alpha_value):
        for i, body in enumerate(vp['bodies']):
            atom_name = atom_names[i]
            color = atom_to_color[atom_name]
            body.set_facecolor(color)
            body.set_edgecolor(color)
            body.set_alpha(alpha_value)

        # Customize the mean lines and other parts
        for partname in ('cmeans', 'cmins', 'cmaxes', 'cbars'):
            if partname in vp:
                items = vp[partname]
                colors = [atom_to_color[atom_names[i]] for i in range(len(atom_names))]
                items.set_color(colors)
                items.set_linewidth(1.5)
                items.set_alpha(alpha_value)

    # Set colors for dataset 1
    set_violin_colors(vp1, alpha_value=0.5)

    # Set colors for dataset 2
    if charges_data2 is not None:
        set_violin_colors(vp2, alpha_value=0.2)

    # Add scatter points for initial charges
    ax.scatter(
        positions1,
        init_charges,
        color="black",
        marker='o',
        label="Initial Charge",
        zorder=5,
        alpha=1,
        s=5
    )

    # Add scatter points for average charges
    ax.scatter(
        positions1,
        avg_charges1,
        color="black",
        marker='^',
        label="Average Charge",
        zorder=5,
        alpha=1,
        s=5
    )

    if charges_data2 is not None:
        # Add scatter points for average charges
        ax.scatter(
            positions2,
            avg_charges2,
            color="black",
            marker='^',
            label="Initial Charge",
            zorder=5,
            alpha=0.5,
            s=5
        )


    # Labeling
    ax.set_xticks(positions)
    ax.set_xticklabels(atom_names, rotation=60)
    ax.set_xlabel('Atom Names')
    ax.set_ylabel('Atomic Partial Charge (e)')
    if avg_sum2 is not None:
        ax.set_title(f'Charge Distribution for Each Atom\n(Average Sum Dataset 1: {round(avg_sum1, 2)}, Dataset 2: {round(avg_sum2, 2)})')
    else:
        ax.set_title(f'Charge Distribution for Each Atom (Average Sum: {round(avg_sum1, 2)})')

    # Create custom legend
    group_colors = list(set(atom_to_color.values()))
    handles = []
    #for color in group_colors:
    #    handles.append(plt.Line2D([], [], color=color, marker='s', linestyle='None', markersize=10, label=color))

    handles.append(plt.Line2D([], [], color='black', marker='o', linestyle='None', markersize=5, label='Initial Charge'))
    handles.append(plt.Line2D([], [], color='black', marker='^', linestyle='None', markersize=5, label='Average Charge'))

    # Add dataset labels to legend if comparing
    if charges_data2 is not None:
        handles.extend([
            plt.Line2D([], [], color='darkblue', alpha=0.5, marker='s', linestyle='None', markersize=10, label='Dataset 1'),
            plt.Line2D([], [], color='darkblue', alpha=0.2, marker='s', linestyle='None', markersize=10, label='Dataset 2'),
        ])

    ax.legend(handles=handles, frameon=False)

    # Save and close the plot
    plt.tight_layout()
    figure_name = "charges_by_symmetry.pdf" if atoms_dict2 is None else "charges_by_symmetry_comparison.pdf"
    figure_path = os.path.join(base_dir, figure_name)
    plt.savefig(figure_path)
    plt.close(fig)  # Close the figure to free up memory

    # Save average charges to a file
    charges_name = "average_charges.chg" if atoms_dict2 is None else "average_charges_comparison.chg"
    charges_path = os.path.join(base_dir, charges_name)
    with open(charges_path, "w") as output:
        for i in avg_charges1:
            output.write(str(round(i, 4)) + '\n')
        if avg_charges2 is not None:
            output.write('\n')
            for i in avg_charges2:
                output.write(str(round(i, 4)) + '\n')