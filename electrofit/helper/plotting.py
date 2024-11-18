import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import re
import os

def plot_charges_by_atom(atoms_dict, initial_charges_dict, base_dir):
    """
    Plot the charge distributions of atoms along with their average and initial charges,
    coloring the violin plots based on atom types.

    Parameters:
    - atoms_dict: Dictionary of atoms and their charges.
    - initial_charges_dict: Dictionary of initial charges for each atom.
    - base_dir: Directory to save the plot and charges.
    """
    atom_names = list(atoms_dict.keys())
    charges_data = [atoms_dict[atom]['charges'] for atom in atom_names]
    avg_charges = [atoms_dict[atom]['average_charge'] for atom in atom_names]
    init_charges = [initial_charges_dict[atom]['charges'] for atom in atom_names]
    avg_sum = sum(avg_charges)

    # Define a color mapping for each atom type
    color_map = {
        'H': 'royalblue',     # Hydrogen
        'O': 'darkred',      # Oxygen
        'P': 'darkorange',   # Phosphorus
        'C': 'darkblue',     # Carbon
        # Add more elements and colors as needed
    }

    # Plotting
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot violin plots for each atom's charge distribution
    vp = ax.violinplot(
        charges_data,
        positions=range(len(atom_names)),
        widths=0.6,
        showmeans=True,
        showmedians=False,
        showextrema=True,
    )

    # Iterate over each violin plot and set colors based on atom type
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

    # Customize the mean lines
    if 'cmeans' in vp:
        means = vp['cmeans']
        # Prepare a list of colors for the mean lines
        mean_colors = []
        for i in range(len(atom_names)):
            atom_name = atom_names[i]
            match = re.match(r'^([A-Z][a-z]?)(\d*)', atom_name)
            if match:
                element = match.group(1)
            else:
                element = 'Unknown'
            color = color_map.get(element, 'black')
            mean_colors.append(color)
        # Set the colors and linewidths for the mean lines
        vp['cmeans'].set_color(mean_colors)
        vp['cmeans'].set_linewidth(2)

    # Set the color of the vertical lines and caps
    for partname in ('cbars', 'cmins', 'cmaxes'):
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
                    if partname == 'cbars':
                        item.set_alpha(0.7)
            else:
                # For newer versions where items is a LineCollection
                # Prepare a list of colors
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

    # Add scatter points for average and initial charges
    #ax.scatter(range(len(atom_names)), avg_charges, color="black", marker='^', label="Average Charge", zorder=5)
    ax.scatter(range(len(atom_names)), init_charges, color="black", marker='o', label="Initial Charge", zorder=5, s=5)

    # Labeling
    ax.set_xticks(range(len(atom_names)))
    ax.set_xticklabels(atom_names, rotation=90)
    ax.set_xlabel('Atom Names')
    ax.set_ylabel('Atomic Partial Charge (e)')
    ax.set_title(f'Charge Distribution for Each Atom (Average Sum: {round(avg_sum, 4)})')

    # Create custom legend
    handles = [
        plt.Line2D([], [], color='royalblue', marker='None', linestyle='-', markersize=10, label='Hydrogen'),
        plt.Line2D([], [], color='darkred', marker='None', linestyle='-', markersize=10, label='Oxygen'),
        plt.Line2D([], [], color='darkorange', marker='None', linestyle='-', markersize=10, label='Phosphorus'),
        plt.Line2D([], [], color='darkblue', marker='None', linestyle='-', markersize=10, label='Carbon'),

        plt.Line2D([], [], color='black', marker='o', linestyle='None', markersize=5, label='Initial Charge'),
    ]
    ax.legend(handles=handles, title='Average Charges', frameon=False)

    # Save and close the plot
    plt.tight_layout()
    figure_path = os.path.join(base_dir, "charges.pdf")
    plt.savefig(figure_path)
    plt.close(fig)  # Close the figure to free up memory

    # Save average charges to a file
    charges_path = os.path.join(base_dir, "average_charges.chg")
    with open(charges_path, "w") as output:
        for i in avg_charges:
            output.write(str(round(i,4)) + '\n')

def create_atom_color_mapping(atom_names, symmetry_groups):
    # List of colors to use for the groups
    group_colors_list = ['darkred', 'darkgreen', 'darkorange', 'purple', 'royalblue', 'lightcoral', 'deepskyblue', 'mediumvioletred']
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

def plot_charges_by_symmetry(atoms_dict, initial_charges_dict, base_dir, symmetry_groups):
    """
    Plot the charge distributions of atoms, coloring the violin plots based on symmetry groups.

    Parameters:
    - atoms_dict: Dictionary of atoms and their charges.
    - initial_charges_dict: Dictionary of initial charges for each atom.
    - base_dir: Directory to save the plot and charges.
    - symmetry_groups: Dictionary mapping representative atoms to lists of equivalent atoms.
    """
    import os
    import matplotlib.pyplot as plt

    atom_names = list(atoms_dict.keys())
    charges_data = [atoms_dict[atom]['charges'] for atom in atom_names]
    avg_charges = [atoms_dict[atom]['average_charge'] for atom in atom_names]
    init_charges = [initial_charges_dict[atom]['charges'] for atom in atom_names]
    avg_sum = sum(avg_charges)

    # Create atom-to-color mapping
    atom_to_color = create_atom_color_mapping(atom_names, symmetry_groups)

    # Plotting
    fig, ax = plt.subplots(figsize=(12, 6))

    # Plot violin plots for each atom's charge distribution
    vp = ax.violinplot(
        charges_data,
        positions=range(len(atom_names)),
        widths=0.6,
        showmeans=True,
        showmedians=False,
        showextrema=True,
    )

    # Customize the violin plots
    for i, body in enumerate(vp['bodies']):
        atom_name = atom_names[i]
        color = atom_to_color[atom_name]
        body.set_facecolor(color)
        body.set_edgecolor(color)
        #body.set_alpha(0.7)

    # Customize the mean lines
    if 'cmeans' in vp:
        mean_colors = [atom_to_color[atom_names[i]] for i in range(len(atom_names))]
        vp['cmeans'].set_color(mean_colors)
        vp['cmeans'].set_linewidth(2)

    # Set the color of the vertical lines and min/max markers (caps)
    for partname in ('cbars', 'cmins', 'cmaxes'):
        if partname in vp:
            items = vp[partname]
            colors = [atom_to_color[atom_names[i]] for i in range(len(atom_names))]
            items.set_color(colors)
            items.set_linewidth(1.5)

    # Add red dots for the initial charge
    ax.scatter(
        range(len(atom_names)),
        init_charges,
        color="black",
        marker='o',
        label="Initial Charge",
        zorder=5,
        alpha=1,
        s=5
    )

    # Labeling
    ax.set_xticks(range(len(atom_names)))
    ax.set_xticklabels(atom_names, rotation=60)
    ax.set_xlabel('Atom Names')
    ax.set_ylabel('Atomic Partial Charge (e)')
    ax.set_title(f'Charge Distribution for Each Atom (Average Sum: {round(avg_sum, 4)})')
    ax.legend(frameon=False)

    # Save and close the plot
    plt.tight_layout()
    figure_path = os.path.join(base_dir, "charges_by_symmetry.pdf")
    plt.savefig(figure_path)
    plt.close(fig)  # Close the figure to free up memory

    # Save average charges to a file
    charges_path = os.path.join(base_dir, "average_charges.chg")
    with open(charges_path, "w") as output:
        for i in avg_charges:
            output.write(str(round(i, 4)) + '\n')


