import os

from electrofit.io.files import find_file_with_extension
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
project_path = PROJECT_PATH

# Define the base process directory
process_dir = os.path.join(project_path, "data/todo")


# Loop through each subdirectory in the process directory
for folder_name in os.listdir(process_dir):
    folder_path = os.path.join(process_dir, folder_name)

    # Check if it's a directory
    if os.path.isdir(folder_path):
        # Define the 'run_final_gmx_simulation' directory within this folder
        run_final_sim_dir = os.path.join(folder_path, "run_final_gmx_simulation")
        run_gau_create_gmx_in_dir = os.path.join(folder_path, "run_gau_create_gmx_in")

        # Check if 'run_final_gmx_simulation' exists
        if os.path.isdir(run_final_sim_dir):
            # Define the destination directory 'analyze_final_sim' within the same working directory
            dest_dir = os.path.join(folder_path, "analyze_final_sim")
            os.makedirs(dest_dir, exist_ok=True)
            os.chdir(dest_dir)

            # Path to your MOL2 file
            os.chdir(run_gau_create_gmx_in_dir)
            mol2_file_name = find_file_with_extension("mol2")
            os.chdir(dest_dir)
            mol2_file = os.path.join(run_gau_create_gmx_in_dir, mol2_file_name)

            # Load the molecule
            molecule = Chem.MolFromMol2File(mol2_file, removeHs=False)

            if molecule is None:
                raise ValueError(f"Failed to load molecule from {mol2_file}")

            print("Molecule loaded successfully!")

            # Set preference to use CoordGen
            rdDepictor.SetPreferCoordGen(True)

            # Generate 2D coordinates
            rdDepictor.Compute2DCoords(molecule)

            # Create custom atom labels
            atom_counters = {}
            atom_labels = {}

            for atom in molecule.GetAtoms():
                symbol = atom.GetSymbol()
                idx = atom.GetIdx()

                # Update the counter for this atom type
                if symbol not in atom_counters:
                    atom_counters[symbol] = 1
                else:
                    atom_counters[symbol] += 1

                # Assign the custom label
                atom_labels[idx] = f"{symbol}{atom_counters[symbol]}"

            # Build a mapping from labels to atom indices
            label_to_atom_idx = {label: idx for idx, label in atom_labels.items()}

            # Set SVG size
            svg_size = 500
            # Draw the molecule with highlighted atoms
            drawer = rdMolDraw2D.MolDraw2DSVG(svg_size, svg_size)
            # Draw every atom in black
            # drawer.drawOptions().updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()})
            opts = drawer.drawOptions()
            opts.addAtomIndices = False
            opts.addBondIndices = False
            # Set font size
            opts.baseFontSize = 0.3

            # Assign custom atom labels individually
            for idx, label in atom_labels.items():
                opts.atomLabels[idx] = label

            # Prepare the molecule for drawing
            rdMolDraw2D.PrepareMolForDrawing(molecule)

            # Draw the molecule
            drawer.DrawMolecule(molecule)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            # Save the  SVG
            with open("molecule.svg", "w") as f:
                f.write(svg)
