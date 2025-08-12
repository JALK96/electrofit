# ── prerequisites ──────────────────────────────────────────────────────────────
# conda install -c conda-forge rdkit
# -----------------------------------------------------------------------------
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

# ── configuration ─────────────────────────────────────────────────────────────
mol2_file = "/home/johannal96/MA/electrofit/process/IP_010101/run_gau_create_gmx_in/IP_010101.mol2"

# charges for 010101
your_charge_list = [
    1.3793,
    1.2067,
    1.3793,
    1.2338,
    1.4762,
    1.2338,
    0.1963,
    -0.1353,
    0.1963,
    -0.0179,
    0.3571,
    -0.0179,
    -0.4746,
    -0.2896,
    -0.4746,
    -0.4063,
    -0.5115,
    -0.4063,
    -1.0151,
    -1.0151,
    -1.0151,
    -0.8687,
    -0.8017,
    -0.8687,
    -1.0151,
    -1.0151,
    -1.0151,
    -0.8843,
    -0.7802,
    -0.8843,
    -1.0511,
    -1.0511,
    -1.0511,
    -0.8843,
    -0.7802,
    -0.8843,
    0.0805,
    0.2086,
    0.0805,
    0.1409,
    0.1392,
    0.1409,
    0.3832,
    0.3832,
    0.3988,
]

# ── molecule loading & sanity checks ──────────────────────────────────────────
mol = Chem.MolFromMol2File(mol2_file, removeHs=False)
if mol is None:
    raise ValueError(f"Failed to load molecule from {mol2_file}")
n_atoms = mol.GetNumAtoms()
if n_atoms != len(your_charge_list):
    raise ValueError(
        f"Charge list has {len(your_charge_list)} values, "
        f"but molecule has {n_atoms} atoms"
    )

print(f"Molecule loaded ({n_atoms} atoms)")

# ── generate 2-D coordinates ─────────────────────────────────────────────────
rdDepictor.SetPreferCoordGen(True)
rdDepictor.Compute2DCoords(mol)

# ── build element-counter labels (C1, C2 …) ───────────────────────────────────
elem_counters = {}
custom_labels = {}
for atom in mol.GetAtoms():
    sym = atom.GetSymbol()
    count = elem_counters.get(sym, 0) + 1
    elem_counters[sym] = count
    custom_labels[atom.GetIdx()] = f"{sym}{count}"

# ── attach charges as atom notes (displayed automatically) ────────────────────
for idx, charge in enumerate(your_charge_list):
    atom = mol.GetAtomWithIdx(idx)
    atom.SetProp("atomNote", f"{charge:+.2f}")  # “+” always shows sign

# ── drawing ───────────────────────────────────────────────────────────────────
drawer = rdMolDraw2D.MolDraw2DSVG(500, 500)
opts = drawer.drawOptions()

# draw every element in black, keep note text the same colour
# opts.updateAtomPalette({k: (0, 0, 0) for k in DrawingOptions.elemDict.keys()})

opts.addAtomIndices = False  # we make our own labels
opts.addBondIndices = False
opts.baseFontSize = 0.3  # tweak to taste
opts.annotationFontScale = 0.8  # charge text a bit smaller than label


# plug in the element-counter labels
for idx, lab in custom_labels.items():
    opts.atomLabels[idx] = lab

# prepare & draw
rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
drawer.FinishDrawing()
svg = drawer.GetDrawingText()

with open("molecule_with_charges.svg", "w") as fh:
    fh.write(svg)
print("SVG written to molecule_with_charges.svg")
