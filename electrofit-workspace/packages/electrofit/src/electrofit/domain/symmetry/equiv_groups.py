"""Equivalence group derivation for MOL2 structures.

Pure domain function `generate_equivalence_groups` returns a list of groups (each
group is a list of atom indices sharing equivalence under current heuristic).

Extraction from legacy `core/equiv_groups.py` (now shimmed). Command-line logic
will be provided separately via a dedicated CLI module.
"""
from __future__ import annotations

from pathlib import Path
from typing import List

from electrofit.io.files import (
    build_connections,
    create_equivalence_groups,
    parse_mol2,
    write_equivalence_groups,
)

__all__ = ["generate_equivalence_groups", "write_equivalence_groups"]


def generate_equivalence_groups(mol2_file: str | Path) -> List[list[int]]:
    """Parse a MOL2 file and compute equivalence groups.

    Parameters
    ----------
    mol2_file : str | Path
        Path to MOL2 file.

    Returns
    -------
    list[list[int]]
        List of equivalence groups (each group is a list of atom indices).
    """
    atoms, bonds = parse_mol2(str(mol2_file))
    build_connections(atoms, bonds)
    equiv_groups = create_equivalence_groups(atoms)
    return equiv_groups or []
