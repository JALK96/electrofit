"""CLI for generating MOL2 atom equivalence groups (symmetry constraints).

Usage:
  electrofit-create-equiv-groups input.mol2 output.json

Parses the MOL2 file, computes equivalence groups and writes them as JSON. Exits
0 even if empty (empty JSON represents absence of symmetry constraints).
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from electrofit.domain.symmetry.equiv_groups import (
    generate_equivalence_groups,
    write_equivalence_groups,
)

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        "electrofit-create-equiv-groups",
        description="Generate MOL2 atom equivalence groups (symmetry constraints) as JSON",
    )
    p.add_argument("mol2", help="Input MOL2 file")
    p.add_argument("output", help="Output JSON file for equivalence groups")
    p.add_argument("--verbose", "-v", action="count", default=0, help="Increase verbosity (-v, -vv)")
    return p


def configure_logging(levels: int):
    if logging.getLogger().handlers:
        return
    if levels >= 2:
        lvl = logging.DEBUG
    elif levels == 1:
        lvl = logging.INFO
    else:
        lvl = logging.WARNING
    logging.basicConfig(level=lvl, format="%(message)s")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    configure_logging(args.verbose)

    mol2_path = Path(args.mol2)
    if not mol2_path.is_file():
        logging.error("Input MOL2 file not found: %s", mol2_path)
        return 1
    groups = generate_equivalence_groups(mol2_path)
    if not groups:
        logging.warning("No equivalence groups found (empty result written anyway).")
    write_equivalence_groups(groups, args.output)
    logging.info("Wrote %d groups -> %s", len(groups), args.output)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
