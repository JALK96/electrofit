"""Deprecated shim for equivalence group creation.

Former location of logic now in `electrofit.domain.symmetry.equiv_groups`.
Use `generate_equivalence_groups` and `write_equivalence_groups` from the domain
module instead of invoking this script. This shim will be removed in a future
minor release.
"""
from __future__ import annotations

import sys
import warnings
from electrofit.domain.symmetry.equiv_groups import (
    generate_equivalence_groups,
    write_equivalence_groups,
)

warnings.warn(
    "equiv_groups module moved to electrofit.domain.symmetry.equiv_groups; "
    "this core shim will be removed in a future release.",
    DeprecationWarning,
    stacklevel=2,
)

def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if len(argv) != 2:
        print("Usage: python -m electrofit.core.equiv_groups <input.mol2> <output.json>")
        return 1
    mol2_path, out_json = argv
    groups = generate_equivalence_groups(mol2_path)
    if not groups:
        print("No equivalence groups found. Please check your MOL2 file's structure.")
        return 0
    write_equivalence_groups(groups, out_json)
    return 0

if __name__ == "__main__":  # pragma: no cover - legacy path
    raise SystemExit(main())
