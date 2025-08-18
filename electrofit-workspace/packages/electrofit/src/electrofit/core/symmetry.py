"""Deprecated shim: moved to `electrofit.domain.symmetry.write_symmetry`.

Keeps backward compatibility for imports `from electrofit.core.symmetry import write_symmetry`.
Will be removed in a future minor release (see refactor plan Phase 4).
"""
from __future__ import annotations
import warnings
from electrofit.domain.symmetry.write_symmetry import write_symmetry  # re-export

__all__ = ["write_symmetry"]


def main() -> int:  # pragma: no cover
    warnings.warn(
        "electrofit.core.symmetry is deprecated; use electrofit.domain.symmetry.write_symmetry",
        DeprecationWarning,
        stacklevel=2,
    )
    # Keep CLI semantics if still invoked directly
    import sys
    if len(sys.argv) != 3:
        print("Usage: python -m electrofit.core.symmetry <input_file> <output_file>")
        return 1
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]
    try:
        out_path = write_symmetry(input_filename, output_filename)
        print(f"Wrote symmetry mapping to '{out_path}'.")
        return 0
    except Exception as e:  # pragma: no cover
        print(f"Error: {e}")
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
