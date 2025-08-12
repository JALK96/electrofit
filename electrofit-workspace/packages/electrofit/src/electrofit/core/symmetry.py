import os
import sys
from contextlib import contextmanager
from typing import Dict, List


@contextmanager
def _pushd(path: str | None):
    """Temporarily change the working directory if ``path`` is provided."""
    if not path:
        # no-op context manager
        yield
        return
    prev = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev)


def write_symmetry(
    input_filename: str,
    output_filename: str,
    cwd: str | None = None,
    atomic_number_mapping: Dict[int, str] | None = None,
) -> str:
    """
    Parse a RESP input file and write a compact symmetry table.

    Parameters
    ----------
    input_filename : str
        Path to the RESP input file (e.g., "ANTECHAMBER_RESP1.IN").
    output_filename : str
        Path to the output text file to write the symmetry mapping.
    cwd : Optional[str]
        If provided, the operation is executed with this directory as the
        current working directory (useful for scratch execution).
    atomic_number_mapping : Optional[Dict[int, str]]
        Mapping from atomic number to element symbol. If omitted, a default
        mapping for common elements is used.

    Returns
    -------
    str
        Absolute path to the written output file.
    """
    if atomic_number_mapping is None:
        atomic_number_mapping = {
            1: "H",
            6: "C",
            8: "O",
            15: "P",
            # Extend as needed
        }

    with _pushd(cwd):
        if not os.path.isfile(input_filename):
            raise FileNotFoundError(f"Input file '{input_filename}' does not exist.")

        # Read data from the input file
        with open(input_filename, "r") as f:
            lines: List[str] = f.readlines()

        # Find the line where the data starts (after the second occurrence)
        start_index = None
        occurrences = 0
        for i, line in enumerate(lines):
            if line.strip() == "Resp charges for organic molecule":
                occurrences += 1
                if occurrences == 2:
                    start_index = i + 2  # Data starts after this line
                    break
        if start_index is None:
            raise ValueError("Data section not found in the input file.")

        # Process data lines
        atoms = []
        element_counts: Dict[str, int] = {}
        line_to_atom = {}

        for idx, line in enumerate(lines[start_index:], start=1):  # idx starts at 1
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            try:
                atomic_number = int(parts[0])
                index = int(parts[1])
            except ValueError:
                # Skip lines that do not match the expected two integers
                continue

            element = atomic_number_mapping.get(atomic_number, f"Element{atomic_number}")
            element_counts[element] = element_counts.get(element, 0) + 1
            atom_label = f"{element}{element_counts[element]}"
            atom_info = {
                "label": atom_label,
                "atomic_number": atomic_number,
                "index": index,
                "line_number": idx,  # Line number starting from 1
            }
            atoms.append(atom_info)
            line_to_atom[idx] = atom_info  # line_to_atom[1] is the first atom

        # Build output lines
        output_lines = []
        for atom in atoms:
            label = atom["label"]
            index = atom["index"]

            if index == 0:
                # Unique atom, no equal charge assignment
                output_line = f"{atom['atomic_number']:>3} {index:>4} \t {label}"
            else:
                # Equal charge to atom at line 'index'
                target_atom = line_to_atom.get(index)
                if target_atom:
                    target_label = target_atom["label"]
                    output_line = f"{atom['atomic_number']:>3} {index:>4} \t {label} = {target_label}"
                else:
                    output_line = (
                        f"{atom['atomic_number']:>3} {index:>4} \t {label} (reference atom not found)"
                    )
            output_lines.append(output_line)

        # Write output
        with open(output_filename, "w") as f:
            f.write("\n".join(output_lines))

        return os.path.abspath(output_filename)


def main() -> None:
    """
    CLI entry point: keep backward compatibility with the previous script-style usage.
    Usage: python write_symmetry.py <input_file> <output_file>
    """
    if len(sys.argv) != 3:
        print("Usage: python write_symmetry.py <input_file> <output_file>")
        print(
            "Example: python write_symmetry.py ANTECHAMBER_RESP1_MOD.IN symmetry_resp_MOD.txt"
        )
        sys.exit(1)

    input_filename = sys.argv[1]
    output_filename = sys.argv[2]

    try:
        out_path = write_symmetry(input_filename, output_filename)
        print(f"Processing complete. Output written to '{out_path}'.")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
