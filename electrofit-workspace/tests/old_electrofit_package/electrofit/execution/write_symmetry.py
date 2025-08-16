import os
import sys

def main():
    """
    Main function to handle command-line arguments and initiate processing.
    """
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python write_symmetry.py <input_file> <output_file>")
        print("Example: python write_symmetry.py ANTECHAMBER_RESP1_MOD.IN symmetry_resp_MOD.txt")
        sys.exit(1)

    # Assign command-line arguments to variables
    input_filename = sys.argv[1]
    output_filename = sys.argv[2]

    # Verify that the input file exists
    if not os.path.isfile(input_filename):
        print(f"Error: Input file '{input_filename}' does not exist.")
        sys.exit(1)

    # Mapping from atomic number to element symbol
    atomic_number_mapping = {
        1: 'H',
        6: 'C',
        8: 'O',
        15: 'P',
        # Add other atomic numbers if necessary
    }

    # Read data from the input file
    try:
        with open(input_filename, 'r') as f:
            lines = f.readlines()
    except Exception as e:
        print(f"Error: Failed to read input file '{input_filename}'.")
        print(f"Details: {e}")
        sys.exit(1)

    # Find the line where the data starts (after the second "Resp charges for organic molecule")
    start_index = None
    occurrences = 0  # Counter for occurrences of the target line
    for i, line in enumerate(lines):
        if line.strip() == "Resp charges for organic molecule":
            occurrences += 1
            if occurrences == 2:
                start_index = i + 2  # Data starts after this line
                break

    if start_index is None:
        print("Error: Data section not found in the input file.")
        sys.exit(1)

    # Initialize variables
    atoms = []
    element_counts = {}
    line_to_atom = {}

    # Process each line of data
    for idx, line in enumerate(lines[start_index:], start=1):  # idx starts from 1
        parts = line.strip().split()
        if len(parts) != 2:
            continue  # Skip lines that don't have exactly two columns
        try:
            atomic_number = int(parts[0])
            index = int(parts[1])
        except ValueError:
            print(f"Warning: Skipping line with non-integer values: '{line.strip()}'")
            continue

        element = atomic_number_mapping.get(atomic_number, f"Element{atomic_number}")
        count = element_counts.get(element, 0) + 1
        element_counts[element] = count

        atom_label = f"{element}{count}"
        atom_info = {
            'label': atom_label,
            'atomic_number': atomic_number,
            'index': index,
            'line_number': idx  # Line number starting from 1
        }
        atoms.append(atom_info)
        line_to_atom[idx] = atom_info  # line_to_atom[1] is the first atom

    # Build the output table
    output_lines = []
    for atom in atoms:
        label = atom['label']
        index = atom['index']
        line_number = atom['line_number']

        if index == 0:
            # Unique atom, no equal charge assignment
            output_line = f"{atom['atomic_number']:>3} {index:>4} \t {label}"
        else:
            # Find the atom this one is equal to
            target_atom = line_to_atom.get(index)
            if target_atom:
                target_label = target_atom['label']
                output_line = f"{atom['atomic_number']:>3} {index:>4} \t {label} = {target_label}"
            else:
                output_line = f"{atom['atomic_number']:>3} {index:>4} \t {label} (reference atom not found)"

        output_lines.append(output_line)

    # Write the output to a file
    try:
        with open(output_filename, 'w') as f:
            f.write('\n'.join(output_lines))
    except Exception as e:
        print(f"Error: Failed to write to output file '{output_filename}'.")
        print(f"Details: {e}")
        sys.exit(1)

    print(f"Processing complete. Output written to '{output_filename}'.")

if __name__ == "__main__":
    main()
    