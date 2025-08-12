import sys

from electrofit.io.files import (
    build_connections,
    create_equivalence_groups,
    parse_mol2,
    write_equivalence_groups,
)


def main():
    if len(sys.argv) != 3:
        print(
            "Usage: python create_equiv_groups.py <input_MOL2_file> <output_JSON_file>"
        )
        print("Example: python create_equiv_groups.py molecule.mol2 equiv_groups.json")
        sys.exit(1)

    mol2_file = sys.argv[1]
    output_json = sys.argv[2]

    atoms, bonds = parse_mol2(mol2_file)
    build_connections(atoms, bonds)
    equiv_groups = create_equivalence_groups(atoms)

    if not equiv_groups:
        print("No equivalence groups found. Please check your MOL2 file's structure.")
    else:
        write_equivalence_groups(equiv_groups, output_json)


if __name__ == "__main__":
    main()
