#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import argparse
import re

# Run this script from the WORKSPACE ROOT: electrofit/reconstruct_ef/
ROOTS = [
    Path("packages/electrofit/src"),
    Path("packages/electrofit-fep/src"),
    Path("packages/electrofit-analysis/src"),
]

# 1) Exact module path renames (handles both `from X import` and `import X`)
RENAMES = {
    # --- legacy bare-package forms ---
    # removed: curlyBrace legacy (now electrofit.viz.curly_brace.draw_curly_brace)
    "electrofit.helper.config_parser": "electrofit.config_parser",
    "helper.config_parser": "electrofit.config_parser",
    "helper.file_manipulation": "electrofit.io.files",
    "helper.setup_finalize_scratch": "electrofit.scratch.manager",
    # removed: old plotting helpers (now electrofit.viz.helpers)
    "helper.set_logging": "electrofit.logging",
    "helper.eqFEP": "electrofit_fep.core.eqfep",
    "helper.neqFEP": "electrofit_fep.core.neqfep",
    "execution.update_mol2": "electrofit.io.mol2",
    "execution.edit_resp": "electrofit.io.resp",
    "execution.write_symmetry": "electrofit.core.symmetry",
    "execution.create_equiv_groups": "electrofit.domain.symmetry.equiv_groups",
    "execution.run_eqFEP": "electrofit_fep.workflows.run_eqfep",
    "execution.run_neqFEP": "electrofit_fep.workflows.run_neqfep",
    "execution.run_FEP_old": "electrofit_fep.workflows.run_fep_old",
    # treat analysis as a prefix
    "execution.analysis.": "electrofit_analysis.fep.",
    "main.gmx_simulation": "electrofit.adapters.gromacs",
    "main.process_initial_structure": "electrofit.core.process_initial_structure",
    "main.process_conform": "electrofit.core.process_conform",
    "commands.run_commands": "electrofit.cli.run_commands",
    "main_executable.analysis.dihedral_density": "electrofit_analysis.structure.dihedral_density",
    "main_executable.analysis.h_bonds": "electrofit_analysis.structure.h_bonds",
    "main_executable.analysis.NA_P_distance": "electrofit_analysis.structure.na_p_distance",
    "main_executable.analysis.NA_P_count": "electrofit_analysis.structure.na_p_count",
    "main_executable.analysis.nojump_center": "electrofit_analysis.structure.nojump_center",
    "main_executable.analysis.plot_charges": "electrofit_analysis.charges.plot_charges",
    "main_executable.analysis.plot_molecule_2d": "electrofit_analysis.structure.plot_molecule_2d",
    "main_executable.analysis.structure_2d_charge_annotation": "electrofit_analysis.structure.structure_2d_charge_annotation",

    # --- names already qualified with `electrofit.` ---
    "electrofit.helper.file_manipulation": "electrofit.io.files",
    "electrofit.helper.setup_finalize_scratch": "electrofit.scratch.manager",
    # removed: old plotting helpers (now electrofit.viz.helpers)
    "electrofit.helper.set_logging": "electrofit.logging",
    "electrofit.helper.eqFEP": "electrofit_fep.core.eqfep",
    "electrofit.helper.neqFEP": "electrofit_fep.core.neqfep",
    "electrofit.execution.update_mol2": "electrofit.io.mol2",
    "electrofit.execution.edit_resp": "electrofit.io.resp",
    "electrofit.execution.write_symmetry": "electrofit.core.symmetry",
    "electrofit.execution.create_equiv_groups": "electrofit.domain.symmetry.equiv_groups",
    "electrofit.execution.run_eqFEP": "electrofit_fep.workflows.run_eqfep",
    "electrofit.execution.run_neqFEP": "electrofit_fep.workflows.run_neqfep",
    "electrofit.execution.run_FEP_old": "electrofit_fep.workflows.run_fep_old",
    # prefix tree
    "electrofit.execution.analysis.": "electrofit_analysis.fep.",
    "electrofit.main.gmx_simulation": "electrofit.adapters.gromacs",
    "electrofit.main.process_initial_structure": "electrofit.core.process_initial_structure",
    "electrofit.main.process_conform": "electrofit.core.process_conform",
    "electrofit.commands.run_commands": "electrofit.cli.run_commands",
    "electrofit.main_executable.analysis.dihedral_density": "electrofit_analysis.structure.dihedral_density",
    "electrofit.main_executable.analysis.h_bonds": "electrofit_analysis.structure.h_bonds",
    "electrofit.main_executable.analysis.NA_P_distance": "electrofit_analysis.structure.na_p_distance",
    "electrofit.main_executable.analysis.NA_P_count": "electrofit_analysis.structure.na_p_count",
    "electrofit.main_executable.analysis.nojump_center": "electrofit_analysis.structure.nojump_center",
    "electrofit.main_executable.analysis.plot_charges": "electrofit_analysis.charges.plot_charges",
    "electrofit.main_executable.analysis.plot_molecule_2d": "electrofit_analysis.structure.plot_molecule_2d",
    "electrofit.main_executable.analysis.structure_2d_charge_annotation": "electrofit_analysis.structure.structure_2d_charge_annotation",
}

# 2) Handle `from helper import X, Y as Z` and `from electrofit.helper import ...`
HELPER_MEMBER_MAP = {
    # removed: curlyBrace legacy symbol
    "config_parser": "electrofit.config_parser",
    "file_manipulation": "electrofit.io.files",
    "setup_finalize_scratch": "electrofit.scratch.manager",
    # removed: plotting legacy symbol
    "set_logging": "electrofit.logging",
    "eqFEP": "electrofit_fep.core.eqfep",
    "neqFEP": "electrofit_fep.core.neqfep",
}

# Build regex replacements for exact renames
PAIR_PATTERNS: list[tuple[re.Pattern, str]] = []
for old, new in RENAMES.items():
    if old.endswith('.'):
        # prefix match (e.g., execution.analysis. or electrofit.execution.analysis.)
        p_from = re.compile(rf"\bfrom\s+{re.escape(old)}", flags=re.MULTILINE)
        r_from = f"from {new}"
        p_imp = re.compile(rf"\bimport\s+{re.escape(old)}", flags=re.MULTILINE)
        r_imp = f"import {new}"
    else:
        p_from = re.compile(rf"\bfrom\s+{re.escape(old)}\b", flags=re.MULTILINE)
        r_from = f"from {new}"
        p_imp = re.compile(rf"\bimport\s+{re.escape(old)}\b", flags=re.MULTILINE)
        r_imp = f"import {new}"
    PAIR_PATTERNS.append((p_from, r_from))
    PAIR_PATTERNS.append((p_imp, r_imp))

# `from helper import X, Y as Z`
FROM_HELPER_RE = re.compile(r"^\s*from\s+helper\s+import\s+([^#\n]+)", flags=re.MULTILINE)
# `from electrofit.helper import X, Y as Z`
FROM_ELECTROFIT_HELPER_RE = re.compile(r"^\s*from\s+electrofit\.helper\s+import\s+([^#\n]+)", flags=re.MULTILINE)

# `import helper.X [as Y]`
IMPORT_HELPER_MEMBER_RE = re.compile(r"^\s*import\s+helper\.([A-Za-z0-9_]+)(?:\s+as\s+([A-Za-z0-9_]+))?\s*$", flags=re.MULTILINE)
# `import electrofit.helper.X [as Y]`
IMPORT_ELECTROFIT_HELPER_MEMBER_RE = re.compile(r"^\s*import\s+electrofit\.helper\.([A-Za-z0-9_]+)(?:\s+as\s+([A-Za-z0-9_]+))?\s*$", flags=re.MULTILINE)


def _from_import_as(mod_fq: str, as_name: str | None) -> str:
    parts = mod_fq.split('.')
    parent, leaf = '.'.join(parts[:-1]), parts[-1]
    if as_name:
        return f"from {parent} import {leaf} as {as_name}"
    else:
        # bind the leaf name so downstream `leaf.func()` keeps working
        return f"from {parent} import {leaf} as {leaf}"


def _rewrite_from_helper_block(items_str: str) -> str:
    items = [s.strip() for s in items_str.split(',')]
    out_lines: list[str] = []
    for item in items:
        if not item:
            continue
        if ' as ' in item:
            name, alias = [p.strip() for p in item.split(' as ', 1)]
        else:
            name, alias = item, None
        target = HELPER_MEMBER_MAP.get(name)
        if target:
            out_lines.append(_from_import_as(target, alias or name))
        else:
            out_lines.append(f"# TODO: unmapped helper member: {item}\nfrom helper import {item}")
    return '\n'.join(out_lines)


def rewrite_text(text: str) -> str:
    new = text

    # Pass 1: exact path renames
    for pat, repl in PAIR_PATTERNS:
        new = pat.sub(repl, new)

    # Pass 2: `from helper import ...` and `from electrofit.helper import ...`
    new = FROM_HELPER_RE.sub(lambda m: _rewrite_from_helper_block(m.group(1)), new)
    new = FROM_ELECTROFIT_HELPER_RE.sub(lambda m: _rewrite_from_helper_block(m.group(1)), new)

    # Pass 3: `import helper.X [as Y]` and `import electrofit.helper.X [as Y]`
    def repl_import_member(m: re.Match) -> str:
        name = m.group(1)
        alias = m.group(2)
        target = HELPER_MEMBER_MAP.get(name)
        if target:
            return _from_import_as(target, alias or name)
        return m.group(0)

    new = IMPORT_HELPER_MEMBER_RE.sub(repl_import_member, new)
    new = IMPORT_ELECTROFIT_HELPER_MEMBER_RE.sub(repl_import_member, new)

    return new


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true", help="Print files that would change, but don't write.")
    args = ap.parse_args()

    changed = 0
    scanned = 0
    for root in ROOTS:
        for py in root.rglob("*.py"):
            scanned += 1
            t = py.read_text(encoding="utf-8")
            new = rewrite_text(t)
            if new != t:
                if args.dry_run:
                    print(f"[DRY] would rewrite: {py}")
                else:
                    py.write_text(new, encoding="utf-8")
                    print(f"rewrote {py}")
                changed += 1
    if args.dry_run:
        print(f"[DRY] {changed} files would change out of {scanned} scanned.")
    else:
        print(f"Done. {changed} files changed out of {scanned} scanned.")


if __name__ == "__main__":
    main()