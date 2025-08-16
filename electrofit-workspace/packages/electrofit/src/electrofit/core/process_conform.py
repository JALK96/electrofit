import logging
import os
import subprocess
import sys
import time
import shutil

from electrofit.cli.run_commands import (
    create_gaussian_input,
    run_command,
    run_espgen,
    run_gaussian_calculation,
    run_resp,
    run_python,
)
from electrofit.io.files import (
    find_file_with_extension,
    modify_gaussian_input,
    pdb_to_mol2,
    replace_charge_in_ac_file,
)
from electrofit.cli.safe_run import ensure_finalized
from electrofit.io.resp import edit_resp_input
from electrofit.scratch.manager import (
    setup_scratch_directory,
    finalize_scratch_directory,
)
from electrofit.io.mol2 import update_mol2_charges
from electrofit.core.symmetry import write_symmetry

PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
project_path = PROJECT_PATH


def process_conform(
    molecule_name: str,
    pdb_file: str,
    base_scratch_dir: str,
    net_charge: int,
    residue_name: str,
    adjust_sym: bool = False,
    ignore_sym: bool = False,
    exit_screen: bool = True,
    protocol: str = "bcc",
):
    """Process a single conformer: PDB -> MOL2, Gaussian ESP, RESP charges.

    Steps:
      1. Convert PDB to MOL2 (standard naming/atom ordering)
      2. Build & run Gaussian single point (or opt+SP if future variant)
      3. Generate ESP surface file (`espgen`)
      4. (bcc protocol) Reconstruct RESP input via bondtype + respgen
      5. Optional symmetry constraint editing (or ignoring) for RESP1
      6. Run RESP stage 1 & 2
      7. Write RESP-charged mol2
    Scratch is auto-finalized (copy back + cleanup) via ensure_finalized.
    """

    # Decide which input files to copy to scratch
    if protocol == "opt":
        respin1_base = "ANTECHAMBER_RESP1_MOD.IN" if adjust_sym else "ANTECHAMBER_RESP1.IN"
        input_files: list[str] = [pdb_file, respin1_base, "ANTECHAMBER_RESP2.IN"]
    else:  # bcc
        if adjust_sym:
            json_file = find_file_with_extension("json")
            input_files = [pdb_file, json_file] if json_file else [pdb_file]
        else:
            input_files = [pdb_file]

    scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch_dir)
    print(f"Following {protocol} protocol.")

    try:
        defer_finalize = os.environ.get("ELECTROFIT_DEBUG_DEFER_FINALIZE") == "1"
        if defer_finalize:
            print(f"[debug-defer] start body (no auto-finalize) scratch={scratch_dir}", flush=True)
        context_mgr = None if defer_finalize else ensure_finalized(
            original_dir=original_dir,
            scratch_dir=scratch_dir,
            input_files=input_files,
            remove_parent_if_empty=False,
        )
        if context_mgr is not None:
            with context_mgr:
                # Step 0: PDB -> MOL2
                mol2_file = f"{molecule_name}.mol2"
                pdb_to_mol2(pdb_file, mol2_file, residue_name, cwd=scratch_dir)
                # Step 1
                create_gaussian_input(mol2_file, molecule_name, net_charge, scratch_dir)
                gaussian_input = f"{molecule_name}.gcrt"
                modify_gaussian_input(gaussian_input)
                # Step 2 (Gaussian) with optional debug cache skip
                gaussian_cache = os.environ.get("ELECTROFIT_DEBUG_GAUSSIAN_CACHE")
                skipped_gaussian = False
                if gaussian_cache:
                    cache_base = os.path.join(gaussian_cache, molecule_name)
                    # Expect precomputed files; copy if present
                    copied = []
                    for ext in ("gcrt", "gcrt.log", "gesp"):
                        src = f"{cache_base}.{ext}"
                        if os.path.isfile(src):
                            dst = os.path.join(scratch_dir, f"{molecule_name}.{ext}")
                            shutil.copy2(src, dst)
                            copied.append(os.path.basename(dst))
                    if any(f.endswith(".gesp") for f in copied):
                        logging.info("[debug-gaussian-cache] Using cached Gaussian artifacts: %s", copied)
                        skipped_gaussian = True
                    else:
                        logging.warning("[debug-gaussian-cache] Cache enabled but missing .gesp for %s; running Gaussian normally", molecule_name)
                if not skipped_gaussian:
                    run_gaussian_calculation(gaussian_input, molecule_name, scratch_dir)
                # Step 3
                gesp_file = f"{molecule_name}.gesp"
                esp_file = f"{molecule_name}.esp"
                run_espgen(gesp_file, esp_file, scratch_dir)
                # Step 4 (bcc)
                if protocol == "bcc":
                    run_command(f"bondtype -i {mol2_file} -o {molecule_name}.ac -f mol2 -j part", cwd=scratch_dir)
                    replace_charge_in_ac_file(file_path=f"{molecule_name}.ac", new_charge_float=net_charge, cwd=scratch_dir)
                    run_command(f"respgen -i {molecule_name}.ac -o ANTECHAMBER_RESP1.IN -f resp1", cwd=scratch_dir)
                    run_command(f"respgen -i {molecule_name}.ac -o ANTECHAMBER_RESP2.IN -f resp2", cwd=scratch_dir)
                    if adjust_sym:
                        json_filename = run_python(find_file_with_extension, "json", cwd=scratch_dir)
                        if not json_filename:
                            raise FileNotFoundError("No *.json symmetry file found in scratch.")
                        json_symmetry_file = json_filename if os.path.isabs(json_filename) else os.path.join(scratch_dir, json_filename)
                        run_python(
                            edit_resp_input,
                            input_file="ANTECHAMBER_RESP1.IN",
                            equiv_groups_file=json_symmetry_file,
                            output_file="ANTECHAMBER_RESP1_MOD.IN",
                            ignore_sym=ignore_sym,
                            cwd=scratch_dir,
                        )
                respin1_file = os.path.join(scratch_dir, "ANTECHAMBER_RESP1_MOD.IN" if adjust_sym else "ANTECHAMBER_RESP1.IN")
                respin2_file = os.path.join(scratch_dir, "ANTECHAMBER_RESP2.IN")
                for fpath in (respin1_file, respin2_file):
                    fname = os.path.basename(fpath)
                    if os.path.isfile(fpath):
                        logging.debug("[resp-check] Searching for %s: FOUND (proceeding)", fname)
                    else:
                        logging.debug("[resp-check] Searching for %s: MISSING", fname)
                        raise FileNotFoundError(f"Missing RESP input file: {fname}")
                logging.info("RESP input files ready.")
                run_python(write_symmetry, respin1_file, "symmetry_resp_MOD.txt" if adjust_sym else "symmetry_resp.txt", cwd=scratch_dir)
                resp_output1 = f"{molecule_name}-resp1.out"; resp_pch1 = f"{molecule_name}-resp1.pch"; resp_chg1 = f"{molecule_name}-resp1.chg"
                if adjust_sym:
                    logging.debug("Sleeping 1 second before first RESP (MOD file in use)")
                    time.sleep(1)
                run_resp(respin1_file, esp_file, resp_output1, resp_pch1, resp_chg1, scratch_dir)
                resp_output2 = f"{molecule_name}-resp2.out"; resp_pch2 = f"{molecule_name}-resp2.pch"; resp_chg2 = f"{molecule_name}-resp2.chg"
                run_resp(respin2_file, esp_file, resp_output2, resp_pch2, resp_chg2, scratch_dir, resp_chg1)
                mol2_resp = f"{molecule_name}_resp.mol2"
                run_python(update_mol2_charges, mol2_file, resp_chg2, mol2_resp, cwd=scratch_dir)
        else:
            # Manual mode without automatic finalize (deferred finalize)
            mol2_file = f"{molecule_name}.mol2"
            pdb_to_mol2(pdb_file, mol2_file, residue_name, cwd=scratch_dir)
            create_gaussian_input(mol2_file, molecule_name, net_charge, scratch_dir)
            gaussian_input = f"{molecule_name}.gcrt"
            modify_gaussian_input(gaussian_input)
            gaussian_cache = os.environ.get("ELECTROFIT_DEBUG_GAUSSIAN_CACHE")
            skipped_gaussian = False
            if gaussian_cache:
                cache_base = os.path.join(gaussian_cache, molecule_name)
                copied = []
                for ext in ("gcrt", "gcrt.log", "gesp"):
                    src = f"{cache_base}.{ext}"
                    if os.path.isfile(src):
                        dst = os.path.join(scratch_dir, f"{molecule_name}.{ext}")
                        shutil.copy2(src, dst)
                        copied.append(os.path.basename(dst))
                if any(f.endswith(".gesp") for f in copied):
                    logging.info("[debug-gaussian-cache] Using cached Gaussian artifacts: %s", copied)
                    skipped_gaussian = True
                else:
                    logging.warning("[debug-gaussian-cache] Cache enabled but missing .gesp for %s; running Gaussian normally", molecule_name)
            if not skipped_gaussian:
                run_gaussian_calculation(gaussian_input, molecule_name, scratch_dir)
            gesp_file = f"{molecule_name}.gesp"; esp_file = f"{molecule_name}.esp"
            run_espgen(gesp_file, esp_file, scratch_dir)
            if protocol == "bcc":
                run_command(f"bondtype -i {mol2_file} -o {molecule_name}.ac -f mol2 -j part", cwd=scratch_dir)
                replace_charge_in_ac_file(file_path=f"{molecule_name}.ac", new_charge_float=net_charge, cwd=scratch_dir)
                run_command(f"respgen -i {molecule_name}.ac -o ANTECHAMBER_RESP1.IN -f resp1", cwd=scratch_dir)
                run_command(f"respgen -i {molecule_name}.ac -o ANTECHAMBER_RESP2.IN -f resp2", cwd=scratch_dir)
                if adjust_sym:
                    json_filename = run_python(find_file_with_extension, "json", cwd=scratch_dir)
                    if not json_filename:
                        raise FileNotFoundError("No *.json symmetry file found in scratch.")
                    json_symmetry_file = json_filename if os.path.isabs(json_filename) else os.path.join(scratch_dir, json_filename)
                    run_python(
                        edit_resp_input,
                        input_file="ANTECHAMBER_RESP1.IN",
                        equiv_groups_file=json_symmetry_file,
                        output_file="ANTECHAMBER_RESP1_MOD.IN",
                        ignore_sym=ignore_sym,
                        cwd=scratch_dir,
                    )
            respin1_file = os.path.join(scratch_dir, "ANTECHAMBER_RESP1_MOD.IN" if adjust_sym else "ANTECHAMBER_RESP1.IN")
            respin2_file = os.path.join(scratch_dir, "ANTECHAMBER_RESP2.IN")
            for fpath in (respin1_file, respin2_file):
                fname = os.path.basename(fpath)
                if os.path.isfile(fpath):
                    logging.debug("[resp-check] Searching for %s: FOUND (proceeding)", fname)
                else:
                    logging.debug("[resp-check] Searching for %s: MISSING", fname)
                    raise FileNotFoundError(f"Missing RESP input file: {fname}")
            logging.info("RESP input files ready.")
            run_python(write_symmetry, respin1_file, "symmetry_resp_MOD.txt" if adjust_sym else "symmetry_resp.txt", cwd=scratch_dir)
            resp_output1 = f"{molecule_name}-resp1.out"; resp_pch1 = f"{molecule_name}-resp1.pch"; resp_chg1 = f"{molecule_name}-resp1.chg"
            if adjust_sym:
                logging.debug("Sleeping 1 second before first RESP (MOD file in use)")
                time.sleep(1)
            run_resp(respin1_file, esp_file, resp_output1, resp_pch1, resp_chg1, scratch_dir)
            resp_output2 = f"{molecule_name}-resp2.out"; resp_pch2 = f"{molecule_name}-resp2.pch"; resp_chg2 = f"{molecule_name}-resp2.chg"
            run_resp(respin2_file, esp_file, resp_output2, resp_pch2, resp_chg2, scratch_dir, resp_chg1)
            mol2_resp = f"{molecule_name}_resp.mol2"
            run_python(update_mol2_charges, mol2_file, resp_chg2, mol2_resp, cwd=scratch_dir)
            print(f"[debug-defer] manual finalize now", flush=True)
            finalize_scratch_directory(original_dir, scratch_dir, input_files, output_files=None, overwrite=True, remove_parent_if_empty=False, reason="defer-finalize")

        # Screen session exit outside context (scratch already finalized)
        print(f"[process-conform-after-finalize] {molecule_name} defer={defer_finalize}", flush=True)
        if exit_screen:
            try:
                os.chdir(original_dir)
            except Exception:
                pass
            sty = os.environ.get("STY")
            if sty:
                try:
                    subprocess.run(["screen", "-S", sty, "-X", "quit"], check=True)
                except subprocess.CalledProcessError:
                    logging.warning("Failed to quit screen session %s", sty)
            else:
                logging.debug("No STY env var; not in screen.")

    except Exception as e:
        logging.error(f"Error processing conform: {e}")
        raise
    # End marker for debugging hard exit issues
    print(f"[process-conform-end] {molecule_name}", flush=True)
