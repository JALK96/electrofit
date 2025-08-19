"""Legacy helpers for historical binary->IP directory naming & mol2 renaming.

These utilities were formerly located in ``electrofit.io.files`` but are
unused by the current refactored pipeline. They are retained here temporarily
for backwards compatibility / potential one-off migration scripts.

Planned removal: after confirming no external workflows depend on them.
"""
from __future__ import annotations

from pathlib import Path
import os
import shutil
import logging

__all__ = [
    'binary_to_ip',
    'copy_and_rename_folders',
    'rename_mol2_binary',
]

# Static mapping preserved verbatim
binary_to_ip = {
    "010101": "IP21",
    "101010": "IP42",
    "101100": "IP44",
    "111000": "IP56",
    "010111": "IP23",
    "101101": "IP45",
    "111001": "IP57",
    "011111": "IP31",
    "111011": "IP59",
    "101111": "IP47",
    "111101": "IP61",
}


def copy_and_rename_folders(source: str, destination: str, nested_folder: str = "run_gau_create_gmx_in") -> None:
    """Copy directories (historical step0 helper) and create nested folder.

    No longer used in the modern code path; kept for legacy compatibility.
    """
    if not os.path.isdir(source):  # pragma: no cover - defensive legacy path
        logging.warning("[legacy.binary_mapping] source missing: %s", source)
        return
    os.makedirs(destination, exist_ok=True)
    for folder_name in os.listdir(source):
        folder_path = os.path.join(source, folder_name)
        if not os.path.isdir(folder_path):
            continue
        nested_dest_path = os.path.join(destination, folder_name, nested_folder)
        os.makedirs(nested_dest_path, exist_ok=True)
        for item in os.listdir(folder_path):
            src_item = os.path.join(folder_path, item)
            dst_item = os.path.join(nested_dest_path, item)
            if os.path.isfile(src_item):
                shutil.copy2(src_item, dst_item)
            elif os.path.isdir(src_item):
                shutil.copytree(src_item, dst_item, dirs_exist_ok=True)
        logging.info("[legacy.binary_mapping] Copied '%s' -> '%s'", folder_path, nested_dest_path)


def rename_mol2_binary(base_dir: str, binary: str) -> None:
    """Rename a *.mol2 inside ``run_gau_create_gmx_in`` folder using the binary mapping.

    Legacy one-off helper. Logs warnings instead of printing to stdout.
    """
    ip_name = binary_to_ip.get(binary)
    if ip_name is None:
        logging.warning("[legacy.binary_mapping] binary '%s' not in mapping", binary)
        return
    folder_path = os.path.join(base_dir, "run_gau_create_gmx_in")
    if not os.path.isdir(folder_path):  # pragma: no cover
        logging.warning("[legacy.binary_mapping] folder missing: %s", folder_path)
        return
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.mol2'):
            old_path = os.path.join(folder_path, file_name)
            new_path = os.path.join(folder_path, f"{ip_name}.mol2")
            try:
                os.rename(old_path, new_path)
                logging.info("[legacy.binary_mapping] Renamed '%s' -> '%s'", old_path, new_path)
            except OSError as e:  # pragma: no cover
                logging.warning("[legacy.binary_mapping] rename failed %s -> %s: %s", old_path, new_path, e)
            break
