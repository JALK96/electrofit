from pathlib import Path
import os

from electrofit.io.files import copy_and_rename_folders

# Resolve project path (from env var or current working directory)
PROJECT_PATH = os.environ.get("ELECTROFIT_PROJECT_PATH", os.getcwd())
project_path = Path(PROJECT_PATH)

# Define source and destination paths
src = project_path / "data" / "input"
dst = project_path / "process"

print("Source path:", src)
print("Destination path:", dst)

# Check if src exists to prevent errors
if not src.exists():
    raise FileNotFoundError(f"Source directory '{src}' does not exist.")

# Try to locate `pis.sh` first in the project, then fall back to workspace `scripts/`
# (which lives at the workspace root alongside `packages/`).
def _find_pis_script(project_root: Path) -> Path | None:
    # 1) Project-local scripts/pis.sh
    candidate = project_root / "scripts" / "pis.sh"
    if candidate.is_file():
        print(f"Using project pis.sh at: {candidate}")
        return candidate

    # 2) Workspace fallback: <workspace>/scripts/pis.sh
    try:
        # step0_setup.py -> workflows -> electrofit -> src -> electrofit -> packages -> <workspace>
        ws_root = Path(__file__).resolve().parents[5]
        fallback = ws_root / "scripts" / "pis.sh"
        if fallback.is_file():
            print(f"Using workspace pis.sh at: {fallback}")
            return fallback
    except Exception:
        pass

    print("Warning: pis.sh not found in project or workspace; continuing without it.")
    return None

bash_src = _find_pis_script(project_path)
print("Bash script source:", bash_src)

# Build arguments for the copy helper. Only pass bash_script_source if we actually found it.
kwargs = {"source": str(src), "destination": str(dst)}
if bash_src is not None:
    kwargs["bash_script_source"] = str(bash_src)

# Execute the function
copy_and_rename_folders(**kwargs)
