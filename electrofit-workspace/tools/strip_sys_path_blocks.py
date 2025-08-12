from __future__ import annotations
from pathlib import Path
import re

ROOT = Path("packages")
FILES = list(ROOT.rglob("*.py"))

# Remove the function+append block that was used to hack sys.path
PAT = re.compile(
    r"(?s)\n?def\s+find_project_root\([^)]*\):.*?^\s*sys\.path\.append\(project_path\)\s*\n",
    flags=re.MULTILINE,
)

removed = 0
for py in FILES:
    txt = py.read_text(encoding="utf-8")
    new = PAT.sub("\n", txt)
    if new != txt:
        py.write_text(new, encoding="utf-8")
        print("stripped", py)
        removed += 1

print(f"Done. Removed blocks from {removed} files.")