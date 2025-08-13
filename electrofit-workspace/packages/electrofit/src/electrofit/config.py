# packages/electrofit/src/electrofit/config.py
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import typing as t

try:
    import tomllib  # py>=3.11
except ModuleNotFoundError:  # pragma: no cover
    import tomli as tomllib  # fallback if ever needed

from .config_parser import ConfigParser as LegacyConfigParser


@dataclass
class ProjectSection:
    name: str | None = None
    molecule_name: str | None = None
    residue_name: str | None = None
    charge: int | None = None
    protocol: str | None = None
    adjust_symmetry: bool = False
    ignore_symmetry: bool = False
    atom_type: str | None = None


@dataclass
class PathsSection:
    mdp_dir: str = "data/MDP"
    base_scratch_dir: str | None = None


@dataclass
class HPCSection:
    enabled: bool = True
    host: str = "qcm04"
    username: str | None = None
    conda_env: str = "AmberTools23"
    shell_init: str = "~/.bashrc"
    use_screen: bool = True
    screen_name_prefix: str = "ef"


@dataclass
class GMXSection:
    script_relpath: str = "scripts/gmx.sh"
    entrypoint: str = "python -m electrofit.external.gromacs"


@dataclass
class Config:
    project_root: Path
    project: ProjectSection = field(default_factory=ProjectSection)
    paths: PathsSection = field(default_factory=PathsSection)
    hpc: HPCSection = field(default_factory=HPCSection)
    gmx: GMXSection = field(default_factory=GMXSection)


def _load_toml(path: Path) -> dict:
    with path.open("rb") as f:
        return tomllib.load(f)


def load_config(
    project_root: t.Union[str, Path],
    config_path: t.Union[str, Path, None] = None,
) -> Config:
    root = Path(project_root).resolve()
    data: dict = {}

    cfg_path = Path(config_path) if config_path else root / "electrofit.toml"
    if cfg_path.is_file():
        data = _load_toml(cfg_path)
    else:
        # --- Fallback to legacy .ef ---
        ef_path: Path | None = None
        # look in project root and data/input convenience paths
        for cand in [root] + [root / "data" / "input"]:
            if cand.is_dir():
                ef_files = list(cand.glob("*.ef"))
                if ef_files:
                    ef_path = ef_files[0]
                    break
        if ef_path is None:
            ef_path = next(root.rglob("*.ef"), None)

        if ef_path:
            legacy = LegacyConfigParser(str(ef_path))
            data = {
                "project": {
                    "name": getattr(legacy, "MoleculeName", None),
                    "molecule_name": getattr(legacy, "MoleculeName", None),
                    "residue_name": getattr(legacy, "ResidueName", None),
                    "charge": getattr(legacy, "Charge", None),
                    "protocol": getattr(legacy, "Protocol", None),
                    "adjust_symmetry": getattr(legacy, "AdjustSymmetry", False),
                    "ignore_symmetry": getattr(legacy, "IgnoreSymmetry", False),
                    "atom_type": getattr(legacy, "AtomType", None),
                },
                "paths": {
                    "base_scratch_dir": getattr(legacy, "BaseScratchDir", None),
                },
            }

    cfg = Config(project_root=root)

    def merge(section_name: str):
        section = getattr(cfg, section_name)
        for k, v in data.get(section_name, {}).items():
            if hasattr(section, k) and v is not None:
                setattr(section, k, v)
        setattr(cfg, section_name, section)

    for sec in ("project", "paths", "hpc", "gmx"):
        merge(sec)
    return cfg


def write_hpc_env_file(cfg: Config, dest_dir: t.Union[str, Path]) -> Path:
    dest = Path(dest_dir)
    env_path = dest / ".electrofit_env"
    lines = [
        f'REMOTE_HOST={cfg.hpc.host}',
        f'REMOTE_USER={cfg.hpc.username or ""}',
        f'REMOTE_SHELL_INIT={cfg.hpc.shell_init}',
        f'REMOTE_CONDA_ENV={cfg.hpc.conda_env}',
        f'USE_SCREEN={"1" if cfg.hpc.use_screen else "0"}',
        f'SCREEN_NAME_PREFIX={cfg.hpc.screen_name_prefix}',
        f'PYTHON_ENTRYPOINT={cfg.gmx.entrypoint}',
    ]
    env_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return env_path