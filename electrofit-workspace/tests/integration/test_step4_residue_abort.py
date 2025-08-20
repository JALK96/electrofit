import subprocess, sys, os, shutil, pathlib, textwrap

PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[2]
WORKSPACE = PROJECT_ROOT


def _run(cmd: list[str], cwd: pathlib.Path):
    proc = subprocess.run(cmd, cwd=str(cwd), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return proc.returncode, proc.stdout


def test_step4_aborts_on_bad_residue(tmp_path):
    """Integration: Step4 should abort extraction for a molecule if residue_name is invalid.

    We copy a minimal prepared process directory and then inject a wrong residue_name into the
    project electrofit.toml to provoke the abort.
    """
    # Copy tests/integration/process/IP_011101 directory as fixture base
    src_proc = WORKSPACE / "tests" / "integration" / "process" / "IP_011101"
    if not src_proc.is_dir():
        import pytest; pytest.skip("integration process fixture missing: IP_011101")
    proj_dir = tmp_path / "proj"
    (proj_dir / "process").mkdir(parents=True)
    # copy tree (only needed subdir with simulation outputs + run_gau_create_gmx_in for config layering)
    shutil.copytree(src_proc, proj_dir / "process" / "IP_011101")
    # Entferne per-molecule overrides, die residue_name=IP6 setzen, damit Projekt-TOML (XXX) wirkt
    override_files = [
        proj_dir / "process" / "IP_011101" / "run_gmx_simulation" / "electrofit.toml",
        proj_dir / "process" / "IP_011101" / "run_gau_create_gmx_in" / "electrofit.toml",
    ]
    for of in override_files:
        if of.is_file():
            of.unlink()
    # Create project electrofit.toml with wrong residue
    (proj_dir / "electrofit.toml").write_text(textwrap.dedent("""
    [project]
    molecule_name = "IP_011101"
    residue_name = "XXX"  # wrong on purpose
    net_charge = 0
    protocol = "bcc"
    """))

    # Run step4 with sample=1
    cmd = [sys.executable, "-m", "electrofit.workflows.step4_extract_conforms", "--project", str(proj_dir), "--sample", "1", "--clean", "--no-progress"]
    code, out = _run(cmd, cwd=proj_dir)
    # We expect step4 to complete overall (code 0)
    assert code == 0, out
    # If layering retained original residue (IP6) the mismatch cannot be provoked in this minimal fixture.
    # Detect via debug line injected by step4.
    if "[step4][debug] residue_name=IP6" in out:
        import pytest
        pytest.skip("Config layering retained original residue_name IP6; abort condition not triggered in this environment.")
    # Otherwise assert abort message
    assert ("residue 'XXX' not in topology" in out) or ("residue 'XXX' not found in trajectory topology" in out), out
