import os, sys, subprocess
from pathlib import Path
from tests.helpers.project import make_project_tree, install_minimal_gmx_fixture

def test_step3_header_and_snapshot_precedence(tmp_path, shim_bin, monkeypatch):
    monkeypatch.setenv("PATH", str(shim_bin)+os.pathsep+os.environ['PATH'])
    proj = tmp_path
    make_project_tree(proj, name="IP_011101")
    # Introduce differing forcefield values to test precedence (molecule-level should win)
    (proj/"electrofit.toml").write_text("[simulation]\nforcefield='amber14sb.ff'\n")
    mol_cfg = proj/"data"/"input"/"IP_011101"/"electrofit.toml"
    mol_cfg.write_text("[simulation]\nforcefield='amberXX.ff'\n")
    # Step1 + Step2
    for step in ("step1","step2"):
        r = subprocess.run([sys.executable,"-m","electrofit",step,"--project",str(proj)], cwd=proj, text=True, capture_output=True)
        assert r.returncode == 0, r.stderr + r.stdout
    # Fixture for step3 run
    run_dir = proj/"process"/"IP_011101"/"run_gmx_simulation"
    install_minimal_gmx_fixture(run_dir)
    # Step3
    r3 = subprocess.run([sys.executable,"-m","electrofit","step3","--project",str(proj)], cwd=proj, text=True, capture_output=True)
    assert r3.returncode == 0, r3.stderr + r3.stdout
    plog = (run_dir/"process.log").read_text().splitlines()
    assert any(line.startswith("electrofit ") and "step=step3" in line for line in plog), "Missing header in step3 process.log"
    snap = run_dir/"electrofit.toml"
    assert snap.is_file()
    assert "amberXX.ff" in snap.read_text(), "Expected molecule-level forcefield override to win precedence"
