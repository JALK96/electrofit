# tests/e2e/test_step2_local.py
import os, sys, subprocess
from pathlib import Path
from tests.helpers.project import make_project_tree

def test_step2_builds_run_dir(tmp_path, shim_bin, monkeypatch):
    monkeypatch.setenv("PATH", str(shim_bin)+os.pathsep+os.environ["PATH"])
    proj = tmp_path
    make_project_tree(proj)

    # step1 first (creates acpype dir via shim)
    r1 = subprocess.run([sys.executable,"-m","electrofit","step1","--project",str(proj)],
                        cwd=proj, text=True, capture_output=True)
    assert r1.returncode == 0, r1.stderr + r1.stdout

    # step2 creates run_gmx_simulation with run.json, copies *_GMX.*
    r2 = subprocess.run([sys.executable,"-m","electrofit","step2","--project",str(proj)],
                        cwd=proj, text=True, capture_output=True)
    assert r2.returncode == 0, r2.stderr + r2.stdout

    run_dir = proj/"process"/"IP_011101"/"run_gmx_simulation"
    assert run_dir.is_dir()
    assert (run_dir/"run.json").is_file()
    # inputs present
    assert any(p.name.endswith(("_GMX.gro","_GMX.itp",".top")) for p in run_dir.iterdir())
    # MDP copied/linked
    assert (run_dir/"MDP").is_dir()