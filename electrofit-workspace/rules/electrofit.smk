# rules/electrofit.smk
import os
from pathlib import Path

PROJECT = config["project_dir"]
RUNS    = [Path(PROJECT) / r for r in config["runs"]]
MARKERS = Path(config.get("markers_dir", ".snakemake_markers"))

def done_marker(run_dir: Path, name: str) -> str:
    rel = run_dir.relative_to(PROJECT)
    return str(MARKERS / rel / f"{name}.done")

rule all:
    input:
        # everything you want completed; start with step1
        expand(lambda r: done_marker(r, "pis"), r=RUNS),
        # later add gmx, pc, etc.
        # expand(lambda r: done_marker(r, "gmx"), r=RUNS),

rule pis:
    input:
        cfg = lambda wildcards, run_dir: str(Path(run_dir) / "electrofit.toml")
    output:
        done = lambda wildcards, run_dir: done_marker(Path(wildcards.run_dir), "pis")
    params:
        run_dir = lambda wildcards: wildcards.run_dir
    conda:
        "../envs/ambertools.yml"
    shell:
        r"""
        set -euo pipefail
        cd {PROJECT}/{params.run_dir}
        # Snakemake conda env is already active here
        electrofit run-process-initial-structure \
          --project {PROJECT} \
          --config  {input.cfg}
        mkdir -p $(dirname {output.done})
        touch {output.done}
        """

# Example: GROMACS step (replacing gmx.sh) if you want to chain it
rule gmx:
    input:
        pis_done = lambda wildcards, run_dir: done_marker(Path(run_dir), "pis")
    output:
        done    = lambda wildcards, run_dir: done_marker(Path(wildcards.run_dir), "gmx")
    params:
        run_dir = lambda wildcards: wildcards.run_dir
    conda:
        "../envs/ambertools.yml"
    shell:
        r"""
        set -euo pipefail
        cd {PROJECT}/{params.run_dir}
        python -m electrofit.external.gromacs
        mkdir -p $(dirname {output.done})
        touch {output.done}
        """