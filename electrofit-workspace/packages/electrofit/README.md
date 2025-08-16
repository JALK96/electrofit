# electrofit (core)

Minimal core tooling to process small molecules / residues, run short MD (GROMACS) and extract conformers for charge workflows.

## Quick Start

Create / activate an environment, install the package (editable example):

```bash
pip install -e ./packages/electrofit
```

Prepare a project folder (contains `electrofit.toml`, `data/input/<MOLNAME>/<files>`):

```bash
electrofit step0 --project /path/to/project
```

Run successive workflow steps (see below). Each step reads the project root you pass with `--project` and optional config overrides via `--config`.

## Configuration (`electrofit.toml`)

Key sections (example):

```toml
[project]
molecule_name = "IP_011101"
residue_name  = "IP6"
charge        = -8
protocol      = "bcc"     # or "opt"

[paths]
mdp_dir = "data/MDP"
base_scratch_dir = "/scratch/${USER}/electrofit"

[gromacs.runtime]
threads = 8
pin = true

[sampling]
method = "maxmin"   # linear | random | maxmin
count = 25
seed = 123
```

CLI arguments override config values. (Früherer ENV-Fallback `ELECTROFIT_CONFIG_PATH` wurde entfernt – nur noch explizit via `--config`.)

## Workflow Steps (0–4)

### step0 – Initial project scaffold
Creates baseline directory structure (e.g. `process/`), validates inputs, copies / normalises initial structural files if required.

Command:
```
electrofit step0 --project <project_root>
```

### step1 – Initial processing
Legacy pre-processing (Gaussian / RESP input staging, symmetry handling). Produces intermediate files placed under `process/<mol>/run_gau_create_gmx_in`.

### step2 – Simulation directory setup
Builds per-molecule GROMACS run directories (`run_gmx_simulation`), inserts MDP templates, generates topology / coordinate starting point.

### step3 – Start short MD production
Runs energy minimisation + equilibration + (short) production with parameterised runtime flags (threads, pinning). Produces `md_center.xtc` and `md.gro` used later.

```bash
electrofit step3 --project <project_root>
```

### step4 – Extract representative conformers
Reads each `process/<mol>/run_gmx_simulation` trajectory and writes sampled conformers into `process/<mol>/extracted_conforms/` (PDB + copied ancillary input files).

```bash
electrofit step4 --project <project_root> \
	--sample 20 --sampling-method maxmin --seed 123 --workers 2
```

## Conformer Sampling (Step 4)

Supported strategies (`--sampling-method`):

* `linear` (default): evenly spaced frames (deterministic)
* `random`: uniform random without replacement (seeded)
* `maxmin`: farthest-point (max–min) diversity in RMSD space

Important options:

* `--sample N` desired number of conformers (fallback 20 or `[sampling].count`)
* `--sampling-method METHOD` selection strategy
* `--seed S` seed for random / starting frame in maxmin
* `--workers K` parallelise across molecule directories (0 = auto heuristic)
* `--no-progress` disable the progress bar

Config section (optional):

```toml
[sampling]
method = "maxmin"
count = 25
seed = 123
```

CLI overrides config. Output line includes method, e.g. `Extracted 25 conformers (method=maxmin) ...`.

Performance notes:

* `maxmin` costs roughly O(k * N) RMSD evaluations (k = selected frames). For very large trajectories try a smaller `--sample` first.
* Parallelisation is at the granularity of molecules, not frames inside one trajectory.

Planned extensions: RMSD cutoff incremental selection, k-medoids clustering, PCA/grid stratified sampling.

## Parallel Execution

`--workers` uses a process pool. If unset / 0 an automatic value (CPU count minus one) capped by number of molecules is used. Sequential mode (workers=1) retains a live progress bar.

## Logging & Scratch

Scratch directories are created under `base_scratch_dir` (environment variables expanded). Finalisation is idempotent. GROMACS step automates interactive selections (e.g. `trjconv`).

## Testing

Unit tests include sampling selection determinism and basic functional checks. Shortened MDP templates keep test runtime low.

## License

## Konfig-Präzedenz & Fill-In

Die endgültige `electrofit.toml` Snapshot-Datei in jedem Arbeitsverzeichnis entsteht über einen klar definierten Layering-Prozess:

Starke Overrides (überschreiben vorhandene Werte in Reihenfolge – spätere gewinnt):
1. Molekül-spezifische Eingabe: `data/input/<MOL>/electrofit.toml`
2. Prozess-spezifische Datei: `process/<MOL>/electrofit.toml` bzw. vorherige Schritt-Ausgabe (z.B. `run_gau_create_gmx_in/electrofit.toml` oder `results/electrofit.toml` je nach Schritt)
3. CLI `--config` (falls angegeben)

Fill-In Ebene (füllt nur fehlende Schlüssel, überschreibt niemals):
4. Projektweite Defaults: `<project_root>/electrofit.toml`

Semantik:
* „Override“ ersetzt bestehende Werte (deep merge, scalars & Subtrees vollständig überschrieben).
* „Fill-In“ ergänzt nur Keys, die noch nicht existieren (rekursiv), lässt vorhandene Werte unverändert.
* Dieses Muster verhindert, dass projektweite Defaults unabsichtlich molekül-spezifische Einstellungen verdrängen, reduziert aber Duplikation bei globalen Parametern.

Beispiel – gewünschte Kraftfeld-Priorität:
* Molekül: `simulation.forcefield = "amber14sb.ff"`
* Projekt setzt keinen Forcefield-Schlüssel → Snapshot übernimmt Molekülwert.
* Falls Projekt später `simulation.forcefield` hätte, würde er NICHT das molekül-spezifische überschreiben (weil Projekt nur Fill-In ist) – gewünschtes Verhalten.

Logging-Markierungen:
* `[config][override] key: old -> new` für Überschreibungen durch starke Layer.
* `[config][fill] key: value` wenn ein fehlender Schlüssel durch Projekt-Defaults ergänzt wurde.
* `[config] no overrides applied ...` / `no fills` wenn keine Änderung.

CLI `--config`:
* Wird als stärkste Override-Ebene eingehängt (nach Molekül & Prozess), um gezielt Werte temporär zu ersetzen.
* Überschreibt vorhandene Werte, führt keine Fill-Ins durch.

Implementierung:
* Zentral in `build_snapshot_with_layers` (`workflows/snapshot.py`).
* Alle relevanten Schritte (1,2,3,4,6,7) benutzen diese Funktion, um doppelte Logik zu vermeiden.

Vorteile:
* Vorhersehbare Priorität ohne nachträgliche „Hack“-Korrekturen.
* Minimierte Redundanz in per-molekül TOMLs (nur spezifisches definieren, globales bleibt im Projektfile).
* Klare Trennung zwischen „Ändern“ (Override) und „Auffüllen“ (Fill-In).

Best Practices:
* Molekül-spezifische Parameter (Ladung, Forcefield, besondere Sampling-Optionen) in der Molekül-Datei definieren.
* Einheitliche globale Pfade / Laufzeit-Defaults nur im Projektroot pflegen.
* `--config` für experimentelle Runs oder CI-spezifische Tuning verwenden (nicht dauerhaft einchecken).

Edge Cases:
* Fehlt ein Molekül-TOML komplett, greifen Prozess/CLI/Projekt in genannter Reihenfolge.
* Mehrere Moleküle: bestimmte restriktive Keys werden geloggt, wenn sie in Multi-Mol Kontext überschrieben würden (siehe `RESTRICTED_DEFAULT`).

TBD (add your license information here).
