# Dependency Inventory

Erzeugt: 2025-08-17T12:56:33.073001+00:00  
Quelle JSON: `tools/inventory/dependency_graph.json`

## Pakete

### electrofit (electrofit)

Version: 0.1.0

Dependencies (deklariert):

- app
- matplotlib
- mdtraj
- numpy
- openbabel-wheel
- openmol
- pandas
- typer
- tomli>=2.0.1; python_version < "3.11"
- acpype
- tqdm
- tomli_w
Scripts:

- electrofit -> electrofit.cli.app:main

### electrofit-fep (electrofit-fep)

Version: 0.1.0

Dependencies (deklariert):

- electrofit>=0.1.0
- typer>=0.12
- rich>=13
Scripts:

- electrofit-fep -> electrofit_fep.cli:main

### electrofit-analysis (electrofit-analysis)

Version: 0.1.0

Dependencies (deklariert):

- MDAnalysis
- alchemlyb
- joblib
- matplotlib
- numpy
- pandas
- rdkit
- scikit-learn
- scipy
- seaborn
- typer
Scripts:

- electrofit-analysis -> electrofit_analysis.cli:main

## Paket-Kanten (interne Cross-Package Imports)

- electrofit_analysis -> electrofit: 13
- electrofit_fep -> electrofit: 12

## Externe Top-Level Imports (Häufigkeit)

- os: 53
- logging: 34
- matplotlib: 32
- argparse: 25
- numpy: 24
- pathlib: 22
- sys: 17
- __future__: 13
- subprocess: 12
- re: 12
- seaborn: 12
- shutil: 10
- typing: 9
- json: 9
- pandas: 9
- rdkit: 9
- alchemlyb: 7
- glob: 5
- scipy: 5
- math: 5
- contextlib: 4
- fnmatch: 4
- concurrent: 4
- pmx: 4
- traceback: 3
- tomli: 3
- tomllib: 3
- hashlib: 3
- warnings: 3
- atexit: 2
- signal: 2
- copy: 2
- time: 2
- uuid: 2
- filecmp: 2
- multiprocessing: 2
- typer: 2
- itertools: 2
- MDAnalysis: 2
- io: 1
- importlib: 1
- runpy: 1
- threading: 1
- dataclasses: 1
- tomli_w: 1
- openmol: 1
- openbabel: 1
- tempfile: 1
- viz: 1
- mdtraj: 1
- random: 1
- tqdm: 1
- faulthandler: 1
- gc: 1
- joblib: 1
- sklearn: 1
- types: 1

## Modulzählung

Anzahl Module: 83

## Hinweise

- Dynamische / optionale Imports sind nicht speziell markiert.
- Nutzung (Symbole) nicht aufgelöst; es handelt sich nur um Modulreferenzen.
- Bei Refactor erneut Skript ausführen und Änderungen committen.