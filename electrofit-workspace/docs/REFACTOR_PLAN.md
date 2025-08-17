# Elektrofit Refactor & Architektur Leitfaden

Stand: 2025-08-17

Dieses Dokument dient als lebende Referenz für die laufende Struktur‑Konsolidierung und Architektur‑Verbesserungen. Ziel: Klare Schichten, minimale Kopplung, eindeutige Verantwortlichkeiten, konsistente Namensgebung und wartbare Public API.

---
## Leitprinzipien
1. Single Source of Truth: Konfiguration & Layering nur an einer Stelle (compose/build Snapshot).
2. Dünne Workflows: Step-* Dateien orchestrieren nur, keine Geschäftslogik.
3. Domain vor Technik: Fachlogik (Symmetrie, Ladungen, Konformer) in `domain/*` statt in generischen Sammelordnern.
4. Infrastruktur kapseln: Logging, Scratch, Config, externe Tool-Adapter getrennt in `infra/` bzw. `adapters/`.
5. Öffentliche vs. Interne API: Explizit via `__all__` + Namenskonvention (`_internal`).
6. Predictable Imports: Pipelines importieren Services, nicht Einzel-Funktionsmodule querbeet.
7. Testbarkeit: Domain-Funktionen ohne Seiteneffekte (kein global logging, kein cwd change) testbar.
8. Idempotenz & Reentrance: Wiederholter Aufruf eines Steps verändert nur notwendige Artefakte.
9. Minimaler Log-Lärm: Keine redundanten Merge-Logs pro Konformer.
10. Optionalität: Plotting & Visualisierung als optionale Extras (lazy import / optional dependency).

---
## Ziel-Schichten (Soll-Struktur)
```
electrofit/
  api/                # stabile Fassaden (optional)
  cli/                # Argument-Parsing + Dispatcher
  pipeline/steps/     # step1.py ... stepN.py (nur Orchestrierung)
  domain/
    symmetry/
    charges/
    conformers/
    prep/
  infra/
    config_snapshot.py
    logging.py
    scratch_manager.py
  adapters/
    gromacs.py
    gaussian.py (potenziell)
    acpype.py (potenziell)
  io/
    mol2.py
    resp.py
    files.py
  viz/ (oder plotting/)
  internal/ (falls nötig für Migrationsphase)
```

---
## Aktueller Ist-Zustand (Kurz)
- `workflows/` enthält Schritte + Snapshot-Builder (gemischt: Orchestrierung & Infrastruktur).
- `core/` mischt charges/symmetry/prep.
- `external/` nur GROMACS Adapter (Name zu generisch).
- `logging.py` top-level; `scratch/manager.py` eigenes Mini-Paket.
- `utils_curly_brace.py` unklare Funktion / Name.

---
## Phasenplan
### Phase 1 – Reorganisation (Low Risk)
- [x] Verschiebe `workflows/snapshot.py` → `infra/config_snapshot.py` (Export: `compose_snapshot` (alias für `build_snapshot_with_layers`), `CONFIG_ARG_HELP`).
- [x] Umbenenne / entferne `external/` → `adapters/` (alter Namespace gelöscht, kein Shim nötig).
- [x] Verschiebe `logging.py` → `infra/logging.py`.
- [x] Verschiebe `scratch/manager.py` → `infra/scratch_manager.py`; entferne Paket `scratch`.
- [ ] Benenne `utils_curly_brace.py` → `templating.py` (oder integriere in `io/`).
- [ ] Re-Exports für alte Pfade mit DeprecationWarning (nur falls extern genutzt; aktuell nicht benötigt für GROMACS Adapter).

### Phase 2 – Domain-Isolation
- [ ] Erzeuge `domain/symmetry/` (move `symmetry.py`, `equiv_groups.py`).
- [ ] Erzeuge `domain/charges/` (move `process_conform.py` + spätere Aggregationslogik).
- [ ] Erzeuge `domain/prep/` (move `process_initial_structure.py`).
- [ ] Ergänze Services: `services/conformer_sampling.py`, `services/config_service.py` (Wrapper um Snapshot + Loader).

### Phase 3 – API & Services
- [ ] `pipeline/steps/stepX.py` einführen; vorhandene `stepX_*.py` umbenennen / verschieben.
- [ ] Öffentlicher Einstieg: `electrofit.api` definiert stabile Funktionen (z.B. `run_step(step:int, project:Path, config:Path|None)`).
- [ ] Interne Module prefix `_`; unveröffentlichte Objekte entfernen aus `__all__`.

### Phase 4 – Bereinigung & Deprecation
- [ ] Entferne Legacy-Dateien / Aliase (`run_process_*`).
- [ ] Entferne `merge_into_snapshot` aus Public Surface (nur intern verfügbar).
- [ ] Einheitliche Schritt-Namen: `step1.py`, `step2.py` ...
- [ ] Dokumentiere neue Modulkarte in README / docs/architecture.md.

### Phase 5 – Optional Enhancements
- [ ] Optionales Extra: `pip install electrofit[viz]` für plotting.
- [ ] Typsicherheit: mypy-Konfiguration + strictere Typen in infra & domain.
- [ ] Structured Logging (JSON) optional für CI.
- [ ] Performance Messpunkte (Timer-Dekoratoren) in Services.

---
## Migrationsstrategie
1. Verschieben mit Re-Exports + DeprecationWarnings.
2. CI Test-Suite sichern (jede Phase grün bevor nächste startet).
3. Release Notizen: „0.2.x – interne Reorganisation, API unverändert / 0.3.0 – neue API-Fassade“.
4. Nach zwei Minor-Releases: Entfernen der Deprecation-Re-Exports.

---
## Naming Guidelines
- Dateien: snake_case, prägnant (keine überflüssigen Präfixe wie `step3_start_sim` → `step3.py`).
- Funktionen: Verben für Aktionen (`compose_snapshot`, `sample_conformers`), Substantive für reine Modelle.
- Keine Doppelungen in Pfad + Modul (kein `symmetry/symmetry.py` → besser `symmetry/core.py` oder `symmetry/ops.py`).
- Interne Implementierungen: `_helper.py` oder `_function` Präfix.

---
## Logging Guidelines
- Workflows: Nur High-Level Fortschritt.
- Domain/Services: Debug-fähige Detail Logs via Logger `electrofit.<sub>`.
- Konfig-Layering Log-Kategorien: `[config][override]`, `[config][fill]`, `[config][warn]`.

---
## Offene Fragen / ToDo klären
- Benötigen wir einen stabilen Python API-Modus außerhalb CLI? (Batch/Library Use)
- Sollen Konformer-Snapshots langfristig wegfallen (nur Root-Snapshot)?
- Einführen eines Settings-Objekts statt direkte weitergabe von Path/Sting?

---
## Aktueller Status (Fortschritt markieren)
| Phase | Aufgabe | Status |
|-------|---------|--------|
| 1 | snapshot -> infra | done |
| 1 | external -> adapters (Namespace entfernt) | done |
| 1 | logging -> infra | done |
| 1 | scratch -> infra | done |
| 1 | utils_curly_brace rename | open |
| 2 | domain symmetry | open |
| 2 | domain charges | open |
| 2 | domain prep | open |
| 3 | pipeline/steps Struktur | open |
| 3 | api facade | open |
| 4 | remove legacy aliases | open |
| 4 | hide merge_into_snapshot | partial (intern markiert) |

(Beim Fortschritt jeweils "open | in-progress | done" einsetzen.)

---
## Pflege
Dieses Dokument wird bei jeder abgeschlossenen Umsetzungs-Phase aktualisiert. Bitte Pull Requests mit Änderungen an Architektur / Struktur verlinken und Status-Tabelle pflegen.

