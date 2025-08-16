from pathlib import Path
import os
import logging
import traceback
from electrofit.config.loader import load_config, dump_config
from electrofit.workflows.snapshot import build_snapshot_with_layers
from electrofit.logging import log_run_header, reset_logging, setup_logging
from electrofit.core.process_conform import process_conform
from electrofit.io.files import find_file_with_extension, strip_extension


def process_one(conf_dir_str: str,
                project_root_str: str,
                override_cfg_path: str | None,
                multi_mol: bool,
                mock: bool,
                verbose: bool) -> tuple[str, bool, str]:
    conf_dir = Path(conf_dir_str)
    project_root = Path(project_root_str)
    prev = os.getcwd()
    try:
        os.chdir(conf_dir)
        # Advanced one-time diagnostics
        try:
            import sys as _sys_d, os as _os_d, faulthandler as _fh_d, signal as _sig_d, atexit as _ax_d, gc as _gc_d
            if not getattr(_sys_d.modules.get(__name__), '_EF_ADV_DIAG', False):
                try:
                    _fh_d.enable(file=_sys_d.stderr)
                except Exception:
                    pass
                def _ax():
                    try:
                        print('[debug-atexit] worker normal shutdown', flush=True)
                    except Exception:
                        pass
                try:
                    _ax_d.register(_ax)
                except Exception:
                    pass
                def _sig_handler(sig, frame):
                    try:
                        print(f'[debug-signal] sig={sig}', flush=True)
                        _fh_d.dump_traceback(file=_sys_d.stderr)
                    except Exception:
                        pass
                for _s in [getattr(_sig_d, n, None) for n in ('SIGTERM','SIGINT','SIGSEGV')]:
                    if _s is not None:
                        try:
                            _sig_d.signal(_s, _sig_handler)
                        except Exception:
                            pass
                try:
                    if hasattr(_gc_d, 'callbacks'):
                        _gc_d.callbacks.append(lambda phase, info: phase=='stop' and print('[debug-gc] cycle', flush=True))  # type: ignore
                except Exception:
                    pass
                setattr(_sys_d.modules[__name__], '_EF_ADV_DIAG', True)
        except Exception:
            pass
        # Monkeypatch sys.exit / os._exit for diagnostics
        import sys as _sys, os as _os
        _orig_sys_exit = _sys.exit
        _orig_os_exit = _os._exit
        def _dbg_sys_exit(code=0):
            try:
                print(f"[debug-sys-exit] code={code}", flush=True)
            except Exception:
                pass
            return _orig_sys_exit(code)
        def _dbg_os_exit(code=0):
            try:
                print(f"[debug-os-exit] code={code}", flush=True)
            except Exception:
                pass
            return _orig_os_exit(code)
        _sys.exit = _dbg_sys_exit  # type: ignore
        _os._exit = _dbg_os_exit   # type: ignore
        # logging per conformer
        reset_logging()
        log_path = Path(os.getcwd()) / "process.log"
        setup_logging(str(log_path), also_console=verbose)
        log_run_header("step5")
        pdb_file = find_file_with_extension("pdb")
        if not pdb_file:
            raise FileNotFoundError("No PDB file found in conformer directory")
        mol_name = conf_dir.parent.parent.name if len(conf_dir.parts) >= 2 else conf_dir.parent.name
        build_snapshot_with_layers(
            conf_dir,
            project_root,
            mol_name,
            multi_molecule=multi_mol,
            log_fn=logging.info,
            upstream=conf_dir.parent / "electrofit.toml",  # parent is extracted_conforms root
            process_cfg=conf_dir.parent.parent / "electrofit.toml" if len(conf_dir.parts) >= 3 else None,
            molecule_input=project_root / "data" / "input" / mol_name / "electrofit.toml",
            project_defaults=project_root / "electrofit.toml",
            extra_override=Path(override_cfg_path) if override_cfg_path else None,
        )
        cfg = load_config(project_root, context_dir=conf_dir)
        dump_config(cfg, header=False, log_fn=logging.info)
        proj = cfg.project
        molecule_name = proj.molecule_name or strip_extension(pdb_file)
        if mock:
            with open("executed.txt", "w") as f:
                f.write(f"run{molecule_name}")
        else:
            print(f"[debug-before-process_conform] pdb={pdb_file}", flush=True)
            process_conform(
                molecule_name=molecule_name,
                pdb_file=pdb_file,
                base_scratch_dir=cfg.paths.base_scratch_dir or "/tmp/electrofit_scratch",
                net_charge=proj.charge or 0,
                residue_name=proj.residue_name or "LIG",
                adjust_sym=getattr(proj, "adjust_symmetry", False),
                protocol=proj.protocol or "bcc",
                ignore_sym=getattr(proj, "ignore_symmetry", False),
            )
            print(f"[debug-after-process_conform] pdb={pdb_file}", flush=True)
        rel = conf_dir.relative_to(project_root)
        print(f"[worker-return] {rel} ok", flush=True)
        return (str(rel), True, "ok")
    except Exception as e:  # pragma: no cover
        try:
            tb = traceback.format_exc()
            logging.error(f"Worker failure {conf_dir}: {e}\n{tb}")
        except Exception:
            pass
        try:
            rel = conf_dir.relative_to(project_root)
            rel_str = str(rel)
        except Exception:
            rel_str = conf_dir_str
        print(f"[worker-exc] {rel_str}: {e}", flush=True)
        return (rel_str, False, str(e))
    finally:
        # Restore patched exits
        try:
            import sys as _sys, os as _os
            if '_orig_sys_exit' in locals():
                _sys.exit = _orig_sys_exit  # type: ignore
            if '_orig_os_exit' in locals():
                _os._exit = _orig_os_exit   # type: ignore
        except Exception:
            pass
        try:
            os.chdir(prev)
        except Exception:
            pass
