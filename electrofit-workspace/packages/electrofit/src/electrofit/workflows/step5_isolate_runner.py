"""Isolated single conformer runner.

Invoked as a subprocess to sandbox potential hard exits (sys.exit / native abort)
that would otherwise terminate the parent step5 driver. It reuses process_one
and prints a single RESULT line parsable by the parent.
"""
from __future__ import annotations
import argparse, json, sys
from electrofit.workflows.step5_worker import process_one


def main(argv: list[str] | None = None):  # pragma: no cover - diagnostic utility
    ap = argparse.ArgumentParser()
    ap.add_argument("--conf", required=True)
    ap.add_argument("--project", required=True)
    ap.add_argument("--override")
    ap.add_argument("--multi-mol", action="store_true")
    ap.add_argument("--mock", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args(argv)
    rel, ok, msg = process_one(
        args.conf,
        args.project,
        args.override,
        args.multi_mol,
        args.mock,
        args.verbose,
    )
    rec = {"rel": rel, "ok": ok, "msg": msg}
    print("RESULT:" + json.dumps(rec), flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":  # pragma: no cover
    main()
