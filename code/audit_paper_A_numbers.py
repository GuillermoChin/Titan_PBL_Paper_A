#!/usr/bin/env python
"""
audit_paper_A_numbers.py -- verify the single source of numbers and the released files.

1. Regenerates the single source of numbers from the analysis outputs under
   build/ (run code/run_all.py first, or pass --full to run it here) into
   build/audit/, without touching numbers/.
2. Compares the regenerated file with numbers/numbers_paper_A.json: byte for
   byte and key by key (full value, printed value, status, uncertainty).
3. Checks the SHA-256 of every file under numbers/, figures/ and tables/
   against the list in PROVENANCE.md (the files of the release).

Exits with code 0 only if nothing changed; any difference is listed and the
exit code is 1 (2 if the analysis outputs are missing).

Usage:
    python code/audit_paper_A_numbers.py [--full]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "numbers" / "numbers_paper_A.json"
AUDIT_DIR = ROOT / "build" / "audit"
PROVENANCE = ROOT / "PROVENANCE.md"
CHECKED_DIRS = ("numbers", "figures", "tables")


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def load_numbers(path: Path) -> dict:
    """The single source as a dict; ValueError if the file is not valid UTF-8 JSON (a BOM included)."""
    try:
        return json.loads(path.read_bytes().decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{path.relative_to(ROOT).as_posix()} is not valid UTF-8 JSON ({exc})") from exc


def compare_numbers(ref: Path, new: Path) -> list[str]:
    try:
        a, b = load_numbers(ref), load_numbers(new)
    except ValueError as exc:
        return [str(exc)]
    ea, eb = {e["key"]: e for e in a["numbers"]}, {e["key"]: e for e in b["numbers"]}
    out = [f"key only in the released file: {k}" for k in sorted(set(ea) - set(eb))]
    out += [f"key only in the regenerated file: {k}" for k in sorted(set(eb) - set(ea))]
    for k in sorted(set(ea) & set(eb)):
        for f in ("value", "print", "status", "uncertainty"):
            if json.dumps(ea[k].get(f)) != json.dumps(eb[k].get(f)):
                out.append(f"{k}.{f}: released {ea[k].get(f)!r} -> regenerated {eb[k].get(f)!r}")
    if not out and sha(ref) != sha(new):
        out.append("same keys and values, but the files differ byte for byte (metadata or text fields)")
    return out


def check_provenance() -> tuple[list[str], int]:
    if not PROVENANCE.is_file():
        return ["PROVENANCE.md not found: release hashes not checked"], 0
    listed = {}
    for m in re.finditer(r"^\|\s*`([^`]+)`\s*\|\s*([0-9a-f]{64})\s*\|", PROVENANCE.read_text(encoding="utf-8"), re.M):
        listed[m.group(1)] = m.group(2)
    out, n = [], 0
    for rel, h in sorted(listed.items()):
        if not rel.startswith(CHECKED_DIRS):
            continue
        n += 1
        p = ROOT / rel
        if not p.is_file():
            out.append(f"missing: {rel}")
        elif sha(p) != h:
            out.append(f"changed: {rel}")
    for d in CHECKED_DIRS:
        for p in sorted((ROOT / d).rglob("*")):
            rel = p.relative_to(ROOT).as_posix()
            if p.is_file() and rel not in listed:
                out.append(f"not in PROVENANCE.md: {rel}")
    return out, n


def main() -> int:
    ap = argparse.ArgumentParser(description="Verify the single source of numbers and the released files.")
    ap.add_argument("--full", action="store_true", help="run code/run_all.py first (needs the external data)")
    args = ap.parse_args()
    if args.full:
        rc = subprocess.run([sys.executable, str(ROOT / "code" / "run_all.py")], cwd=ROOT, timeout=3600).returncode
        if rc != 0:
            print("run_all.py failed")
            return 2
    needed = [ROOT / "build" / "source_a" / "numbers_r1.json", ROOT / "build" / "r2" / "numbers_r2.json",
              ROOT / "build" / "r3" / "numbers_r3_staging.json"]
    if not all(p.is_file() for p in needed):
        print("The analysis outputs are missing (build/). Run `python code/run_all.py` first, or use --full.")
        return 2
    AUDIT_DIR.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PAPER_A_NUMBERS_OUT_DIR=str(AUDIT_DIR), PYTHONDONTWRITEBYTECODE="1")
    rc = subprocess.run([sys.executable, str(ROOT / "code" / "build_numbers_paper_A.py")], cwd=ROOT, env=env,
                        capture_output=True, text=True, timeout=600).returncode
    if rc != 0:
        print("build_numbers_paper_A.py failed")
        return 2
    diffs = compare_numbers(REFERENCE, AUDIT_DIR / "numbers_paper_A.json")
    n_keys = len(load_numbers(AUDIT_DIR / "numbers_paper_A.json")["numbers"])
    prov, n_files = check_provenance()
    print(f"single source: {n_keys} keys regenerated and compared -> "
          f"{'identical' if not diffs else f'{len(diffs)} differences'}")
    for d in diffs:
        print("  ", d)
    print(f"release files: {n_files} files of numbers/, figures/ and tables/ checked against PROVENANCE.md -> "
          f"{'identical' if not prov else f'{len(prov)} issues'}")
    for d in prov:
        print("  ", d)
    return 0 if not diffs and not prov else 1


if __name__ == "__main__":
    raise SystemExit(main())
