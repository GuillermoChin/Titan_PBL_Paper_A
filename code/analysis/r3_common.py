#!/usr/bin/env python
"""
r3_common.py -- shared helpers of the stage-r3 analyses and of the single source.

Paths of the package, the import guard (re-applied through r2_common), the
mtWRF constants inferred from the archive, and small I/O helpers. Everything
is written under build/, except the single source (numbers/) and the figures
(figures/).
"""

from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT / "code" / "titan_pbl_subset") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "code" / "titan_pbl_subset"))

import titan_pbl  # noqa: E402

# Guard de R1 (P3 de R3): aborta si titan_pbl no viene del repo.
if (REPO_ROOT / "code" / "titan_pbl_subset") not in Path(titan_pbl.__file__).resolve().parents:
    raise ImportError(f"titan_pbl imported from {titan_pbl.__file__}, not from code/titan_pbl_subset.")

AUDIT = REPO_ROOT / "code" / "analysis"
if str(AUDIT) not in sys.path:
    sys.path.insert(0, str(AUDIT))

import numpy as np  # noqa: E402

import r2_common as C2  # noqa: E402

R1_OUT = C2.R1_OUT
R2_OUT = C2.R2_OUT
R3_OUT = REPO_ROOT / "build" / "r3"
QUICKLOOK = R3_OUT / "quicklook"
DATA_R1 = C2.DATA_R1
DATA_R2 = C2.DATA_R2
DATA_REV = REPO_ROOT / "build" / "numbers"
FIG_DIR = REPO_ROOT / "figures"
INPUTS = C2.INPUTS
NUMBERS_R1 = C2.NUMBERS_R1
NUMBERS_R2 = C2.NUMBERS_R2
NUMBERS_R3 = R3_OUT / "numbers_r3_staging.json"
NUMBERS_DIR = REPO_ROOT / "numbers"
NUMBERS_REV = NUMBERS_DIR / "numbers_paper_A.json"
ARCHIVE_RUNS = C2.ARCHIVE_RUNS

# Indices de registros por Tsol en los wrfout (los mismos de R2/T3b).
SOL_INDICES = {2: (87, 231), 3: (231, 375), 4: (375, 519), 5: (519, 663)}
CORE = slice(86, 664)  # incluye un registro a cada lado para diferencias centradas

# Constantes de R1/R2 (no se cambian para comparabilidad).
GRAVITY = 1.352
R_GAS_R1 = 296.8
CP_HARTWICK = 1101.0
SOL_S = 1378080.0

dump_json = C2.dump_json
clean_nan = C2.clean_nan
rel = C2.rel
guard_info = C2.guard_info


def git_head() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO_ROOT,
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:  # pragma: no cover - solo informativo
        return "unknown"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = []
        for r in rows:
            for k in r:
                if k not in fieldnames:
                    fieldnames.append(k)
    with path.open("w", newline="", encoding="utf-8") as fh:
        # Fin de linea LF: el paquete exige salidas identicas byte a byte en cualquier sistema.
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: _fmt(v) for k, v in r.items()})


def _fmt(v):
    # float(...) evita que los np.float64 se escriban como "np.float64(...)".
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (float, np.floating)):
        return repr(float(v))
    if isinstance(v, np.integer):
        return int(v)
    return v


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def numbers(path: Path) -> list[dict]:
    return json.loads(path.read_text(encoding="utf-8"))["numbers"]


def stage_numbers(entries: list[dict], prefix: str, generated_by: str) -> None:
    """Add/replace R3 entries (keys starting with ``prefix``) in the R3 staging file.

    The staging file is consolidated with numbers_r1/r2 into the single source
    in T6. Each entry: key, value, uncertainty, units, lw_boundary, sw_variant,
    source_file, source_function, note.
    """
    if NUMBERS_R3.is_file():
        data = json.loads(NUMBERS_R3.read_text(encoding="utf-8"))
    else:
        data = {"generated_by": [], "numbers": []}
    data["numbers"] = [n for n in data["numbers"] if not n["key"].startswith(prefix)] + entries
    if generated_by not in data["generated_by"]:
        data["generated_by"].append(generated_by)
    keys = [n["key"] for n in data["numbers"]]
    if len(keys) != len(set(keys)):
        raise SystemExit("Claves duplicadas en numbers_r3_staging.json")
    dump_json(NUMBERS_R3, data)


def r3_entry(key, value, units, source_file, source_function, *, uncertainty=None,
             boundary="n/a", variant="n/a", note="") -> dict:
    return clean_nan({
        "key": key, "value": value, "uncertainty": uncertainty, "units": units,
        "lw_boundary": boundary, "sw_variant": variant,
        "source_file": source_file, "source_function": source_function, "note": note,
    })
