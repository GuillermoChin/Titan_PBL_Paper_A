#!/usr/bin/env python
"""
r2_common.py -- shared helpers of the archive and TAM analyses (stage r2).

Paths of the package, the import guard (titan_pbl must come from
code/titan_pbl_subset), the TAM file list with the MD5 published by Zenodo
(the analyses stop if a file differs), and loaders for the mtWRF archive runs.

External data (not in the package):
  TITAN_DATA_DIR  directory that contains "Dry Simulations/" of the archive of
                  Hartwick et al. (2026), doi:10.5281/zenodo.17246644
                  (default: <package>/data).
  TITAN_TAM_DIR   directory with the eight atmos_diurnal_*.nc files of Han and
                  Lora (2025), doi:10.5281/zenodo.16791001
                  (default: $TITAN_DATA_DIR/TAM Simulations_zenodo_16791001).
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT / "code" / "titan_pbl_subset") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "code" / "titan_pbl_subset"))

import titan_pbl  # noqa: E402

if (REPO_ROOT / "code" / "titan_pbl_subset") not in Path(titan_pbl.__file__).resolve().parents:
    raise ImportError(f"titan_pbl imported from {titan_pbl.__file__}, not from code/titan_pbl_subset.")

import numpy as np  # noqa: E402

AUDIT = REPO_ROOT / "code" / "analysis"
INPUTS = AUDIT / "inputs"
BUILD = REPO_ROOT / "build"
R1_OUT = BUILD / "r1"
R2_OUT = BUILD / "r2"
QUICKLOOK = R2_OUT / "quicklook"
DATA_R1 = BUILD / "source_a"
DATA_R2 = BUILD / "r2"
NUMBERS_R1 = DATA_R1 / "numbers_r1.json"
NUMBERS_R2 = DATA_R2 / "numbers_r2.json"

DATA_ROOT = Path(os.environ.get("TITAN_DATA_DIR", str(REPO_ROOT / "data")))
TAM_DIR = Path(os.environ.get("TITAN_TAM_DIR", str(DATA_ROOT / "TAM Simulations_zenodo_16791001")))
TAM_CASES = [f"case{i}_{site}" for site in ("huygens", "dragonfly") for i in (1, 2, 3, 4)]
TAM_DAY_S = 86400.0  # "days since 0000-00-00": dia de FMS = 86400 s

ARCHIVE_DIR = DATA_ROOT / "Dry Simulations"


def wrfout(run: str) -> Path:
    """wrfout of an archive run, under its published name or the NTFS-safe one (':' -> '_')."""
    d = ARCHIVE_DIR / run
    published = d / "wrfout_d01_2007-00001_09:30:00"
    return published if published.exists() else d / "wrfout_d01_2007-00001_09_30_00"


ARCHIVE_RUNS = {
    "TMN91": wrfout("1D_Huygens_TMN91_5sol"),
    "baseline5sol": wrfout("1D_Huygens_baseline_5sol"),
    "0.1flhc": wrfout("1D_Huygens_0.1flhc_5sol"),
}

SIGMAS = {
    "CODATA_5.670374419e-8": 5.670374419e-8,
    "5.67e-8": 5.67e-8,
    "FMS_5.6734e-8": 5.6734e-8,
}


def guard_info() -> dict:
    return {"titan_pbl_file": str(Path(titan_pbl.__file__).resolve()), "repo_root": str(REPO_ROOT)}


def rel(path: Path) -> str:
    """Path relative to the package; external data as $TITAN_DATA_DIR/... or $TITAN_TAM_DIR/... (no absolute paths)."""
    p = Path(path).resolve()
    for base, label in ((TAM_DIR, "$TITAN_TAM_DIR"), (DATA_ROOT, "$TITAN_DATA_DIR"), (REPO_ROOT, None)):
        try:
            r = str(p.relative_to(Path(base).resolve())).replace("\\", "/")
        except ValueError:
            continue
        return r if label is None else f"{label}/{r}"
    return p.name


def tam_path(case: str) -> Path:
    return TAM_DIR / f"atmos_diurnal_{case}.nc"


def tam_registered_md5() -> dict[str, str]:
    """MD5 published by Zenodo for record 16791001 (inputs/tam_16791001_md5.txt)."""
    text = (INPUTS / "tam_16791001_md5.txt").read_text(encoding="utf-8")
    out = {}
    for m in re.finditer(r"^(atmos_diurnal_\S+\.nc) md5=([0-9a-f]{32})", text, re.M):
        out[m.group(1)] = m.group(2)
    return out


def check_tam_md5() -> dict:
    """Stop if any TAM file differs from the MD5 published by Zenodo."""
    registered = tam_registered_md5()
    result = {}
    for case in TAM_CASES:
        path = tam_path(case)
        md5 = hashlib.md5(path.read_bytes()).hexdigest()
        ok = registered.get(path.name) == md5
        result[path.name] = {"md5": md5, "registered_R1": registered.get(path.name), "match": ok}
        if not ok:
            raise SystemExit(f"MD5 of {path.name} differs from the one published by Zenodo ({md5} vs {registered.get(path.name)})")
    return result


def load_figure_module():
    """Import code/analysis/source_a_core.py unchanged, writing its outputs to a temp dir."""
    os.environ["PAPER_A_OUT"] = tempfile.mkdtemp(prefix="r2_fig_")
    spec = importlib.util.spec_from_file_location(
        "make_paper_A_figures_r2", REPO_ROOT / "code" / "analysis" / "source_a_core.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_r1_flux_module():
    """Import code/analysis/r1_flux_comparison.py unchanged (interface reconstruction)."""
    spec = importlib.util.spec_from_file_location("r1_flux_comparison", AUDIT / "r1_flux_comparison.py")
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def numbers_r1() -> dict:
    return json.loads(NUMBERS_R1.read_text(encoding="utf-8"))


def number_r1(key: str):
    for n in numbers_r1()["numbers"]:
        if n["key"] == key:
            return n["value"]
    raise KeyError(key)


def update_numbers_r2(entries: list[dict], prefix: str | tuple[str, ...], generated_by: str) -> None:
    """Add/replace entries in numbers_r2.json (same schema as numbers_r1.json).

    Entries whose key starts with ``prefix`` are replaced. Keys already present
    in numbers_r1.json are rejected.
    """
    DATA_R2.mkdir(parents=True, exist_ok=True)
    r1_keys = {n["key"] for n in numbers_r1()["numbers"]}
    clash = [e["key"] for e in entries if e["key"] in r1_keys]
    if clash:
        raise SystemExit(f"Claves de numbers_r2 que repiten claves de R1: {clash}")
    if NUMBERS_R2.is_file():
        data = json.loads(NUMBERS_R2.read_text(encoding="utf-8"))
    else:
        data = {"generated_by": [], "sol_length_s": numbers_r1()["sol_length_s"], "numbers": []}
    data["numbers"] = [n for n in data["numbers"] if not n["key"].startswith(prefix)] + entries
    if generated_by not in data["generated_by"]:
        data["generated_by"].append(generated_by)
    keys = [n["key"] for n in data["numbers"]]
    if len(keys) != len(set(keys)):
        raise SystemExit("Claves duplicadas en numbers_r2.json")
    NUMBERS_R2.write_text(json.dumps(data, indent=2, ensure_ascii=False, default=_json_default), encoding="utf-8")


def entry(key, value, units, boundary, variant, source_file, source_function, sig, why) -> dict:
    if isinstance(value, (np.floating, np.integer)):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        value = None
    return {
        "key": key,
        "value": value,
        "units": units,
        "lw_boundary": boundary,
        "sw_variant": variant,
        "source_file": source_file,
        "source_function": source_function,
        "print_sig_figs": sig,
        "print_sig_figs_reason": why,
    }


def _json_default(obj):
    if isinstance(obj, (np.floating,)):
        v = float(obj)
        return v if np.isfinite(v) else None
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, Path):
        return str(obj)
    raise TypeError(type(obj))


def dump_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean_nan(obj), indent=2, ensure_ascii=False, default=_json_default),
                    encoding="utf-8")


def clean_nan(obj):
    """Replace non-finite floats by None recursively (strict JSON)."""
    if isinstance(obj, dict):
        return {k: clean_nan(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean_nan(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return clean_nan(obj.tolist())
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    if isinstance(obj, np.integer):
        return int(obj)
    return obj
