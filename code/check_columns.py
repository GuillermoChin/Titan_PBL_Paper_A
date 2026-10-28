#!/usr/bin/env python
"""Verify the single-column assumption of the archived 1D runs.

The archived wrfout files carry `west_east = 2`, because WRF always allocates a
minimum horizontal stencil even for a single-column integration. Every analysis
in Paper A reads index `[..., 0, 0]`, which is only legitimate if the two
columns are identical, or if the second is a halo that carries no independent
solution.

This script checks that assumption explicitly rather than assuming it, and
exits non-zero if the columns diverge.

    python code/check_columns.py [path_to_wrfout]

The archive is located with TITAN_DATA_DIR (default: <package>/data); the
output files may carry the published name or the NTFS-safe one (':' -> '_').
"""
from __future__ import annotations
import sys
from pathlib import Path
import os
import numpy as np
import xarray as xr

DATA_DIR = os.environ.get("TITAN_DATA_DIR", str(Path(__file__).resolve().parents[1] / "data"))

RUNS = [
    "Dry Simulations/1D_Huygens_TMN91_5sol",
    "Dry Simulations/1D_Huygens_baseline_5sol",
    "Dry Simulations/1D_Huygens_0.1flhc_5sol",
]

# Variables que el analisis del Paper A lee efectivamente.
CHECKED = ["TSK", "T_PHY", "PBLH", "HFX", "FLHC", "FLUX_LW_D", "FLUX_LW_U", "P", "PB", "PH", "PHB", "PSFC"]


def resolve(relative: str) -> Path | None:
    directory = Path(DATA_DIR) / relative
    if not directory.is_dir():
        return None
    matches = sorted(directory.glob("wrfout_d01*"))
    return matches[0] if matches else None


def check(path: Path) -> tuple[int, int]:
    """Return (variables compared, variables whose two columns differ)."""
    ds = xr.open_dataset(path, decode_times=None)
    compared = differing = 0
    for name in CHECKED:
        if name not in ds.variables:
            continue
        array = ds[name]
        if "west_east" not in array.dims or ds.sizes.get("west_east", 1) < 2:
            continue
        values = np.asarray(array.values, dtype=np.float64)
        axis = array.dims.index("west_east")
        first = np.take(values, 0, axis=axis)
        second = np.take(values, 1, axis=axis)
        compared += 1
        if not np.allclose(first, second, rtol=0, atol=0):
            differing += 1
            worst = float(np.max(np.abs(first - second)))
            print(f"    {name:<12} columns DIFFER, max |delta| = {worst:.6g}")
        else:
            print(f"    {name:<12} columns identical")
    ds.close()
    return compared, differing


def main() -> int:
    targets = [Path(sys.argv[1])] if len(sys.argv) > 1 else [
        p for p in (resolve(r) for r in RUNS) if p is not None
    ]
    if not targets:
        print(f"No archive run found under {DATA_DIR}.")
        print("Set TITAN_DATA_DIR or pass the path of a wrfout file.")
        return 2

    total_differing = 0
    for path in targets:
        print(f"\n{path.parent.name}")
        compared, differing = check(path)
        total_differing += differing
        print(f"    -> {compared} variables compared, {differing} differing")

    print()
    if total_differing:
        print(f"FAIL: {total_differing} variables with differing columns. "
              "Reading [..., 0, 0] is not equivalent to reading the physical column.")
        return 1
    print("OK: the two columns are identical in every variable checked; "
          "reading index 0 is reading the single physical column.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
