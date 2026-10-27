#!/usr/bin/env python
"""
run_all.py -- reproduce the single source of numbers, the figures and the tables.

Runs every analysis in order, from scratch (build/ is deleted first unless
--keep-build is given):

  stage source_a  framework diagnostics of Sections 2-3 (code/analysis/source_a_core.py,
                  r1_regenerate_source_a.py, r1_threshold_verification.py) and the
                  flux comparison on the archive (r1_flux_comparison.py)
  stage r2        archive and TAM analyses (r2_t*.py)
  stage r3        lower boundary of the MYJ scheme and column closure, flux comparison
                  in three runs, robustness, radiative timescales (r3_t*.py)
  numbers         code/build_numbers_paper_A.py  -> numbers/
  figures         code/analysis/make_figure_data.py, code/make_figures.py -> figures/
  tables          code/make_tables.py -> tables/

Environment variables:
  TITAN_DATA_DIR  directory that contains "Dry Simulations/" (Hartwick et al. 2026,
                  doi:10.5281/zenodo.17246644). Default: <package>/data.
  TITAN_TAM_DIR   directory with the eight atmos_diurnal_*.nc files of Han and Lora (2025),
                  doi:10.5281/zenodo.16791001. Default: $TITAN_DATA_DIR/TAM Simulations_zenodo_16791001.
  WRF_SRC         optional: directory with <tag>/module_*.F of the public WRF source
                  (V3.0, V3.6, V3.9.1.1, v4.8.0). Without it, the recorded transcription is used.
  SKIP_ARCHIVE=1  run only what does not need the archive or the TAM files (framework
                  diagnostics and threshold verification) and stop with a notice.

Usage:
    python code/run_all.py [--keep-build] [--only STEP[,STEP...]]
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
AN = CODE / "analysis"
BUILD = ROOT / "build"
TIMEOUT_S = 600

# (nombre, script, argumentos, necesita archivo/TAM)
STEPS = [
    # source_a_core.py honra SKIP_ARCHIVE por si mismo (solo omite las lecturas del archivo).
    ("source_a_full_absorption", AN / "source_a_core.py",
     ["--lw-boundary", "full_absorption", "--outdir", str(BUILD / "source_a" / "full_absorption")], False),
    ("source_a_kirchhoff", AN / "source_a_core.py",
     ["--lw-boundary", "kirchhoff", "--outdir", str(BUILD / "source_a" / "kirchhoff")], False),
    ("r1_regenerate_source_a", AN / "r1_regenerate_source_a.py", [], True),
    ("r1_threshold_verification", AN / "r1_threshold_verification.py", [], False),
    ("r1_flux_comparison", AN / "r1_flux_comparison.py", [], True),
    ("r2_t1_tam_boundary", AN / "r2_t1_tam_boundary.py", [], True),
    ("r2_t3_mayor3_budget", AN / "r2_t3_mayor3_budget.py", [], True),
    ("r2_t2_energy_partition", AN / "r2_t2_energy_partition.py", [], True),
    ("r2_t4_model_comparison", AN / "r2_t4_model_comparison.py", [], True),
    ("r2_t5_sw_dtreq", AN / "r2_t5_sw_dtreq.py", [], True),
    ("r2_t6_literature_contrasts", AN / "r2_t6_literature_contrasts.py", [], True),
    ("r2_t8_tau_rad", AN / "r2_t8_tau_rad.py", [], True),
    ("r2_t3b_column_energy", AN / "r2_t3b_column_energy.py", [], True),
    # Orden de registro de la etapa r3 (fija el orden de las entradas de la fuente unica).
    ("r3_t4_robustness", AN / "r3_t4_robustness.py", [], True),
    ("r3_t5_tau_rad", AN / "r3_t5_tau_rad.py", [], True),
    ("r3_t1_h33_wrf", AN / "r3_t1_h33_wrf.py", [], True),
    ("r3_t3_flux_comparison_all", AN / "r3_t3_flux_comparison_all.py", [], True),
    ("build_numbers", CODE / "build_numbers_paper_A.py", [], True),
    ("make_figure_data", AN / "make_figure_data.py", [], True),
    ("make_figures", CODE / "make_figures.py", [], True),
    ("make_tables", CODE / "make_tables.py", [], True),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--keep-build", action="store_true", help="do not delete build/ first")
    ap.add_argument("--only", default="", help="comma-separated step names")
    args = ap.parse_args()
    skip_archive = os.environ.get("SKIP_ARCHIVE", "0") == "1"
    only = [s for s in args.only.split(",") if s]
    if not args.keep_build and not only and BUILD.exists():
        shutil.rmtree(BUILD)
    (BUILD / "logs").mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PAPER_A_OUT=str(BUILD / "source_a" / "_v1_outputs"))
    t_all = time.time()
    for name, script, extra, needs_data in STEPS:
        if only and name not in only:
            continue
        if skip_archive and needs_data:
            print(f"[skip] {name}: needs the archive or the TAM files (SKIP_ARCHIVE=1)")
            continue
        t0 = time.time()
        log = BUILD / "logs" / f"{name}.log"
        with log.open("w", encoding="utf-8") as fh:
            try:
                rc = subprocess.run([sys.executable, str(script), *extra], cwd=ROOT, env=env, stdout=fh,
                                    stderr=subprocess.STDOUT, timeout=TIMEOUT_S).returncode
            except subprocess.TimeoutExpired:
                rc = "timeout"
        print(f"[{'ok' if rc == 0 else 'FAIL'}] {name} ({time.time() - t0:.0f} s) -> {log.relative_to(ROOT)}")
        if rc != 0:
            print(log.read_text(encoding="utf-8", errors="replace")[-3000:])
            return 1
    if skip_archive:
        print("NOTICE: SKIP_ARCHIVE=1. Only the framework diagnostics that do not read the archive of "
              "Hartwick et al. (2026) or the TAM output of Han and Lora (2025) were run. The single source of "
              "numbers, the figures and the tables need both; set TITAN_DATA_DIR (and TITAN_TAM_DIR) and run again.")
    print(f"total {time.time() - t_all:.0f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
