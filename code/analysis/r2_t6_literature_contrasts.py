#!/usr/bin/env python
"""
r2_t6_literature_contrasts.py -- surface temperature quantities of the literature (stage r2).

Package version. In the working repository this step read the original
submission and the source PDFs to verify every quoted value; neither is
redistributed here. The literature values, with their page and verification
status, are the curated input inputs/literature_quantities.csv. The two
contrasts read from the archive of Hartwick et al. (2026) are not typed: they
are the maxima over sols 2-5 computed by r1_regenerate_source_a.py
(numbers_r1.json, keys archive_contrast_max_sols2to5_*).

Salida: build/r2/t6_literature_contrasts.csv (columnas que usa
build_numbers_paper_A.py: id, estudio, contraste_K, naturaleza, altura_Ta).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import r2_common as C  # noqa: E402

ARCHIVE_ROWS = [
    ("archive_baseline", "Archive of Hartwick et al. 2026 (baseline TMN91)", "archive_contrast_max_sols2to5_baseline"),
    ("archive_suppressed", "Archive of Hartwick et al. 2026 (0.1flhc)", "archive_contrast_max_sols2to5_suppressed"),
]


def main() -> int:
    with (C.INPUTS / "literature_quantities.csv").open(newline="", encoding="utf-8") as fh:
        lit = list(csv.DictReader(fh))
    rows = []
    for r in lit:
        rows.append({"id": r["id"], "estudio": r["study"], "contraste_K": r["value_K"], "naturaleza": r["nature"],
                     "altura_Ta": r["Ta_height"], "estado_cita": r["citation_status"], "fuente_pagina": r["source_page"],
                     "nota": r["note"]})
    for rid, study, key in ARCHIVE_ROWS:
        rows.append({"id": rid, "estudio": study, "contraste_K": repr(float(C.number_r1(key))),
                     "naturaleza": "surface-to-air contrast (TSK - T_PHY[0])", "altura_Ta": "lowest mass level (~2.7 m)",
                     "estado_cita": "archive reading (Zenodo 17246644)", "fuente_pagina": "",
                     "nota": "read from the archive, not from the text; maximum over sols 2-5"})
    order = ["v0_generic_low", "v0_generic_high", "hanlora_amp_low", "hanlora_amp_high", "hartwick_gt1K",
             "archive_baseline", "archive_suppressed", "archive_dTs", "charnay_theta_10m", "charnay_Ts_amp",
             "tokano2006_qual"]
    rows.sort(key=lambda r: order.index(r["id"]))
    C.R2_OUT.mkdir(parents=True, exist_ok=True)
    with (C.R2_OUT / "t6_literature_contrasts.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"{len(rows)} filas -> {C.rel(C.R2_OUT / 't6_literature_contrasts.csv')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
