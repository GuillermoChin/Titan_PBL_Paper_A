#!/usr/bin/env python
"""Audit every number quoted in the Paper A manuscript against its CSV source.

Run this after `make_paper_A_figures.py` and before submitting. Each row
declares what the manuscript says, where the value must come from, and the
tolerance at which the two are considered to agree. The script exits non-zero
if any number in the text has drifted away from the output that produced it.

    python data/audit_paper_A_numbers.py

Why this exists: the manuscript quotes about twenty numbers, all of them
generated. Every regeneration of the figures is an opportunity for the text to
fall out of step with them, and a reviewer who recomputes one number and finds
a mismatch has no reason to trust the other nineteen.
"""
from __future__ import annotations
import sys
from pathlib import Path
# El subconjunto de titan_pbl vive junto a este script
sys.path.insert(0, str(Path(__file__).resolve().parent / "titan_pbl_subset"))
import numpy as np
import pandas as pd

# Las rutas se derivan de la ubicacion del script para que el paquete de
# reproducibilidad funcione desde cualquier directorio.
DATA_DIR = Path(__file__).resolve().parent
REPO_ROOT = DATA_DIR.parent
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else DATA_DIR / "paper_A"

thresholds = pd.read_csv(OUT / "table_thresholds.csv")
mixing = pd.read_csv(OUT / "table_mixing_times.csv")
table3 = pd.read_csv(OUT / "table3_validation.csv")
table4 = pd.read_csv(OUT / "table4_contrasts.csv")
table2 = pd.read_csv(OUT / "table2_parameter_sweep.csv")



def thr(name: str) -> float:
    return float(thresholds.loc[thresholds["quantity"].str.contains(name), "contrast_K"].iloc[0])


def mix(contrast: float, column: str) -> float:
    row = mixing.loc[np.isclose(mixing["contrast_K"], contrast)]
    return float(row[column].iloc[0])


CHECKS = [
    # (etiqueta, valor declarado en el texto, valor del CSV, tolerancia, ubicacion)
    ("Umbral analitico", 1.20, thr("Analytic"), 0.005, "§2.2 L156, §3.3 tabla"),
    ("Umbral de calentamiento", 2.15, thr("Numerical heating"), 0.005, "§3.3 tabla"),
    ("Umbral de desestabilizacion", 2.53, thr("Numerical destabilisation"), 0.005, "§3.3 tabla, §3.5, §4.4"),
    ("Tasa a 3.0 K", 0.0183, mix(3.0, "destabilisation_K_per_sol"), 5e-5, "§3.5 tabla"),
    ("Tasa a 3.5 K", 0.0382, mix(3.5, "destabilisation_K_per_sol"), 5e-5, "§3.5 tabla"),
    ("Tasa a 4.5 K", 0.0790, mix(4.5, "destabilisation_K_per_sol"), 5e-5, "§3.5 tabla"),
    ("Mezcla a 3.0 K", 35.7, mix(3.0, "mixing_time_sols"), 0.05, "§3.5 tabla"),
    ("Mezcla a 3.5 K", 17.1, mix(3.5, "mixing_time_sols"), 0.05, "§3.5 tabla"),
    ("Mezcla a 4.5 K", 8.3, mix(4.5, "mixing_time_sols"), 0.05, "§3.5 tabla"),
    ("Mezcla a 3.0 K, redondeo entero", 36.0, round(mix(3.0, "mixing_time_sols")), 0.5, "§4.4, §6 punto 5"),
    ("Concordancia LW capa baja", 0.90, float(table3["difference_percent"].iloc[0]), 0.01, "§3.6 (nueva), §6"),
    ("Contraste maximo suprimido, Tsol 3", 0.67, float(table4["contrast_suppressed_max_K"].iloc[1]), 0.005, "§4.5, punto abierto 1"),
    ("Ascenso de piel, Tsol 5", 1.22, float(table4["delta_Ts_between_runs_max_K"].iloc[3]), 0.005, "punto abierto 1"),
    ("HFX suprimido pico, Tsol 3", 1.06, float(table4["HFX_suppressed_peak_W_m2"].iloc[1]), 0.005, "§4.5"),
    ("HFX baseline pico, Tsol 3", 1.79, float(table4["HFX_baseline_peak_W_m2"].iloc[1]), 0.005, "§4.5"),
    ("Reduccion de flujo", 41.0, (1 - float(table4["HFX_ratio"].iloc[1])) * 100, 0.6, "§4.5, §6 punto 6"),
    ("Deriva PBL suprimida", 55.3, (float(table4["PBL_suppressed_max_m"].iloc[3]) / float(table4["PBL_suppressed_max_m"].iloc[0]) - 1) * 100, 0.1, "§4.5 model drift"),
    ("Maximo calentamiento del barrido", 0.021, float(table2["max_heating_K_per_sol"].max()), 0.0005, "§3.4"),
]

print(f"{'cantidad':<38} {'texto':>10} {'CSV':>12} {'ok':>4}  ubicacion")
print("-" * 100)
failures = []
for label, declared, computed, tol, where in CHECKS:
    ok = abs(declared - computed) <= tol
    if not ok:
        failures.append((label, declared, computed, where))
    print(f"{label:<38} {declared:>10.4f} {computed:>12.6f} {'si' if ok else 'NO':>4}  {where}")

# Delta-Theta viene del perfil HASI, no de un CSV: se recalcula aparte.
print("\n--- verificacion aparte: Delta-Theta 0-4 km desde el perfil HASI ---")
sys.path.insert(0, str(REPO_ROOT))
from titan_pbl.config import load_config
from titan_pbl.core import build_grid, initialize_state

cfg = load_config("huygens_baseline")
state = initialize_state(cfg, build_grid(cfg))
_, z = state.geopotential_heights()
theta = state.potential_temperature()
top = int(np.argmin(np.abs(z - 4000.0)))
delta = float(theta[top] - theta[0])
declared = 0.652
ok = abs(delta - declared) <= 0.002
print(f"{'Delta-Theta 0-4 km':<38} {declared:>10.4f} {delta:>12.6f} {'si' if ok else 'NO':>4}  §3.5, §4.5")
if not ok:
    failures.append(("Delta-Theta 0-4 km", declared, delta, "§3.5"))

print("\n" + "=" * 100)
if failures:
    print(f"{len(failures)} DISCREPANCIAS:")
    for label, d, c, where in failures:
        print(f"  {label}: texto {d}, CSV {c:.6f}  ({where})")
    sys.exit(1)

print(f"Los {len(CHECKS) + 1} numeros del manuscrito coinciden con los CSV.")
sys.exit(0)
