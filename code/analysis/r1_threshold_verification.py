#!/usr/bin/env python
"""
r1_threshold_verification.py -- stage r1, analysis T4.

Verifies the exact zero of the layer-0 heating as a function of the
surface-to-air contrast, against a numerical root found with ``brentq``, for
both surface longwave boundaries and a sweep of radiative parameters.

Derivacion (Seccion 2.3 del manuscrito, comprobada contra gray.py). Con
t = exp(-dtau0), B = sigma*Ta^4 y los flujos del esquema,

    net[1] - net[0] = (1 - t) * [F_down(1) + F_up(0) - 2B] + S,

donde S es la divergencia de onda corta de la capa 0, que no depende de Ts, y
F_down(1) tampoco. El cero exacto es

    O1: eps*sigma*Ts*^4 = 2B - F_down(1) - S/(1 - t)
    K : eps*sigma*Ts*^4 = 2B - F_down(1) - S/(1 - t) - (1 - eps)*F_down(0),
        con F_down(0) = F_down(1)*t + B*(1 - t),

y Delta_T* = Ts* - Ta. Formas linealizadas, solo LW, con D = B - F_down(1):

    O1: Delta_T* ~ [(1 - eps)*B + D] / (4*eps*sigma*Ta^3)
    K : Delta_T* ~ (2 - eps)*D / (4*eps*sigma*Ta^3)       (supone t ~ 1)

Como informacion adicional se registra la linealizacion K sin suponer t ~ 1,
(1 + (1 - eps)*t)*D / (4*eps*sigma*Ta^3).

CRITERIOS DECLARADOS antes de ejecutar:
- |Delta_T_num - Delta_T_exact| <= 1e-6 K en TODAS las combinaciones con raiz
  (si no: condicion de paro P3, no se ajusta nada).
- Isotermo con kirchhoff: Delta_T_exact = 0 con tolerancia 1e-9 K. Se aplica a
  las filas solo LW: con onda corta, S != 0 desplaza el cero por construccion.
- El error de linealizacion (Delta_T_lin - Delta_T_exact) se reporta, no se
  juzga.
- Si brentq no encuentra cambio de signo en [-5, 15] K, la fila se marca
  "sin raiz en el intervalo"; si el cero exacto cae dentro del intervalo, eso
  cuenta como discrepancia.

Salidas: build/r1/t4_thresholds.csv y build/r1/t4_summary.json.
"""

from __future__ import annotations

import csv
import itertools
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import numpy as np  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT / "code" / "titan_pbl_subset") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "code" / "titan_pbl_subset"))

import titan_pbl  # noqa: E402
from titan_pbl.config import load_config  # noqa: E402
from titan_pbl.config.physical_constants import STEFAN_BOLTZMANN  # noqa: E402
from titan_pbl.core import build_grid, initialize_state  # noqa: E402
from titan_pbl.schemes.radiation import GrayRadiation  # noqa: E402

if (REPO_ROOT / "code" / "titan_pbl_subset") not in Path(titan_pbl.__file__).resolve().parents:
    raise ImportError(f"titan_pbl imported from {titan_pbl.__file__}, not from code/titan_pbl_subset.")

OUTPUT_DIR = REPO_ROOT / "build" / "r1"
SIGMA = STEFAN_BOLTZMANN
T_ISO = 93.65
BRACKET = (-5.0, 15.0)
TOL_NUM_EXACT = 1.0e-6
TOL_ISO_KIRCHHOFF = 1.0e-9

PROFILES = [("hasi", 0.0), ("hasi", 2.0), ("hasi", 8.0), ("isothermal", None)]
EMISSIVITIES = (0.90, 0.95, 1.00)
TAUS = (5.0, 10.0, 20.0, 30.0)
WEIGHTS = (0.15, 0.50)
RADIATION = ("lw", "lw+sw_noon")
BOUNDARIES = ("full_absorption", "kirchhoff")


def build(profile, i_ir, eps, tau, f, radiation, boundary):
    overrides = {
        "radiation.surface_emissivity": eps,
        "radiation.optical_depth_ir": tau,
        "radiation.ir_weighting_coefficient": f,
        "radiation.surface_lw_boundary": boundary,
    }
    if profile == "isothermal":
        overrides["atmosphere.initialization.profile_source"] = "isothermal"
        overrides["atmosphere.initialization.isothermal_temperature"] = T_ISO
        overrides["radiation.ir_top_flux"] = SIGMA * T_ISO**4
    else:
        overrides["radiation.ir_top_flux"] = i_ir
    if radiation == "lw":
        overrides["radiation.solar_constant"] = 0.0
    config = load_config("huygens_baseline", overrides=overrides)
    state = initialize_state(config, build_grid(config))
    # Mediodia del Tsol 1, igual que load_framework del script de figuras.
    state.time = (12.0 - config.site.start_local_time) / 24.0 * config.planet.sol_length
    scheme = GrayRadiation(config)
    scheme.initialize(config, state)
    return config, state, scheme


def evaluate(profile, i_ir, eps, tau, f, radiation, boundary) -> dict:
    config, state, scheme = build(profile, i_ir, eps, tau, f, radiation, boundary)
    ta = float(state.atmospheric_temperature[0])

    def layer0(delta: float) -> float:
        state.skin_temperature = ta + delta
        return float(scheme.heating_rate(state)[0])

    # Raiz numerica, imponiendo Delta_T solo via Ts (como en la Fig. 2).
    lo, hi = BRACKET
    f_lo, f_hi = layer0(lo), layer0(hi)
    if np.sign(f_lo) == np.sign(f_hi):
        dt_num = float("nan")
        note = "sin raiz en el intervalo"
    else:
        dt_num = float(brentq(layer0, lo, hi, xtol=1.0e-13, rtol=4.0 * np.finfo(float).eps, maxiter=500))
        note = ""

    # Cantidades del cero exacto. Los flujos descendentes y la onda corta no
    # dependen de Ts, de modo que basta un estado cualquiera.
    state.skin_temperature = ta
    lw_down, _ = scheme.longwave_fluxes(state)
    sw_down, sw_up = scheme.shortwave_fluxes(state)
    sw_net = sw_down - sw_up
    tau_ir = scheme.infrared_optical_depth(state)
    d_tau0 = float(tau_ir[0] - tau_ir[1])
    t = float(np.exp(-d_tau0))
    b = SIGMA * ta**4
    d1 = float(lw_down[1])
    d0 = float(lw_down[0])
    s = float(sw_net[1] - sw_net[0])
    rhs = 2.0 * b - d1 - s / (1.0 - t)
    if boundary == "kirchhoff":
        rhs -= (1.0 - eps) * d0
    dt_exact = (rhs / (eps * SIGMA)) ** 0.25 - ta if rhs > 0.0 else float("nan")

    big_d = b - d1
    denom = 4.0 * eps * SIGMA * ta**3
    if radiation == "lw":
        if boundary == "full_absorption":
            dt_lin = ((1.0 - eps) * b + big_d) / denom
            dt_lin_t = dt_lin
        else:
            dt_lin = (2.0 - eps) * big_d / denom
            dt_lin_t = (1.0 + (1.0 - eps) * t) * big_d / denom
    else:
        dt_lin = dt_lin_t = float("nan")

    diff = abs(dt_num - dt_exact) if np.isfinite(dt_num) and np.isfinite(dt_exact) else float("nan")
    if not np.isfinite(dt_num) and np.isfinite(dt_exact) and lo <= dt_exact <= hi:
        note = "DISCREPANCIA: brentq sin raiz pero el cero exacto cae en el intervalo"

    return {
        "profile": profile,
        "I_IR_W_m2": float(config.radiation.ir_top_flux),
        "emissivity": eps,
        "tau_IR": tau,
        "f": f,
        "radiation": radiation,
        "boundary": boundary,
        "Ta_K": ta,
        "d_tau0": d_tau0,
        "F_down_0_W_m2": d0,
        "F_down_1_W_m2": d1,
        "D_W_m2": big_d,
        "S_sw_divergence_W_m2": s,
        "dT_num_K": dt_num,
        "dT_exact_K": dt_exact,
        "abs_num_minus_exact_K": diff,
        "dT_lin_K": dt_lin,
        "lin_minus_exact_K": dt_lin - dt_exact if np.isfinite(dt_lin) else float("nan"),
        "dT_lin_K_kirchhoff_with_t": dt_lin_t,
        "note": note,
    }


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for (profile, i_ir), eps, tau, f, radiation, boundary in itertools.product(
        PROFILES, EMISSIVITIES, TAUS, WEIGHTS, RADIATION, BOUNDARIES
    ):
        rows.append(evaluate(profile, i_ir, eps, tau, f, radiation, boundary))

    with (OUTPUT_DIR / "t4_thresholds.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    with_root = [r for r in rows if np.isfinite(r["abs_num_minus_exact_K"])]
    discrepancies = [r for r in rows if r["note"].startswith("DISCREPANCIA")]
    max_diff = max(r["abs_num_minus_exact_K"] for r in with_root)
    iso_k = [r for r in rows if r["profile"] == "isothermal" and r["boundary"] == "kirchhoff" and r["radiation"] == "lw"]
    iso_k_max = max(abs(r["dT_exact_K"]) for r in iso_k)
    lin_rows = [r for r in rows if np.isfinite(r["lin_minus_exact_K"])]

    def lin_stats(boundary):
        sel = [r["lin_minus_exact_K"] for r in lin_rows if r["boundary"] == boundary]
        return {"n": len(sel), "min": min(sel), "max": max(sel), "mean": float(np.mean(sel))}

    reference = {
        (r["boundary"], r["radiation"]): r for r in rows
        if r["profile"] == "hasi" and r["I_IR_W_m2"] == 2.0 and r["emissivity"] == 0.95
        and r["tau_IR"] == 10.0 and r["f"] == 0.15
    }
    summary = {
        "task": "R1-T4",
        "n_combinations": len(rows),
        "n_with_root": len(with_root),
        "n_without_root_in_bracket": len(rows) - len(with_root),
        "max_abs_num_minus_exact_K": max_diff,
        "criterion_num_exact_pass": bool(max_diff <= TOL_NUM_EXACT and not discrepancies),
        "discrepancies": [{k: r[k] for k in ("profile", "I_IR_W_m2", "emissivity", "tau_IR", "f", "radiation", "boundary")} for r in discrepancies],
        "isothermal_kirchhoff_lw_max_abs_dT_exact_K": iso_k_max,
        "criterion_isothermal_kirchhoff_pass": bool(iso_k_max <= TOL_ISO_KIRCHHOFF),
        "linearisation_error_K": {b: lin_stats(b) for b in BOUNDARIES},
        "preset_reference_HASI_IIR2_eps095_tau10_f015": {
            f"{b}|{rad}": {k: reference[(b, rad)][k] for k in ("dT_num_K", "dT_exact_K", "dT_lin_K", "dT_lin_K_kirchhoff_with_t")}
            for (b, rad) in reference
        },
    }
    (OUTPUT_DIR / "t4_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=1))
    return 0 if summary["criterion_num_exact_pass"] and summary["criterion_isothermal_kirchhoff_pass"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
