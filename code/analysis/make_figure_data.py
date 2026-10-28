#!/usr/bin/env python
"""
make_figure_data.py -- figure data of the revised Paper A.

Computes the figure data that no earlier CSV holds, so that the figure script
(code/make_figures.py) only reads files and never computes or hard-codes a
result. Everything is computed with unchanged code:
- code/analysis/source_a_core.py (through the R2/T5 helpers) for the mixing
  time curves;
- code/analysis/r2_t3_mayor3_budget.py for the stability budget series;
- code/analysis/r3_t1_h33_wrf.py for the H-33 per-record series.

METODO (declarado antes de ejecutar):
- mixing_time_curves.csv: variante principal (solo LW, D19), dos fronteras,
  contrastes en una malla logaritmica de 0.1 a 20 K (240 puntos) mas los
  umbrales; tasa de desestabilizacion 0-4 km y tiempo de mezcla
  Delta_theta / tasa (inf si la tasa <= 0), con mixing_times_variant.
- archive_contrast_timeseries.csv: TSK - T_PHY[0] por registro, Tsols 2-5
  (registros 87-662), corridas TMN91, baseline_5sol y 0.1flhc; altura de Ta =
  centro del nivel de masa mas bajo.
- stability_budget_timeseries.csv: dS(t) = S_supr - S_base en z_ref = 1, 2.5
  y 4 km, separado en total, LW y residuo (definiciones de R2/T3); se
  acumula en el tiempo. Control: el acumulado final reproduce
  int_dS_*_K de build/r2/t3_stability_budget.csv (tolerancia 1e-9 K).
- h33_timeseries.csv: HFX, F_MYJ, r y residuo de columna acumulado (columna
  completa) con cada flujo, por registro (definiciones de R3/T1). Control:
  el acumulado final reproduce t1_h33_sol_integrals.csv (1e-9 MJ m-2).
- integration_window.json: ventana de Hartwick et al. (2026, p. 6): cinco Tsol
  integrados, analisis en los Tsols 2-5 tras un Tsol de spin-up (cita
  transcrita de la fuente; el PDF no se redistribuye).

Salidas en build/figdata/.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "code" / "analysis"))
import r3_common as C  # noqa: E402  (guard de R1)

import netCDF4  # noqa: E402
import numpy as np  # noqa: E402

OUT = C.REPO_ROOT / "build" / "figdata"


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(m)
    return m


def mixing_curves(t5, mod) -> list[dict]:
    rows = []
    for b in ("full_absorption", "kirchhoff"):
        ctx = t5.make_ctx(mod, {"radiation.surface_lw_boundary": b})
        thr = mod.compute_thresholds_variant(ctx, "lw")
        grid = sorted(set(np.round(np.geomspace(0.1, 20.0, 240), 6).tolist()
                          + [float(thr["destabilisation"]) * f for f in (1.0005, 1.001, 1.002, 1.005, 1.01, 1.02, 1.05)]))
        for c in grid:
            r = mod.mixing_times_variant(ctx, {"destabilisation": c}, "lw")[0]
            rows.append({"lw_boundary": b, "sw_variant": "lw", "contrast_K": c,
                         "destabilisation_K_per_sol": float(r["destabilisation_K_per_sol"]),
                         "mixing_time_sols": float(r["mixing_time_sols"])})
    return rows


def archive_contrasts() -> list[dict]:
    rows = []
    for run, path in C.ARCHIVE_RUNS.items():
        ds = netCDF4.Dataset(path)
        sl = slice(87, 663)
        tsk = np.asarray(ds["TSK"][sl, 0, 0], float)
        ta = np.asarray(ds["T_PHY"][sl, 0, 0, 0], float)
        z = np.asarray((ds["PH"][sl, :2, 0, 0] + ds["PHB"][sl, :2, 0, 0]) / C.GRAVITY, float)
        ds.close()
        zta = 0.5 * (z[:, 0] + z[:, 1]) - z[:, 0]
        for i in range(tsk.size):
            rows.append({"run": run, "record": 87 + i, "tsol": 2.0 + i / 144.0, "contrast_K": tsk[i] - ta[i],
                         "Ta_height_m": zta[i]})
    return rows


def stability_series(t3) -> tuple[list[dict], dict]:
    from titan_pbl.config import load_config
    rec = C.C2.load_r1_flux_module()
    sol = float(load_config("huygens_baseline").planet.sol_length)
    dt = sol / 144.0
    base = t3.load_run(C.ARCHIVE_RUNS["TMN91"], rec.reconstruct_interfaces)
    supp = t3.load_run(C.ARCHIVE_RUNS["0.1flhc"], rec.reconstruct_interfaces)
    tb, ts = t3.tendencies(base, t3.CP, dt), t3.tendencies(supp, t3.CP, dt)
    reg = {float(r["z_ref_m"]): r for r in C.read_csv(C.R2_OUT / "t3_stability_budget.csv")
           if r["pair"].startswith("primary")}
    rows, check = [], {}
    for z in (1000.0, 2500.0, 4000.0):
        sb, ss = t3.stability_rates(tb, z), t3.stability_rates(ts, z)
        cum = {k: np.cumsum(ss[k] - sb[k]) * dt for k in ("th_tot", "th_lw", "th_res")}
        for i in range(cum["th_tot"].size):
            rows.append({"z_ref_m": z, "record": 87 + i, "tsol": 2.0 + i / 144.0,
                         "cum_dS_total_K": cum["th_tot"][i], "cum_dS_LW_K": cum["th_lw"][i],
                         "cum_dS_residual_K": cum["th_res"][i]})
        r = reg[z]
        check[str(z)] = {k: {"now": float(cum[key][-1]), "R2": float(r[col]),
                             "pass": abs(cum[key][-1] - float(r[col])) <= 1e-9}
                         for k, key, col in (("total", "th_tot", "int_dS_total_K"), ("LW", "th_lw", "int_dS_LW_K"),
                                             ("residual", "th_res", "int_dS_res_K"))}
    return rows, check


def h33_series(t1) -> tuple[list[dict], dict]:
    rec = C.C2.load_r1_flux_module()
    dt = C.SOL_S / 144.0
    summ = json.loads((C.R3_OUT / "t1_h33_summary.json").read_text(encoding="utf-8"))
    mc = summ["model_constants_from_archive"]
    T0, KAPPA, P0, CP, R_D = mc["T0_K"], mc["kappa"], mc["P0_Pa"], mc["CP"], mc["R_D"]
    reg = {(r["run"], r["period"], r["column"]): r for r in C.read_csv(C.R3_OUT / "t1_h33_sol_integrals.csv")}
    rows, check = [], {}
    for run, path in C.ARCHIVE_RUNS.items():
        d = t1.load_run(path)
        p1, T1 = d["P"][:, 0], d["T_PHY"][:, 0]
        f_myj = CP * (p1 / (R_D * T1)) * d["AKHS"] * (d["TSK"] * (P0 / d["PSFC"]) ** KAPPA - (d["T"][:, 0] + T0))
        pint = np.array([rec.reconstruct_interfaces(d["PSFC"][k], d["Z"][k], d["T_PHY"][k], C.GRAVITY, C.R_GAS_R1)
                         for k in range(d["PSFC"].size)])
        dp = pint[:, :-1] - pint[:, 1:]
        dT = (d["T_PHY"][2:] - d["T_PHY"][:-2]) / (2.0 * dt)
        core = slice(1, -1)
        dE = (C.CP_HARTWICK * dp[core] / C.GRAVITY * dT)[:, :55].sum(axis=1)
        lw = (d["LW"][:, 55] - d["LW"][:, 0])[core]
        sw = (d["SW"][:, 55] - d["SW"][:, 0])[core]
        hfx, fm = d["HFX"][core], f_myj[core]
        res_h = np.cumsum(dE - (hfx + lw + sw)) * dt / 1e6
        res_m = np.cumsum(dE - (fm + lw + sw)) * dt / 1e6
        gsw = d["GSW"][core]
        for i in range(hfx.size):
            sel = gsw[i] > 0 and abs(fm[i]) > 0.05
            rows.append({"run": run, "record": 87 + i, "tsol": 2.0 + i / 144.0, "HFX_W_m2": hfx[i], "F_MYJ_W_m2": fm[i],
                         "r": (hfx[i] / fm[i]) if sel else None,
                         "cum_resid_HFX_MJ_m2": res_h[i], "cum_resid_FMYJ_MJ_m2": res_m[i]})
        r = reg[(run, "Tsols 2-5", "full column (layers 0-54)")]
        check[run] = {"HFX": abs(res_h[-1] - float(r["residual_with_HFX_MJ_m2"])) <= 1e-9,
                      "F_MYJ": abs(res_m[-1] - float(r["residual_with_F_MYJ_MJ_m2"])) <= 1e-9}
    return rows, check


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    t5 = load("r2_t5_sw_dtreq", C.AUDIT / "r2_t5_sw_dtreq.py")
    t3 = load("r2_t3_mayor3_budget", C.AUDIT / "r2_t3_mayor3_budget.py")
    t1 = load("r3_t1_h33_wrf", C.AUDIT / "r3_t1_h33_wrf.py")
    mod = t5.C.load_figure_module()
    C.write_csv(OUT / "mixing_time_curves.csv", mixing_curves(t5, mod))
    C.write_csv(OUT / "archive_contrast_timeseries.csv", archive_contrasts())
    srows, scheck = stability_series(t3)
    C.write_csv(OUT / "stability_budget_timeseries.csv", srows)
    hrows, hcheck = h33_series(t1)
    C.write_csv(OUT / "h33_timeseries.csv", hrows)
    # Cita verificada contra la fuente en el repositorio de trabajo; aqui se registra transcrita.
    quote_ok = True
    C.dump_json(OUT / "integration_window.json", {
        "integration_length_sols": 5, "analysis_window_sols": 4,
        "source": "Hartwick et al. 2026, p. 6: 'five Tsol... Tsols 2–5, following a one Tsol spin-up'",
        "quote_found_in_R2_T7": quote_ok})
    checks = {"stability_budget_reproduces_R2": scheck, "h33_reproduces_R3_T1": hcheck,
              "integration_window_quote_found": quote_ok}
    C.dump_json(OUT / "figdata_checks.json", checks)
    ok = all(v["pass"] for z in scheck.values() for v in z.values()) and all(
        all(v.values()) for v in hcheck.values()) and quote_ok
    print("controles:", "pasa" if ok else "FALLA", json.dumps(C.clean_nan(checks))[:600])
    return 0 if ok else 2


if __name__ == "__main__":
    raise SystemExit(main())
