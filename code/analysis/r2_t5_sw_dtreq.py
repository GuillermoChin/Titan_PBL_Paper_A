#!/usr/bin/env python
"""
r2_t5_sw_dtreq.py -- stage r2, analysis T5.

Shortwave variants, required contrast Delta T_req(N), real diurnal cycle,
one-at-a-time sensitivity, most favourable corner and HASI uncertainty.
Everything is computed with the functions of code/analysis/source_a_core.py,
imported unchanged.

---------------------------------------------------------------------------
CRITERIOS Y METODO (declarados antes de ejecutar)
---------------------------------------------------------------------------
1. Matriz principal, 2 fronteras x 3 variantes SW (lw, noon, diurnal):
   - umbral de calentamiento de la capa 0 y umbral de desestabilizacion
     0-4 km, con compute_thresholds_variant (misma definicion que
     compute_thresholds);
   - tiempo de mezcla a los contrastes {maximo del archivo (clave
     archive_contrast_max_sols2to5_suppressed de numbers_r1.json), cada
     contraste de la literatura (T6), 3.0, 4.5 K};
   - Delta T_req(N), N en {4, 5} sols: raiz con brentq de mixing_time(dT) = N
     con el metodo de mezcla existente (mixing_times_variant, sin modificar).
   Implementacion de mixing_time(dT): mixing_times_variant(ctx,
   {"destabilisation": dT}, variant)[0]; su primera fila evalua el metodo en
   dT. Como el metodo devuelve inf si la tasa es <= 0, la raiz se busca en la
   forma equivalente rate(dT) = Delta_theta / N (rate = la tasa que el mismo
   metodo devuelve), en [-5, 20] K, xtol = 1e-6; despues se verifica que
   mixing_time(dT_req) = N con el metodo, a 1e-6 relativo.
   Contrastes de la literatura (T6.1, extraidos de v0 antes de ejecutar):
     - "of order 0.1 to 2 K" (v0 par. 70, Sec. 3.3; estudios de la Sec. 4):
       extremos 0.1 y 2.0 K;
     - Han & Lora 2025, "surface diurnal amplitudes from 0.12 to 0.21 K"
       (v0 par. 97, Sec. 4.3): 0.12 y 0.21 K (amplitud de Ts, no contraste);
     - Hartwick et al. 2026, "warms by more than one kelvin" (v0 par. 100,
       101 y 112): 1.0 K (es el aumento de Ts entre corridas, segun v0 par. 103);
     - archivo, "0.22 K ... and 0.67 K" (v0 par. 103), sustituidos por los
       valores D4 de R1: 0.231 y 0.686 K.
2. D12: variante principal = la que minimiza Delta T_req(4) con kirchhoff;
   se registra como headline_sw_variant en numbers_r2.json.
3. Ciclo diurno real: dT(t) = TSK - T_PHY[0] del archivo (TMN91 y 0.1flhc,
   Tsols 2-5) con la SW de la hora local de cada registro (Times "_HH:MM"),
   sobre el estado HASI inicial congelado: state.time = ((LT - 9.5) mod 24) /
   24 * Tsol (Tsol 1 de titan_pbl); tasa = destabilisation_rate (existente).
   Tasa media por sol y mezcla = Delta_theta / tasa media; si <= 0, "nunca".
   Ta del archivo en el nivel de masa mas bajo (~2.7 m) frente a la capa 0 de
   titan_pbl (5.38 m). Las dos fronteras.
4. Sensibilidad una variable a la vez, alrededor del preset, con kirchhoff y
   la variante principal: alpha {0.1..0.5}, eps {0.90, 0.95, 1.00}, tau_IR
   {5, 10, 20, 30}, f {0.15, 0.50}, I_IR {0, 2, 8}, factor S0 {0.9, 1.0, 1.1}.
   Salidas: los dos umbrales, la mezcla a 3 K y Delta T_req(4).
   Esquina mas favorable: para cada factor, el valor que minimiza
   Delta T_req(4) en el barrido; se evalua la combinacion de todos.
   Limitacion declarada: alpha solo actua aqui por la SW reflejada hacia la
   atmosfera; su efecto principal sobre Ts (balance superficial) queda fuera,
   porque los contrastes se toman de corridas publicadas.
5. Incertidumbre de Delta_theta_HASI (0.652 K, 0-4 km): columnas del L4 sin
   incertidumbre (etiqueta PDS3 del producto); Fulchignoni et al. 2005, Fig. 5
   (p. 787): +-0.25 K entre 60 y 110 K. Propagacion: errores independientes,
   sigma_dtheta = sqrt(s1^2 + s2^2) con s_k = 0.25 K * (theta_k / T_k); error
   sistematico: se cancela. La mezcla es proporcional a Delta_theta; se da
   ademas Delta T_req(4) con Delta_theta +- sigma.

Salidas: build/r2/t5_matrix.csv, t5_dTreq.csv, t5_real_diurnal.csv,
t5_sensitivity_oat.csv, t5_corner.json, t5_hasi_uncertainty.json y claves en
build/r2/numbers_r2.json.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT / "code" / "titan_pbl_subset") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "code" / "titan_pbl_subset"))

import titan_pbl  # noqa: E402

if (REPO_ROOT / "code" / "titan_pbl_subset") not in Path(titan_pbl.__file__).resolve().parents:
    raise ImportError(f"titan_pbl imported from {titan_pbl.__file__}, not from code/titan_pbl_subset.")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import netCDF4  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

import r2_common as C  # noqa: E402

BOUNDARIES = ("full_absorption", "kirchhoff")
NS = (4, 5)
SIGMA_T_HASI = 0.25
SOL_INDICES = {2: (87, 231), 3: (231, 375), 4: (375, 519), 5: (519, 663)}
OAT = {
    "radiation.surface_albedo": (0.1, 0.2, 0.3, 0.4, 0.5),
    "radiation.surface_emissivity": (0.90, 0.95, 1.00),
    "radiation.optical_depth_ir": (5.0, 10.0, 20.0, 30.0),
    "radiation.ir_weighting_coefficient": (0.15, 0.50),
    "radiation.ir_top_flux": (0.0, 2.0, 8.0),
    "S0_factor": (0.9, 1.0, 1.1),
}
PRESET = {"radiation.surface_albedo": 0.3, "radiation.surface_emissivity": 0.95, "radiation.optical_depth_ir": 10.0,
          "radiation.ir_weighting_coefficient": 0.15, "radiation.ir_top_flux": 2.0, "S0_factor": 1.0}


def literature_contrasts() -> list[tuple[str, float, str]]:
    return [
        ("archivo_max_D4_suprimida", float(C.number_r1("archive_contrast_max_sols2to5_suppressed")),
         "numbers_r1: archive_contrast_max_sols2to5_suppressed (v0 par. 103 imprimio 0.67)"),
        ("archivo_D4_baseline", float(C.number_r1("archive_contrast_max_sols2to5_baseline")),
         "numbers_r1: archive_contrast_max_sols2to5_baseline (v0 par. 103 imprimio 0.22)"),
        ("v0_rango_0.1K", 0.1, "v0 par. 70 (Sec. 3.3): 'of order 0.1 to 2 K'"),
        ("v0_rango_2K", 2.0, "v0 par. 70 (Sec. 3.3): 'of order 0.1 to 2 K'"),
        ("HanLora_amp_0.12K", 0.12, "v0 par. 97 (Sec. 4.3): amplitud diurna de Ts, no contraste"),
        ("HanLora_amp_0.21K", 0.21, "v0 par. 97 (Sec. 4.3): amplitud diurna de Ts, no contraste"),
        ("Hartwick_mas_de_1K", 1.0, "v0 par. 100-101, 112: aumento de Ts entre corridas, no contraste"),
        ("fijo_3K", 3.0, "fixed value"),
        ("fijo_4.5K", 4.5, "fixed value"),
    ]


def make_ctx(mod, overrides: dict):
    mod.LW_BOUNDARY_OVERRIDES.clear()
    mod.LW_BOUNDARY_OVERRIDES.update(overrides)
    return mod.load_framework()


def mixing_time(mod, ctx, dT: float, variant: str) -> float:
    return float(mod.mixing_times_variant(ctx, {"destabilisation": dT}, variant)[0]["mixing_time_sols"])


def rate_of(mod, ctx, dT: float, variant: str) -> float:
    return float(mod.mixing_times_variant(ctx, {"destabilisation": dT}, variant)[0]["destabilisation_K_per_sol"])


def dT_req(mod, ctx, n: float, variant: str, dtheta: float | None = None) -> float:
    target = (mod.DELTA_THETA_4KM if dtheta is None else dtheta) / n
    f = lambda d: rate_of(mod, ctx, d, variant) - target  # noqa: E731
    lo, hi = -5.0, 20.0
    if f(lo) * f(hi) > 0:
        return float("nan")
    root = brentq(f, lo, hi, xtol=1e-6)
    if dtheta is None:
        mt = mixing_time(mod, ctx, root, variant)
        if not abs(mt - n) <= 1e-5 * n:
            raise SystemExit(f"mixing_time(dT_req) = {mt} != {n}")
    return float(root)


def main() -> int:
    C.R2_OUT.mkdir(parents=True, exist_ok=True)
    mod = C.load_figure_module()
    sol = float(mod.SOL_LENGTH)
    dtheta = float(mod.DELTA_THETA_4KM)
    contrasts = literature_contrasts()

    # 1. Matriz -------------------------------------------------------------
    matrix_rows, dtreq_rows, results = [], [], {}
    for b in BOUNDARIES:
        ctx = make_ctx(mod, {"radiation.surface_lw_boundary": b})
        for v in mod.SW_VARIANTS:
            thr = mod.compute_thresholds_variant(ctx, v)
            row = {"lw_boundary": b, "sw_variant": v,
                   "heating_threshold_layer0_K": thr["heating"],
                   "destabilisation_threshold_0_4km_K": thr["destabilisation"]}
            for label, c, _ in contrasts:
                row[f"mix_{label}_sols"] = mixing_time(mod, ctx, c, v)
            req = {n: dT_req(mod, ctx, n, v) for n in NS}
            row.update({f"dT_req_{n}sol_K": req[n] for n in NS})
            matrix_rows.append(row)
            results[(b, v)] = {"thr": thr, "req": req, "row": row}
            for n in NS:
                for label, c, src in contrasts:
                    dtreq_rows.append({"lw_boundary": b, "sw_variant": v, "N_sols": n, "dT_req_K": req[n],
                                       "contrast_label": label, "contrast_K": c,
                                       "ratio_contrast_over_dT_req": c / req[n],
                                       "contrast_reaches_dT_req": bool(c >= req[n]), "contrast_source": src})
    with (C.R2_OUT / "t5_matrix.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(matrix_rows[0]))
        writer.writeheader()
        writer.writerows(matrix_rows)
    with (C.R2_OUT / "t5_dTreq.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(dtreq_rows[0]))
        writer.writeheader()
        writer.writerows(dtreq_rows)

    # 2. D12 ----------------------------------------------------------------
    k4 = {v: results[("kirchhoff", v)]["req"][4] for v in mod.SW_VARIANTS}
    headline = min(k4, key=k4.get)

    # 3. Ciclo diurno real ----------------------------------------------------
    real_rows, real = [], {}
    arch = {}
    for run in ("TMN91", "0.1flhc"):
        ds = netCDF4.Dataset(C.ARCHIVE_RUNS[run])
        tsk = np.asarray(ds["TSK"][:, 0, 0], dtype=np.float64)
        ta = np.asarray(ds["T_PHY"][:, 0, 0, 0], dtype=np.float64)
        z = np.asarray((ds["PH"][:, :2, 0, 0] + ds["PHB"][:, :2, 0, 0]) / mod.GRAVITY, dtype=np.float64)
        times = [str(netCDF4.chartostring(ds["Times"][i])) for i in range(len(tsk))]
        ds.close()
        arch[run] = {"dT": tsk - ta, "times": times, "z_mass0": 0.5 * (z[:, 0] + z[:, 1]) - z[:, 0]}
    for b in BOUNDARIES:
        ctx = make_ctx(mod, {"radiation.surface_lw_boundary": b})
        start = float(ctx["config"].site.start_local_time)
        z0_tp = float(ctx["altitude"][0])
        for run in ("TMN91", "0.1flhc"):
            a = arch[run]
            per_sol = {}
            all_rates = []
            for s, (i, j) in SOL_INDICES.items():
                rates = []
                for r in range(i, j):
                    hh, mm = a["times"][r][-8:-3].split(":")
                    lt = int(hh) + int(mm) / 60.0
                    ctx["state"].time = ((lt - start) % 24.0) / 24.0 * sol
                    rates.append(mod.destabilisation_rate(ctx, float(a["dT"][r])))
                rates = np.array(rates)
                all_rates.append(rates)
                mean_rate = float(rates.mean())
                per_sol[s] = mean_rate
                real_rows.append({"lw_boundary": b, "run": run, "sol": s,
                                  "mean_contrast_K": float(a["dT"][i:j].mean()),
                                  "max_contrast_K": float(a["dT"][i:j].max()),
                                  "sol_mean_destabilisation_K_per_sol": mean_rate,
                                  "max_instant_destabilisation_K_per_sol": float(rates.max()),
                                  "mixing_time_sols": dtheta / mean_rate if mean_rate > 0 else "nunca"})
            allr = np.concatenate(all_rates)
            mr = float(allr.mean())
            real_rows.append({"lw_boundary": b, "run": run, "sol": "2-5",
                              "mean_contrast_K": float(a["dT"][87:663].mean()),
                              "max_contrast_K": float(a["dT"][87:663].max()),
                              "sol_mean_destabilisation_K_per_sol": mr,
                              "max_instant_destabilisation_K_per_sol": float(allr.max()),
                              "mixing_time_sols": dtheta / mr if mr > 0 else "nunca"})
            real[(b, run)] = {"mean_rate_2_5": mr, "per_sol": per_sol}
        ctx["state"].time = (12.0 - start) / 24.0 * sol
    with (C.R2_OUT / "t5_real_diurnal.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(real_rows[0]))
        writer.writeheader()
        writer.writerows(real_rows)
    height_note = {"archive_Ta_height_m_mean": float(np.mean(arch["TMN91"]["z_mass0"][87:663])),
                   "titan_pbl_layer0_height_m": z0_tp}

    # 4. Sensibilidad OAT y esquina -----------------------------------------
    base_cfg = make_ctx(mod, {"radiation.surface_lw_boundary": "kirchhoff"})["config"]
    s0_preset = float(base_cfg.radiation.solar_constant)

    def overrides_for(values: dict) -> dict:
        ov = {"radiation.surface_lw_boundary": "kirchhoff"}
        for k, val in values.items():
            if k == "S0_factor":
                ov["radiation.solar_constant"] = s0_preset * val
            else:
                ov[k] = val
        return ov

    def evaluate(values: dict, variant: str = headline) -> dict:
        ctx = make_ctx(mod, overrides_for(values))
        thr = mod.compute_thresholds_variant(ctx, variant)
        return {"heating_threshold_layer0_K": thr["heating"],
                "destabilisation_threshold_0_4km_K": thr["destabilisation"],
                "mix_3K_sols": mixing_time(mod, ctx, 3.0, variant),
                "dT_req_4sol_K": dT_req(mod, ctx, 4, variant)}

    oat_rows, best, ties = [], {}, {}
    for factor, values in OAT.items():
        scores = {}
        for val in values:
            vals = dict(PRESET)
            vals[factor] = val
            res = evaluate(vals)
            scores[val] = res["dT_req_4sol_K"]
            oat_rows.append({"factor": factor, "value": val, "is_preset": val == PRESET[factor],
                             "sw_variant": headline, "role": "barrido declarado", **res})
        finite = {k: s for k, s in scores.items() if np.isfinite(s)}
        best[factor] = min(finite, key=finite.get)  # empate: el primero de la lista
        ties[factor] = bool(np.ptp(list(finite.values())) == 0.0)
    # Complemento descriptivo anadido tras ver que, con la variante principal
    # "lw", alpha y S0 no actuan: el mismo barrido con la variante "noon".
    if headline == "lw":
        for factor in ("radiation.surface_albedo", "S0_factor"):
            for val in OAT[factor]:
                vals = dict(PRESET)
                vals[factor] = val
                oat_rows.append({"factor": factor, "value": val, "is_preset": val == PRESET[factor],
                                 "sw_variant": "noon", "role": "complemento descriptivo (no entra en la esquina)",
                                 **evaluate(vals, "noon")})
    with (C.R2_OUT / "t5_sensitivity_oat.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(oat_rows[0]))
        writer.writeheader()
        writer.writerows(oat_rows)
    corner_res = evaluate(best)
    corner_res["dT_req_5sol_K"] = dT_req(mod, make_ctx(mod, overrides_for(best)), 5, headline)
    preset_res = evaluate(dict(PRESET))
    corner = {"task": "R2-T5.4", "lw_boundary": "kirchhoff", "sw_variant": headline,
              "rule": "para cada factor, el valor que minimiza dT_req(4) en el barrido OAT; se evalua la combinacion",
              "chosen_values": best, "factor_without_effect_tie": ties,
              "tie_rule": "si todos los valores del factor dan el mismo dT_req(4), se toma el primero de la lista (sin efecto)",
              "corner_results": corner_res, "preset_results": preset_res,
              "limitation": ("alpha solo actua por la SW reflejada hacia la atmosfera; su efecto principal sobre Ts "
                             "via el balance superficial queda fuera, porque los contrastes se toman de corridas publicadas"),
              "max_literature_contrasts_K": {"v0_rango_2K": 2.0,
                                             "archivo_max_D4": float(C.number_r1("archive_contrast_max_sols2to5_suppressed"))}}
    C.dump_json(C.R2_OUT / "t5_corner.json", corner)

    # 5. Incertidumbre HASI ---------------------------------------------------
    ctx = make_ctx(mod, {"radiation.surface_lw_boundary": "kirchhoff"})
    state = ctx["state"]
    theta = state.potential_temperature()
    top = int(np.argmin(np.abs(ctx["altitude"] - 4000.0)))
    s0 = SIGMA_T_HASI * float(theta[0] / state.atmospheric_temperature[0])
    s1 = SIGMA_T_HASI * float(theta[top] / state.atmospheric_temperature[top])
    sig = float(np.sqrt(s0**2 + s1**2))
    unc = {
        "task": "R2-T5.5",
        "delta_theta_0_4km_K": float(theta[top] - theta[0]), "delta_theta_used_K": dtheta,
        "hasi_l4_file": "titan_pbl/data/huygens/HASI_L4_ATMO_PROFILE_DESCEN_TAB_v2.txt (identico a data/HASI_L4_ATMO_PROFILE_DESCEN.TAB_v2.txt)",
        "hasi_l4_columns": "Time, Altitude, Pressure, Temperature, Density (etiqueta PDS: DESCRIPTION = NULL en las 5); sin columna de incertidumbre",
        "label_source": "PDS3 label of the product HASI_L4_ATMO_PROFILE_DESCEN (doi:10.17189/bcm3-4061)",
        "fulchignoni_2005_quote": "Fig. 5 caption (Nature 438, p. 787): 'Temperature uncertainty is ^0.25 K in the range from 60 to 110 K, and^1 K above 110 K.' (capa de texto; '^' = +-)",
        "sigma_T_K": SIGMA_T_HASI,
        "theta_over_T": {"layer0": float(theta[0] / state.atmospheric_temperature[0]),
                         "layer_4km": float(theta[top] / state.atmospheric_temperature[top]), "layer_4km_index": top},
        "independent": {"sigma_delta_theta_K": sig, "relative": sig / dtheta,
                        "mixing_time_relative_uncertainty": sig / dtheta,
                        "note": "el tiempo de mezcla es proporcional a Delta_theta"},
        "systematic": {"sigma_delta_theta_K": 0.0, "note": "un error comun a todos los puntos se cancela en la diferencia"},
    }
    req_unc = {}
    for b in BOUNDARIES:
        ctxb = make_ctx(mod, {"radiation.surface_lw_boundary": b})
        req_unc[b] = {f"dtheta={dtheta + s:.3f}": dT_req(mod, ctxb, 4, headline, dtheta + s) for s in (-sig, 0.0, sig)}
    unc["dT_req_4sol_K_with_dtheta_pm_sigma_headline_variant"] = req_unc
    C.dump_json(C.R2_OUT / "t5_hasi_uncertainty.json", unc)

    # numbers_r2.json ---------------------------------------------------------
    fn = "code/analysis/r2_t5_sw_dtreq.py::main"
    nums = [C.entry("headline_sw_variant", headline, "texto", "kirchhoff", headline, C.rel(C.R2_OUT / "t5_dTreq.csv"),
                    f"{fn} (D12: argmin dT_req(4) con kirchhoff)", None, "regla D12")]
    sig_thr = (2, "como los umbrales de R1: eps con 2 cifras")
    for (b, v), r in results.items():
        for n in NS:
            nums.append(C.entry(f"dTreq_{n}sol_{b}_{v}", r["req"][n], "K", b, v, C.rel(C.R2_OUT / "t5_dTreq.csv"), fn, *sig_thr))
        for label, c, _ in contrasts:
            if label.startswith("v0_rango") or label.startswith("Hartwick") or label.startswith("archivo_max"):
                val = r["row"][f"mix_{label}_sols"]
                nums.append(C.entry(f"swmatrix_mix_{label}_{b}_{v}", val, "Titan sols (inf = nunca)", b, v,
                                    C.rel(C.R2_OUT / "t5_matrix.csv"), fn, 2, "Delta_theta (3 cifras) / tasa de 2 cifras"))
    for (b, run), rr in real.items():
        tag = run.replace(".", "p")
        nums.append(C.entry(f"realdiurnal_mean_rate_{tag}_{b}", rr["mean_rate_2_5"], "K per Titan sol", b, "ciclo diurno real",
                            C.rel(C.R2_OUT / "t5_real_diurnal.csv"), fn, 2, "media de 576 registros"))
    nums.append(C.entry("corner_dTreq_4sol", corner_res["dT_req_4sol_K"], "K", "kirchhoff", headline,
                        C.rel(C.R2_OUT / "t5_corner.json"), fn, 2, "esquina declarada"))
    nums.append(C.entry("corner_dTreq_5sol", corner_res["dT_req_5sol_K"], "K", "kirchhoff", headline,
                        C.rel(C.R2_OUT / "t5_corner.json"), fn, 2, "esquina declarada"))
    nums.append(C.entry("hasi_unc_sigma_dtheta_independent", sig, "K", "n/a", "n/a",
                        C.rel(C.R2_OUT / "t5_hasi_uncertainty.json"), fn, 2, "propagacion de +-0.25 K"))
    C.update_numbers_r2(nums, ("headline_sw_variant", "dTreq_", "swmatrix_", "realdiurnal_", "corner_", "hasi_unc_"), fn)

    out = {"headline_sw_variant": headline, "dT_req_kirchhoff_4": k4,
           "matrix": [{k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()} for r in matrix_rows],
           "real_diurnal": {f"{b}|{run}": rr for (b, run), rr in real.items()}, "height_note": height_note,
           "corner": corner, "hasi": {k: unc[k] for k in ("independent", "dT_req_4sol_K_with_dtheta_pm_sigma_headline_variant")}}
    C.dump_json(C.R2_OUT / "t5_summary.json", out)
    print(json.dumps(C.clean_nan(out), indent=1, ensure_ascii=False, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
