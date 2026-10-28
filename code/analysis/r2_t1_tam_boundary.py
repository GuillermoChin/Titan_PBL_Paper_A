#!/usr/bin/env python
"""
r2_t1_tam_boundary.py -- stage r2, analysis T1.

Classifies the lower longwave boundary of TAM (Han & Lora 2025) from the 8
diurnal files of Zenodo 16791001 (cases 1-4, Huygens and Dragonfly).

---------------------------------------------------------------------------
CRITERIOS (declarados antes de ejecutar)
---------------------------------------------------------------------------
Para cada registro:
  eps_O = lwup / (sigma * tsurf^4)
  eps_K = (lwup - lwdn) / (sigma * tsurf^4 - lwdn), solo en registros con
          |sigma * tsurf^4 - lwdn| >= 0.05 W m-2 (filtro declarado; se
          reporta el numero de registros que quedan).
sigma en tres valores, sin asumir ninguno: 5.670374419e-8, 5.67e-8 y
5.6734e-8 (constante de FMS).

CV = coeficiente de variacion por archivo (std / |media|).
  "eps = 1, sin reflexion": |mediana(eps_O) - 1| <= 1e-3 y CV(eps_O) <= 1e-3.
  "O": CV(eps_O) <= 1e-3, |mediana(eps_O) - 1| > 1e-3 y CV(eps_K) >= 10 CV(eps_O).
  "K": CV(eps_K) <= 1e-3 y CV(eps_O) >= 10 CV(eps_K).
  En cualquier otro caso: "No determinable".
(Las reglas se evaluan en ese orden; si CV(eps_K) no existe porque el filtro
 deja menos de 2 registros, "O" y "K" no pueden cumplirse.)

PRUEBA AUXILIAR (existen sw_sfc, sens y hf): residuo del balance superficial
  R_full = sw_sfc + lwdn       - lwup - sens - hf
  R_eps  = sw_sfc + eps * lwdn - lwup - sens - hf
con eps = mediana(eps_O) del archivo (sigma de FMS) y, como sensibilidad,
eps = 0.95. Convenciones de signo segun long_name: sw_sfc "SW net flux down
at surface" (+ hacia abajo, neto), lwup/lwdn magnitudes, hf "Net flux into
the ground" (+ hacia el suelo), sens "Sensible heat flux" (+ hacia arriba,
convencion de FMS; se verifica con el signo del residuo). Una variante se
prefiere si su RMS es <= 0.5 veces el de la otra y |RMS_a - RMS_b| > 0.05
W m-2; si no, "sin preferencia". Si eps = 1, las variantes son identicas y
la prueba "no discrimina".

ADEMAS (para Ts - Ta, T4): altura del nivel mas bajo (zg) y su temperatura
(temp); duracion del dia de Titan en TAM por la periodicidad de sw_sfc:
  metodo A: minimos cuadrados de sw_sfc contra media + 3 armonicos de periodo
            P, barrido fino de P en [12, 20] dias FMS; P* = argmin RSS;
  metodo B: regresion lineal de los instantes de "amanecer" (cruce ascendente
            de sw_sfc por el 5 % de su maximo, interpolado linealmente) contra
            el numero de dia; P = pendiente.
Unidad de tiempo: dia FMS = 86400 s ("days since 0000-00-00").

AGREGACION Y DESCRIPTIVOS (anadidos DESPUES de ver la primera corrida; no
alteran la clasificacion por archivo):
- El criterio no declara como agregar las 8 clases. Regla adoptada: si las 8
  coinciden, esa es la clase de TAM; si no, "No determinable", con el conteo
  por clase.
- Prueba de balance con emision explicita (anadida tras ver que Han & Lora
  2025, p. 3, declaran "surface thermal emissivity of 0.95" y que su IR neto
  en superficie, p. 5, es ~0.07 W m-2, incompatible con lwup_sfc - lwdn_sfc):
    R_O1(0.95) = sw_sfc + lwdn - 0.95 sigma Ts^4 - sens - hf
    R_K(0.95)  = sw_sfc + 0.95 lwdn - 0.95 sigma Ts^4 - sens - hf
  y el eps que mejor cierra sw + lwdn - eps sigma Ts^4 - sens - hf = 0. Es
  descriptiva: no cambia la clasificacion por criterio, que se basa en
  lwup_sfc.
- Descriptivos por archivo (sigma de FMS): ajuste por minimos cuadrados de
  lwup = a sigma Ts^4 y de lwup = a sigma Ts^4 + b lwdn, con su RMS; CV(eps_O)
  desplazando tsurf +-1 y +-2 registros; medias de cada termino del balance.

Nota: todas las variables son medias temporales sobre el intervalo de salida
(cell_methods "time: mean"); sigma*mean(Ts)^4 difiere de mean(sigma Ts^4) en
~6 var(Ts)/Ts^2, del orden de 1e-6, despreciable frente a 1e-3.

Salidas: build/r2/t1_tam_boundary.csv y .json
"""

from __future__ import annotations

import csv
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
from scipy.optimize import minimize_scalar  # noqa: E402

import r2_common as C  # noqa: E402

FILTER_K = 0.05
TOL = 1.0e-3
VARS = ("lwup_sfc", "lwdn_sfc", "tsurf", "sw_sfc", "sens", "hf", "zg", "temp", "ps", "time", "pfull", "phalf")
R_GAS = 296.8   # solo para la estimacion hipsometrica de la altura del nivel mas bajo
GRAVITY = 1.352


def cv(x: np.ndarray) -> float:
    if x.size < 2:
        return float("nan")
    return float(np.std(x, ddof=1) / abs(np.mean(x)))


def classify(med_o, cv_o, cv_k) -> str:
    if abs(med_o - 1.0) <= TOL and cv_o <= TOL:
        return "eps=1, sin reflexion"
    k_ok = np.isfinite(cv_k)
    if cv_o <= TOL and abs(med_o - 1.0) > TOL and k_ok and cv_k >= 10.0 * cv_o:
        return "O"
    if k_ok and cv_k <= TOL and cv_o >= 10.0 * cv_k:
        return "K"
    return "No determinable"


def harmonic_rss(t, y, period, n_harm=3):
    cols = [np.ones_like(t)]
    for h in range(1, n_harm + 1):
        w = 2.0 * np.pi * h * t / period
        cols += [np.cos(w), np.sin(w)]
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ coef
    return float(r @ r)


def day_length(t, sw):
    grid = np.arange(12.0, 20.0, 0.001)
    rss = np.array([harmonic_rss(t, sw, p) for p in grid])
    p0 = grid[int(np.argmin(rss))]
    res = minimize_scalar(lambda p: harmonic_rss(t, sw, p), bounds=(p0 - 0.002, p0 + 0.002), method="bounded",
                          options={"xatol": 1e-7})
    period_a = float(res.x)
    # Metodo B: cruces ascendentes por el 5 % del maximo.
    thr = 0.05 * float(np.max(sw))
    up = np.where((sw[:-1] < thr) & (sw[1:] >= thr))[0]
    tc = t[up] + (thr - sw[up]) / (sw[up + 1] - sw[up]) * (t[up + 1] - t[up])
    n = np.arange(tc.size)
    slope, intercept = np.polyfit(n, tc, 1)
    resid = tc - (slope * n + intercept)
    se = float(np.sqrt(np.sum(resid**2) / (n.size - 2) / np.sum((n - n.mean()) ** 2)))
    return {
        "method_A_period_days": period_a,
        "method_A_grid_step_days": 0.001,
        "method_B_period_days": float(slope),
        "method_B_stderr_days": se,
        "method_B_n_crossings": int(tc.size),
        "method_B_threshold_W_m2": thr,
        "record_span_days": float(t[-1] - t[0] + (t[1] - t[0])),
    }


def main() -> int:
    C.R2_OUT.mkdir(parents=True, exist_ok=True)
    md5 = C.check_tam_md5()  # P3
    rows, files = [], {}
    for case in C.TAM_CASES:
        path = C.tam_path(case)
        ds = netCDF4.Dataset(path)
        meta = {v: {"dims": list(ds[v].dimensions),
                    "units": getattr(ds[v], "units", None),
                    "long_name": getattr(ds[v], "long_name", None),
                    "cell_methods": getattr(ds[v], "cell_methods", None)} for v in VARS}
        d = {v: np.asarray(ds[v][:], dtype=np.float64) for v in VARS}
        ds.close()
        up, dn, ts = d["lwup_sfc"], d["lwdn_sfc"], d["tsurf"]
        info = {"file": str(path), "md5": md5[path.name]["md5"], "variables": meta,
                "n_records": int(ts.size), "per_sigma": {}}
        for sname, sigma in C.SIGMAS.items():
            bb = sigma * ts**4
            eps_o = up / bb
            den = bb - dn
            mask = np.abs(den) >= FILTER_K
            eps_k = (up[mask] - dn[mask]) / den[mask]
            med_o, cv_o = float(np.median(eps_o)), cv(eps_o)
            med_k = float(np.median(eps_k)) if eps_k.size else float("nan")
            cv_k = cv(eps_k)
            cls = classify(med_o, cv_o, cv_k)
            res = {
                "sigma": sigma,
                "median_eps_O": med_o, "mean_eps_O": float(np.mean(eps_o)), "cv_eps_O": cv_o,
                "min_eps_O": float(eps_o.min()), "max_eps_O": float(eps_o.max()),
                "n_K_after_filter": int(mask.sum()),
                "median_eps_K": med_k, "cv_eps_K": cv_k,
                "abs_den_K_range_W_m2": [float(np.abs(den).min()), float(np.abs(den).max())],
                "class": cls,
            }
            info["per_sigma"][sname] = res
            rows.append({"case": case, "sigma_name": sname, "sigma": sigma, **{k: v for k, v in res.items()
                         if k not in ("sigma", "abs_den_K_range_W_m2")},
                         "n_records": int(ts.size)})
        classes = {r["class"] for r in info["per_sigma"].values()}
        info["class_all_sigmas"] = classes.pop() if len(classes) == 1 else "depende de sigma"

        # Prueba auxiliar del balance superficial.
        eps_med = info["per_sigma"]["FMS_5.6734e-8"]["median_eps_O"]
        aux = {}
        for label, eps in (("eps_median_O_FMS", eps_med), ("eps_0.95", 0.95)):
            r_full = d["sw_sfc"] + dn - up - d["sens"] - d["hf"]
            r_eps = d["sw_sfc"] + eps * dn - up - d["sens"] - d["hf"]
            rms_f, rms_e = float(np.sqrt(np.mean(r_full**2))), float(np.sqrt(np.mean(r_eps**2)))
            if abs(eps - 1.0) < 1e-12:
                pref = "no discrimina (eps = 1)"
            elif rms_f <= 0.5 * rms_e and abs(rms_f - rms_e) > 0.05:
                pref = "absorcion de lwdn completo"
            elif rms_e <= 0.5 * rms_f and abs(rms_f - rms_e) > 0.05:
                pref = "absorcion de eps*lwdn"
            else:
                pref = "sin preferencia"
            aux[label] = {"eps": eps, "rms_R_full_W_m2": rms_f, "rms_R_eps_W_m2": rms_e,
                          "mean_R_full_W_m2": float(np.mean(r_full)), "mean_R_eps_W_m2": float(np.mean(r_eps)),
                          "preferred": pref}
            # signo alternativo de sens (por si fuera + hacia abajo), solo descriptivo
            r_alt = d["sw_sfc"] + dn - up + d["sens"] - d["hf"]
            aux[label]["rms_R_full_sens_sign_flipped_W_m2"] = float(np.sqrt(np.mean(r_alt**2)))
        info["auxiliary_surface_balance"] = aux

        # Descriptivos (no son criterio).
        s_fms = C.SIGMAS["FMS_5.6734e-8"]
        bb = s_fms * ts**4
        a1, *_ = np.linalg.lstsq(bb[:, None], up, rcond=None)
        r1 = up - bb * a1[0]
        A2 = np.column_stack([bb, dn])
        a2, *_ = np.linalg.lstsq(A2, up, rcond=None)
        r2 = up - A2 @ a2
        lag_cv = {}
        for lag in (-2, -1, 0, 1, 2):
            e = up[3:-3] / (s_fms * np.roll(ts, lag)[3:-3] ** 4)
            lag_cv[str(lag)] = cv(e)
        info["descriptive_FMS_sigma"] = {
            "fit_lwup_eq_a_bb": {"a": float(a1[0]), "rms_W_m2": float(np.sqrt(np.mean(r1**2)))},
            "fit_lwup_eq_a_bb_plus_b_lwdn": {"a": float(a2[0]), "b": float(a2[1]),
                                             "rms_W_m2": float(np.sqrt(np.mean(r2**2)))},
            "cv_eps_O_with_tsurf_shifted_by_records": lag_cv,
            "balance_term_means_W_m2": {k: float(np.mean(d[k])) for k in ("sw_sfc", "sens", "hf")}
            | {"lwdn_sfc": float(np.mean(dn)), "lwup_sfc": float(np.mean(up)), "sigma_Ts4": float(np.mean(bb))},
            "R_full_std_W_m2": float(np.std(d["sw_sfc"] + dn - up - d["sens"] - d["hf"])),
        }
        bal = {}
        rest = d["sw_sfc"] - d["sens"] - d["hf"]
        for sname, sigma in C.SIGMAS.items():
            bbs = sigma * ts**4
            r_o = rest + dn - 0.95 * bbs
            r_k = rest + 0.95 * dn - 0.95 * bbs
            eps_fit = float(np.sum((rest + dn) * bbs) / np.sum(bbs * bbs))
            r_fit = rest + dn - eps_fit * bbs
            bal[sname] = {
                "R_O1_eps095_mean_W_m2": float(np.mean(r_o)), "R_O1_eps095_rms_W_m2": float(np.sqrt(np.mean(r_o**2))),
                "R_K_eps095_mean_W_m2": float(np.mean(r_k)), "R_K_eps095_rms_W_m2": float(np.sqrt(np.mean(r_k**2))),
                "eps_best_fit_full_lwdn_absorption": eps_fit, "rms_at_best_fit_W_m2": float(np.sqrt(np.mean(r_fit**2))),
            }
        info["balance_test_explicit_emission"] = bal
        info["net_IR_up_mean_W_m2"] = {"lwup_minus_lwdn": float(np.mean(up - dn)),
                                       "0.95_sigmaFMS_Ts4_minus_lwdn": float(np.mean(0.95 * s_fms * ts**4 - dn))}

        # Nivel mas bajo y contraste Ts - Ta.
        z_low, t_low = d["zg"][:, -1], d["temp"][:, -1]
        p_low_sigma = d["pfull"][-1] / d["phalf"][-1]  # sigma del nivel mas bajo (supuesto)
        dz_hyps = R_GAS * t_low / GRAVITY * np.log(1.0 / p_low_sigma)
        contrast = ts - t_low
        info["lowest_level"] = {
            "index": int(d["temp"].shape[1] - 1),
            "pfull_hPa": float(d["pfull"][-1]),
            "phalf_bottom_hPa": float(d["phalf"][-1]),
            "zg_mean_m": float(np.mean(z_low)), "zg_min_m": float(z_low.min()), "zg_max_m": float(z_low.max()),
            "temp_mean_K": float(np.mean(t_low)),
            "hypsometric_height_above_surface_m_mean": float(np.mean(dz_hyps)),
            "hypsometric_note": "R = 296.8, g = 1.352, p_k/p_s = pfull/phalf_bottom (supuesto sigma); solo estimacion",
            "second_level_zg_mean_m": float(np.mean(d["zg"][:, -2])),
        }
        info["Ts_minus_Ta"] = {"max_K": float(contrast.max()), "mean_K": float(contrast.mean()),
                               "min_K": float(contrast.min())}
        info["Ts_diurnal"] = {"tsurf_mean_K": float(ts.mean())}
        info["day_length"] = day_length(d["time"], d["sw_sfc"])
        files[case] = info
        print(case, {s: info["per_sigma"][s]["class"] for s in C.SIGMAS}, "| aux:",
              aux["eps_median_O_FMS"]["preferred"], "| day A/B:",
              round(info["day_length"]["method_A_period_days"], 4), round(info["day_length"]["method_B_period_days"], 4))

    with (C.R2_OUT / "t1_tam_boundary.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    per_file = [c["class_all_sigmas"] for c in files.values()]
    counts = {k: per_file.count(k) for k in sorted(set(per_file))}
    overall = set(per_file)
    periods_a = np.array([f["day_length"]["method_A_period_days"] for f in files.values()])
    periods_b = np.array([f["day_length"]["method_B_period_days"] for f in files.values()])
    report = {
        "task": "R2-T1",
        "guard": C.guard_info(),
        "criteria": {"tol": TOL, "filter_K_W_m2": FILTER_K, "sigmas": C.SIGMAS},
        "md5_check_P3": md5,
        "classification_TAM": overall.pop() if len(overall) == 1 else "No determinable",
        "classification_counts_by_file": counts,
        "aggregation_rule": "unanime -> esa clase; si no, No determinable (regla anadida tras ver resultados)",
        "balance_test_summary_FMS": {
            "eps_best_fit_range": [min(f["balance_test_explicit_emission"]["FMS_5.6734e-8"]["eps_best_fit_full_lwdn_absorption"] for f in files.values()),
                                   max(f["balance_test_explicit_emission"]["FMS_5.6734e-8"]["eps_best_fit_full_lwdn_absorption"] for f in files.values())],
            "R_O1_eps095_mean_range_W_m2": [min(f["balance_test_explicit_emission"]["FMS_5.6734e-8"]["R_O1_eps095_mean_W_m2"] for f in files.values()),
                                            max(f["balance_test_explicit_emission"]["FMS_5.6734e-8"]["R_O1_eps095_mean_W_m2"] for f in files.values())],
            "R_K_eps095_mean_range_W_m2": [min(f["balance_test_explicit_emission"]["FMS_5.6734e-8"]["R_K_eps095_mean_W_m2"] for f in files.values()),
                                           max(f["balance_test_explicit_emission"]["FMS_5.6734e-8"]["R_K_eps095_mean_W_m2"] for f in files.values())],
            "note": "descriptivo, posterior: el balance superficial cierra con emision 0.95 sigma Ts^4 y absorcion de lwdn completo (O1)",
        },
        "day_length_summary": {
            "method_A_mean_days": float(periods_a.mean()), "method_A_range_days": [float(periods_a.min()), float(periods_a.max())],
            "method_B_mean_days": float(periods_b.mean()), "method_B_range_days": [float(periods_b.min()), float(periods_b.max())],
            "method_A_mean_s": float(periods_a.mean() * C.TAM_DAY_S),
            "day_unit": "dia FMS = 86400 s",
        },
        "files": files,
    }
    C.dump_json(C.R2_OUT / "t1_tam_boundary.json", report)
    print("CLASIFICACION TAM:", report["classification_TAM"])
    print("DIA DE TITAN EN TAM:", report["day_length_summary"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
