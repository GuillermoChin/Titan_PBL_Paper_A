#!/usr/bin/env python
"""
r2_t2_energy_partition.py -- stage r2, analysis T2.

Surface energy partition per sol (mtWRF archive) and per Titan day (TAM).

---------------------------------------------------------------------------
7.1 ARCHIVO mtWRF (declarado antes de ejecutar)
---------------------------------------------------------------------------
Corridas: baseline TMN91 (principal, H-23), baseline alternativa
(baseline_5sol) y 0.1flhc. Integrales por sol, en MJ m-2, Tsols 2-5, con los
mismos limites de sol que R1/T5 (SOL_INDICES del script de figuras:
{2: (87, 231), 3: (231, 375), 4: (375, 519), 5: (519, 663)}, 144 registros
por Tsol, del registro de 00:00 al de 23:50).
Regla de integracion (declarada): rectangulos sobre los 144 registros del
sol, dt = Tsol / 144 = 9570 s (XTIME avanza 10 "minutos" del Tsol por
registro; 1440 minutos = 1 Tsol = 1 378 080 s).
Variables:
  GSW; GLW (D); U = FLUX_LW_U[:,0]; eps sigma Ts^4 con eps = 0.95 y
  sigma = 5.67051e-8 (STBOLT de WRF, el mejor ajuste de R0/T3); XINET;
  R_O1 = (1 - eps) GLW (absorcion extra de la superficie respecto de
  Kirchhoff); HFX y LH; G (flujo al suelo).
G: se pide SOILCOND con el gradiente de TSLB si hay profundidades; si
no, residuo del balance marcado como tal. Las baselines no traen ZS ni una
conductividad; SOILCOND se verifica antes contra el balance. Si no lo cierra,
G = GSW + XINET - HFX - LH (RESIDUO), y se reporta ademas su diferencia con
FLUX[:,0] ("CONDUCTIVE HEAT FLUX", interfaz superior del suelo).
Convencion de signo: se cita la description de cada variable en el JSON.
Se reporta el cociente R_O1/HFX por sol.

---------------------------------------------------------------------------
7.2 TAM, casos Huygens 1-4 (declarado antes de ejecutar)
---------------------------------------------------------------------------
Integrales por dia de Titan para los 30 dias centrales (de 50); como
sensibilidad, los 50. Variables: sw_sfc, lwdn, lwup, sens, hf.
Termino equivalente a R_O1: (1 - eps) lwdn si T1 clasifico la frontera como
"O" (eps = mediana(eps_O) con sigma de FMS); 0 si "eps = 1" o "K". Si T1 dio
"No determinable" el termino es "No determinable"; el valor que tendria bajo
la hipotesis O se da en una columna aparte, marcado como condicional.
Anadido tras T1 (descriptivo): el balance superficial de TAM cierra con
emision 0.95 sigma Ts^4 y absorcion de lwdn completo; se da tambien el valor
condicional (1 - 0.95) lwdn, que corresponde a esa lectura del balance.
Ventanas (declaradas): dia de Titan de T1 (P = 15.95 dias FMS); cada registro
es una media temporal y se supone centrado en su intervalo de 0.6667 dias;
dia k = [T0 + k P, T0 + (k+1) P) con T0 = time[0] - dt/2; la integral pondera
cada registro por su solape con la ventana. Dias centrales: k = 10..39
(se descartan 10 al principio y 10 al final).

Salidas: build/r2/t2_energy_partition.csv (formato largo), t2_energy_partition.json
         y claves "energy_*" en build/r2/numbers_r2.json.
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

import r2_common as C  # noqa: E402

EPS = 0.95
STBOLT = 5.67051e-8
SOL_INDICES = {2: (87, 231), 3: (231, 375), 4: (375, 519), 5: (519, 663)}
MJ = 1.0e6
ARCHIVE_VARS = ("GSW", "GLW", "XINET", "HFX", "LH", "SOILCOND", "TSK")
TAM_VARS = ("sw_sfc", "lwdn_sfc", "lwup_sfc", "sens", "hf")
TAM_HUYGENS = ("case1_huygens", "case2_huygens", "case3_huygens", "case4_huygens")
N_DAYS, CENTRAL = 50, range(10, 40)


def archive_partition(sol_length: float):
    rows, meta = [], {}
    dt = sol_length / 144.0
    for run, path in C.ARCHIVE_RUNS.items():
        ds = netCDF4.Dataset(path)
        desc = {v: {"description": getattr(ds[v], "description", ""), "units": getattr(ds[v], "units", "")}
                for v in ARCHIVE_VARS + ("FLUX_LW_U", "FLUX", "TSLB")}
        xt = np.asarray(ds["XTIME"][:], dtype=np.float64)
        d = {v: np.asarray(ds[v][:, 0, 0], dtype=np.float64) for v in ARCHIVE_VARS}
        d["U"] = np.asarray(ds["FLUX_LW_U"][:, 0, 0, 0], dtype=np.float64)
        d["FLUX0"] = np.asarray(ds["FLUX"][:, 0, 0, 0], dtype=np.float64)
        has_zs = "ZS" in ds.variables
        ds.close()
        d["EMIS"] = EPS * STBOLT * d["TSK"] ** 4
        d["R_O1"] = (1.0 - EPS) * d["GLW"]
        d["G_residual"] = d["GSW"] + d["XINET"] - d["HFX"] - d["LH"]
        sl = slice(SOL_INDICES[2][0], SOL_INDICES[5][1])
        soilcond_res = d["GSW"] + d["XINET"] - d["HFX"] - d["LH"] + d["SOILCOND"]
        meta[run] = {
            "file": C.rel(path),
            "descriptions": desc,
            "xtime_step_minutes": float(np.median(np.diff(xt))),
            "dt_record_s": dt,
            "has_ZS": has_zs,
            "check_SOILCOND_closes_balance_rms_W_m2": float(np.sqrt(np.mean(soilcond_res[sl] ** 2))),
            "check_SOILCOND_closes_balance_with_minus_sign_rms_W_m2": float(np.sqrt(np.mean(
                (d["GSW"] + d["XINET"] - d["HFX"] - d["LH"] - d["SOILCOND"])[sl] ** 2))),
            "check_G_residual_minus_FLUX0_max_abs_W_m2": float(np.max(np.abs(d["G_residual"][sl] - d["FLUX0"][sl]))),
            "check_XINET_minus_(GLW-U)_max_abs_W_m2": float(np.max(np.abs(d["XINET"][sl] - (d["GLW"][sl] - d["U"][sl])))),
            "LH_max_abs_W_m2": float(np.max(np.abs(d["LH"][sl]))),
        }
        for sol, (i, j) in SOL_INDICES.items():
            vals = {}
            for v in ("GSW", "GLW", "U", "EMIS", "XINET", "R_O1", "HFX", "LH", "G_residual", "SOILCOND"):
                vals[v] = float(np.sum(d[v][i:j]) * dt / MJ)
            vals["R_O1_over_HFX"] = vals["R_O1"] / vals["HFX"] if vals["HFX"] != 0 else float("nan")
            for v, val in vals.items():
                rows.append({"source": "mtWRF", "run": run, "period": f"Tsol {sol}", "period_index": sol,
                             "variable": v, "value": val,
                             "units": "-" if v == "R_O1_over_HFX" else "MJ m-2 per Tsol"})
        meta[run]["per_sol"] = {sol: {r["variable"]: r["value"] for r in rows
                                      if r["run"] == run and r["period_index"] == sol} for sol in SOL_INDICES}
    return rows, meta


def tam_windows(time: np.ndarray, period: float):
    dt = float(np.median(np.diff(time)))
    starts = time - dt / 2.0
    ends = time + dt / 2.0
    t0 = starts[0]
    weights = []
    for k in range(N_DAYS):
        a, b = t0 + k * period, t0 + (k + 1) * period
        overlap = np.clip(np.minimum(ends, b) - np.maximum(starts, a), 0.0, None)
        weights.append(overlap)
    return np.array(weights), dt


def tam_partition(t1: dict):
    rows, meta = [], {}
    period = t1["day_length_summary"]["method_A_mean_days"]
    for case in TAM_HUYGENS:
        ds = netCDF4.Dataset(C.tam_path(case))
        d = {v: np.asarray(ds[v][:], dtype=np.float64) for v in TAM_VARS}
        time = np.asarray(ds["time"][:], dtype=np.float64)
        ds.close()
        w, dt = tam_windows(time, period)
        cls = t1["files"][case]["class_all_sigmas"]
        eps = t1["files"][case]["per_sigma"]["FMS_5.6734e-8"]["median_eps_O"]
        r_cond = (1.0 - eps) * d["lwdn_sfc"]
        if cls == "O":
            r_o1, status = r_cond, "O: (1-eps) lwdn, eps = mediana(eps_O), sigma FMS"
        elif cls in ("eps=1, sin reflexion", "K"):
            r_o1, status = np.zeros_like(r_cond), f"{cls}: 0"
        else:
            r_o1, status = None, "No determinable (T1)"
        daily = {}
        for v in TAM_VARS:
            daily[v] = w @ d[v] * C.TAM_DAY_S / MJ
        daily["R_O1_equiv"] = None if r_o1 is None else w @ r_o1 * C.TAM_DAY_S / MJ
        daily["R_O1_equiv_if_O_conditional"] = w @ r_cond * C.TAM_DAY_S / MJ
        daily["R_O1_equiv_if_O_eps095_balance"] = w @ ((1.0 - EPS) * d["lwdn_sfc"]) * C.TAM_DAY_S / MJ
        meta[case] = {"T1_class": cls, "eps_median_O_FMS": eps, "R_O1_status": status,
                      "window_weight_sums_days": [float(x) for x in w.sum(axis=1)], "record_dt_days": dt,
                      "period_days": period}
        for scope, idx in (("30 centrales (k=10..39)", list(CENTRAL)), ("50 dias", list(range(N_DAYS)))):
            summ = {}
            for v, arr in daily.items():
                if arr is None:
                    summ[v] = None
                    continue
                summ[v] = {"mean": float(np.mean(arr[idx])), "min": float(np.min(arr[idx])),
                           "max": float(np.max(arr[idx]))}
                rows.append({"source": "TAM", "run": case, "period": scope, "period_index": "",
                             "variable": v, "value": float(np.mean(arr[idx])),
                             "units": "MJ m-2 per Titan day (media de los dias)"})
            if summ["R_O1_equiv"] is None:
                rows.append({"source": "TAM", "run": case, "period": scope, "period_index": "",
                             "variable": "R_O1_equiv", "value": "No determinable", "units": ""})
            sens = summ["sens"]["mean"]
            ratio = None if summ["R_O1_equiv"] is None else summ["R_O1_equiv"]["mean"] / sens
            rows.append({"source": "TAM", "run": case, "period": scope, "period_index": "",
                         "variable": "R_O1_equiv_over_sens", "value": "No determinable" if ratio is None else ratio,
                         "units": "-"})
            rows.append({"source": "TAM", "run": case, "period": scope, "period_index": "",
                         "variable": "R_O1_equiv_if_O_conditional_over_sens",
                         "value": summ["R_O1_equiv_if_O_conditional"]["mean"] / sens, "units": "- (condicional)"})
            rows.append({"source": "TAM", "run": case, "period": scope, "period_index": "",
                         "variable": "R_O1_equiv_if_O_eps095_balance_over_sens",
                         "value": summ["R_O1_equiv_if_O_eps095_balance"]["mean"] / sens,
                         "units": "- (condicional, lectura del balance)"})
            meta[case][scope] = summ
        meta[case]["daily_series_MJ_m2"] = {v: (None if arr is None else [float(x) for x in arr])
                                            for v, arr in daily.items()}
    return rows, meta


def main() -> int:
    from titan_pbl.config import load_config

    sol_length = float(load_config("huygens_baseline").planet.sol_length)
    t1 = json.loads((C.R2_OUT / "t1_tam_boundary.json").read_text(encoding="utf-8"))
    a_rows, a_meta = archive_partition(sol_length)
    t_rows, t_meta = tam_partition(t1)
    rows = a_rows + t_rows
    with (C.R2_OUT / "t2_energy_partition.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["source", "run", "period", "period_index", "variable", "value", "units"])
        writer.writeheader()
        writer.writerows(rows)
    sign = {
        "GSW": "+ hacia la superficie (neto): 'NET SHORT WAVE FLUX AT GROUND SURFACE'",
        "GLW": "+ hacia abajo: 'DOWNWARD LONG WAVE FLUX AT GROUND SURFACE'",
        "U": "+ hacia arriba: FLUX_LW_U[:,0], 'Longwave upwelling flux' (indice 0 = superficie, R0)",
        "EMIS": "+ hacia arriba: 0.95 * 5.67051e-8 * TSK^4",
        "XINET": "+ hacia abajo: 'Net downward longwave radiative flux at the surface'",
        "R_O1": "+ hacia la superficie: (1-eps) GLW",
        "HFX": "+ hacia arriba: 'UPWARD HEAT FLUX AT THE SURFACE'",
        "LH": "'LATENT HEAT FLUX AT THE SURFACE'; identicamente 0 en las tres corridas (signo no verificable)",
        "G_residual": "+ hacia el suelo: RESIDUO GSW + XINET - HFX - LH; coincide con FLUX[:,0] 'CONDUCTIVE HEAT FLUX'",
        "SOILCOND": "'Soil Conduction Flux to the surface'; NO cierra el balance con ningun signo (solo descriptivo)",
        "TAM": {"sw_sfc": "+ hacia abajo, neto: 'SW net flux down at surface'",
                "lwdn_sfc": "+ hacia abajo: 'LW flux down at surface'",
                "lwup_sfc": "+ hacia arriba: 'LW flux up at surface'",
                "sens": "'Sensible heat flux' (+ hacia arriba; convencion de FMS, consistente con T1)",
                "hf": "+ hacia el suelo: 'Net flux into the ground'"},
    }
    report = {"task": "R2-T2", "guard": C.guard_info(), "sol_length_s": sol_length, "eps": EPS, "stbolt": STBOLT,
              "sign_conventions": sign, "archive": a_meta, "tam": t_meta}
    C.dump_json(C.R2_OUT / "t2_energy_partition.json", report)

    # numbers_r2.json
    src_csv = C.rel(C.R2_OUT / "t2_energy_partition.csv")
    fn = "code/analysis/r2_t2_energy_partition.py"
    nums = []
    why = "integral por sol de lecturas del archivo; 2 cifras bastan para comparar con HFX"
    for run in C.ARCHIVE_RUNS:
        for sol in SOL_INDICES:
            ps = a_meta[run]["per_sol"][sol]
            tag = run.replace(".", "p")
            nums.append(C.entry(f"energy_archive_HFX_sol{sol}_{tag}", ps["HFX"], "MJ m-2 per Tsol", "n/a", "n/a",
                                src_csv, f"{fn}::archive_partition", 2, why))
            nums.append(C.entry(f"energy_archive_R_O1_sol{sol}_{tag}", ps["R_O1"], "MJ m-2 per Tsol", "full_absorption",
                                "n/a", src_csv, f"{fn}::archive_partition", 2, why))
            nums.append(C.entry(f"energy_archive_R_O1_over_HFX_sol{sol}_{tag}", ps["R_O1_over_HFX"], "-", "full_absorption",
                                "n/a", src_csv, f"{fn}::archive_partition", 2, "cociente de dos integrales de 2 cifras"))
            nums.append(C.entry(f"energy_archive_GSW_sol{sol}_{tag}", ps["GSW"], "MJ m-2 per Tsol", "n/a", "n/a",
                                src_csv, f"{fn}::archive_partition", 2, why))
    for case in TAM_HUYGENS:
        summ = t_meta[case]["30 centrales (k=10..39)"]
        nums.append(C.entry(f"energy_tam_sens_daily_central30_{case}", summ["sens"]["mean"], "MJ m-2 per Titan day",
                            "n/a", "n/a", src_csv, f"{fn}::tam_partition", 2, "media de 30 dias de una lectura de TAM"))
        nums.append(C.entry(f"energy_tam_R_O1_equiv_daily_central30_{case}",
                            None if summ["R_O1_equiv"] is None else summ["R_O1_equiv"]["mean"],
                            "MJ m-2 per Titan day", "n/a", "n/a", src_csv, f"{fn}::tam_partition", 2,
                            "None = No determinable en T1"))
        nums.append(C.entry(f"energy_tam_R_O1_equiv_eps095_balance_daily_central30_{case}",
                            summ["R_O1_equiv_if_O_eps095_balance"]["mean"], "MJ m-2 per Titan day", "n/a", "n/a", src_csv,
                            f"{fn}::tam_partition", 2, "condicional: lectura O1 eps=0.95 del balance superficial (T1, descriptivo)"))
    C.update_numbers_r2(nums, "energy_", fn)

    for run in C.ARCHIVE_RUNS:
        print(run, {k: v for k, v in a_meta[run].items() if k.startswith("check") or k in ("dt_record_s", "has_ZS")})
        for sol, ps in a_meta[run]["per_sol"].items():
            print("  Tsol", sol, {k: round(v, 4) for k, v in ps.items()})
    for case in TAM_HUYGENS:
        s = t_meta[case]["30 centrales (k=10..39)"]
        print(case, t_meta[case]["T1_class"], {k: (None if v is None else round(v["mean"], 4)) for k, v in s.items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
