#!/usr/bin/env python
"""
r2_t3_mayor3_budget.py -- stage r2, analysis T3.

Thermal budget of the archived mtWRF runs,
baseline against suppressed sensible heat flux (0.1flhc), and offline
attribution of the longwave change with the titan_pbl LW scheme verified in
R1/T6.

Cantidades que calcula:
 (a) cuanto aumenta el calentamiento radiativo de la atmosfera baja por la
     mayor emision superficial;
 (b) si su estructura vertical desestabiliza;
 (c) si alcanza a erosionar la estabilidad de fondo en escala diurna;
 (d) si la PBL mas profunda se explica por la radiacion o por otro termino
     del esquema radiativo o turbulento.

---------------------------------------------------------------------------
CRITERIOS Y METODO (declarados antes de ejecutar)
---------------------------------------------------------------------------
1. Emparejamiento: se comparan TODOS los atributos globales entre 0.1flhc y
   cada baseline y se transcribe literalmente toda diferencia, incluida TMN.
   Par principal: 0.1flhc - TMN91. Sensibilidad: 0.1flhc - baseline_5sol.
2. Geometria: reconstruccion de interfaces de R1/T6 (reconstruct_interfaces de
   r1_flux_comparison.py, importada sin modificar). P5: debe reproducir los
   valores registrados en build/r1/t6_flux_comparison.json (TMN91, registros
   87-662) a tolerancia relativa 1e-12.
3. theta en niveles de masa con p = P + PB: theta = T (p_ref/p)^kappa,
   kappa = R/cp, R = 296.8, p_ref = 145000 Pa. cp: cita literal de Hartwick
   (Tabla 2, p. 5: "Atmosphere (N 2 ) 1101.0 (E. W. Lemmon & R. T. Jacobsen
   1998)", columna "Specific Heat Capacity (J kg-1 K-1)") => cp = 1101.0, sin
   marca. Sensibilidad cp +-5 %.
4. Tendencias por nivel y registro, Tsols 2-5 (registros 87-662):
   total = diferencia centrada en el tiempo de theta (dt = 9570 s);
   LW = divergencia del flujo neto LW entre interfaces, g (F[k+1]-F[k]) /
        (cp dp_k), convertida a theta con (p_ref/p_k)^kappa; capas 0-54 (la 55
        no tiene interfaz de tope);
   residuo = total - LW (incluye SW, turbulencia y efectos numericos).
   El archivo SI trae perfiles SW (FLUX_SW_D/U, interfaces 0-55, D0-U0 = GSW):
   la tendencia SW se calcula igual y se reporta aparte; la definicion de
   "residuo" declarado (total - LW) se mantiene para Q_S. Se documenta la
   diferencia de GSW entre corridas.
5. Atribucion offline (pregunta a), maquina LW de R1/T6, full_absorption:
     dLW_emis = H(Ts_s, T_b) - H(Ts_b, T_b)
     dLW_prof = H(Ts_s, T_s) - H(Ts_s, T_b)
     cierre: dLW_emis + dLW_prof = H(Ts_s, T_s) - H(Ts_b, T_b) a precision de
     maquina. Se repite con kirchhoff.
   Malla comun (declarada): la del registro de la baseline (interfaces
   reconstruidas, P+PB y PSFC); T_s se asigna capa a capa. Asi la identidad
   de cierre es exacta. Como control se evalua ademas H(Ts_s, T_s) en la malla
   propia de la suprimida y se compara con el LW archivado.
6. Estabilidad, z_ref en {1000, 2500, 4000} m: S = theta(z_ref) - theta(z1),
   z1 = nivel de masa mas bajo; theta(z_ref) por interpolacion lineal en la
   altura de los niveles de masa (media de los niveles escalonados vecinos,
   sobre la superficie). Se descompone la tasa de cambio de dS entre corridas
   en parte LW y residuo; se integra en los Tsols 2-5.
     Q_S = |int dS_LW dt| / |int dS_total dt|
     Q_S >= 0.5 -> APOYA (el cambio radiativo es el contribuyente principal);
     Q_S <= 0.1 -> NO APOYA (origen no radiativo);
     0.1 < Q_S < 0.5 -> PARCIAL.
   Agregacion entre z_ref (declarada aqui, antes de ejecutar): si los tres
   z_ref dan el mismo veredicto, ese es el veredicto; si no, "depende de
   z_ref", con los tres.
7. Estructura vertical (pregunta b): signo de d(dLW)/dz en 0-1 km (pendiente
   de minimos cuadrados sobre los niveles de masa con z <= 1000 m del perfil
   medio de Tsols 2-5), separando dia (GSW de la baseline > 0) y noche.
   Negativo (mas calentamiento abajo) = desestabiliza. Se evalua para
   dLW_emis offline (primario para (b)) y para el cambio LW archivado.
8. Escala de erosion (pregunta c): t_eros = S_base / <dS_LW> en sols, con
   S_base a las 09:30 LT de cada sol (registro "_09:30:00"). Se reporta el
   valor literal con su signo: <dS_LW> < 0 erosiona (desestabiliza); > 0
   estabiliza y entonces no hay erosion. Con erosion: |t_eros| <= 1 sol ->
   puede actuar en escala diurna; >= 10 sols -> no puede; entre ambos ->
   intermedio.
9. PBLH por sol en las dos corridas, solo descriptivo (cuantizada en niveles).
10. Presupuesto en energia por capa: dE_X = cp dp/g (dT_X/dt) integrado por
    sol, X = LW, total, residuo, en J m-2 sol-1 (diferencia supr - base,
    media de los Tsols 2-5).
11. Three values of Section 5 (added in R4B), staged with the prefix "sec5_":
    - sec5_dTSK_end_sol5: TSK(0.1flhc) - TSK(TMN91) at the last record of
      sol 5 (record 662), in K.
    - sec5_S1km_oscillation_halfamp: half-amplitude of the within-sol
      oscillation of the total stability difference at z_ref = 1 km, i.e. of
      dS_tot(t) = cumulative time integral of [dS/dt(0.1flhc) - dS/dt(TMN91)]
      with S = theta(1 km) - theta(z_1) (the "total" curve of Figure 6a).
      Statistic: for each sol 2-5 (records of SOL_INDICES), (max - min)/2 of
      dS_tot within that sol; the value staged is the MAXIMUM over the four
      sols; the MEDIAN over the sols and the per-sol values are given in the
      note and in t3_summary.json ("sec5").
    - sec5_tau1_height: height above the surface, in km, at which the
      infrared optical depth accumulated from the surface, tau(p_s) - tau(p),
      reaches 1, with the preset huygens_baseline on its 56-layer grid and the
      HASI initial state (heights from state.geopotential_heights); linear
      interpolation in height between the two bracketing interfaces.

Salidas: build/r2/t3_run_pairing.csv, t3_tendency_profiles.csv,
t3_attribution.csv, t3_stability_budget.csv, t3_summary.json, quicklook/t3_*.png
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
from titan_pbl.config import load_config  # noqa: E402
from titan_pbl.core import build_grid, initialize_state  # noqa: E402
from titan_pbl.schemes.radiation import GrayRadiation  # noqa: E402

R_GAS = 296.8
CP = 1101.0
P_REF = 145000.0
GRAVITY = 1.352
FIRST, LAST = 87, 663
SOL_INDICES = {2: (87, 231), 3: (231, 375), 4: (375, 519), 5: (519, 663)}
Z_REFS = (1000.0, 2500.0, 4000.0)
NLW = 55  # capas con tendencia LW (0..54)
PAIRS = {"primary_0.1flhc-TMN91": ("0.1flhc", "TMN91"), "sensitivity_0.1flhc-baseline5sol": ("0.1flhc", "baseline5sol")}
BOUNDARIES = ("full_absorption", "kirchhoff")


# ---------------------------------------------------------------------------
# Lectura y geometria
# ---------------------------------------------------------------------------

def load_run(path: Path, recon) -> dict:
    ds = netCDF4.Dataset(path)
    sl = slice(FIRST - 1, LAST + 1)  # vecinos para diferencias centradas

    def g2(v):
        return np.asarray(ds[v][sl, :, 0, 0], dtype=np.float64)

    def g1(v):
        return np.asarray(ds[v][sl, 0, 0], dtype=np.float64)

    d = {
        "T": g2("T_PHY"), "P": g2("P") + g2("PB"), "Z": g2("PH") / GRAVITY + g2("PHB") / GRAVITY,
        "LWD": g2("FLUX_LW_D"), "LWU": g2("FLUX_LW_U"), "SWD": g2("FLUX_SW_D"), "SWU": g2("FLUX_SW_U"),
        "PSFC": g1("PSFC"), "TSK": g1("TSK"), "GSW": g1("GSW"), "PBLH": g1("PBLH"), "HFX": g1("HFX"),
        "Times": [str(netCDF4.chartostring(ds["Times"][i])) for i in range(FIRST - 1, LAST + 1)],
    }
    d["Z"] = np.asarray((ds["PH"][sl, :, 0, 0] + ds["PHB"][sl, :, 0, 0]) / GRAVITY, dtype=np.float64)
    extra = {}
    for v in ("RTHRATLW", "RTHRATSW"):
        if v in ds.variables:
            extra[v] = g2(v)
    if "ZNW" in ds.variables:
        extra["ZNW"] = np.asarray(ds["ZNW"][sl], dtype=np.float64)
        extra["P_TOP"] = np.asarray(ds["P_TOP"][sl], dtype=np.float64)
    d["extra"] = extra
    d["attrs"] = {k: ds.getncattr(k) for k in ds.ncattrs()}
    d["TMN0"] = float(ds["TMN"][0, 0, 0])
    d["variables"] = sorted(ds.variables)
    ds.close()
    n = d["T"].shape[0]
    d["p_int"] = np.array([recon(d["PSFC"][r], d["Z"][r], d["T"][r], GRAVITY, R_GAS) for r in range(n)])
    d["z_mass"] = 0.5 * (d["Z"][:, :-1] + d["Z"][:, 1:]) - d["Z"][:, :1]
    d["dp"] = d["p_int"][:, :-1] - d["p_int"][:, 1:]
    return d


def p5_check(recon) -> dict:
    """Stop condition P5: the R1/T6 reconstruction reproduces its registered values."""
    reg = json.loads((C.R1_OUT / "t6_flux_comparison.json").read_text(encoding="utf-8"))["geometry"]
    ds = netCDF4.Dataset(C.ARCHIVE_RUNS["TMN91"])
    sl = slice(FIRST, LAST)
    P = np.asarray(ds["P"][sl, :, 0, 0] + ds["PB"][sl, :, 0, 0], dtype=np.float64)
    Z = np.asarray((ds["PH"][sl, :, 0, 0] + ds["PHB"][sl, :, 0, 0]) / GRAVITY, dtype=np.float64)
    T = np.asarray(ds["T_PHY"][sl, :, 0, 0], dtype=np.float64)
    PSFC = np.asarray(ds["PSFC"][sl, 0, 0], dtype=np.float64)
    ds.close()
    p_int = np.array([recon(PSFC[r], Z[r], T[r], GRAVITY, R_GAS) for r in range(PSFC.size)])
    p_eta = np.empty_like(p_int)
    p_eta[:, 0] = PSFC
    for k in range(P.shape[1]):
        p_eta[:, k + 1] = 2.0 * P[:, k] - p_eta[:, k]
    now = {
        "p_int_top_range_Pa": [float(p_int[:, -1].min()), float(p_int[:, -1].max())],
        "hydro_vs_eta_max_abs_diff_Pa": float(np.max(np.abs(p_int - p_eta))),
        "dp_layer0_range_Pa": [float((p_int[:, 0] - p_int[:, 1]).min()), float((p_int[:, 0] - p_int[:, 1]).max())],
        "mass_levels_bracketed_by_interfaces": bool(np.all((p_int[:, :-1] > P) & (P > p_int[:, 1:]))),
    }
    ok = True
    detail = {}
    for k, v in now.items():
        r = reg[k]
        if isinstance(v, list):
            rel = max(abs(a - b) / abs(b) for a, b in zip(v, r))
            good = rel <= 1e-12
        elif isinstance(v, bool):
            rel, good = None, v == r
        else:
            rel = abs(v - r) / abs(r)
            good = rel <= 1e-12
        detail[k] = {"now": v, "registered_R1": r, "rel_diff": rel, "pass": good}
        ok &= good
    if not ok:
        raise SystemExit(f"P5: la reconstruccion de interfaces no reproduce R1: {detail}")
    return {"pass": ok, "detail": detail}


# ---------------------------------------------------------------------------
# Tendencias
# ---------------------------------------------------------------------------

def tendencies(d: dict, cp: float, dt: float) -> dict:
    kappa = R_GAS / cp
    exner = (P_REF / d["P"]) ** kappa
    theta = d["T"] * exner
    # registros 1..n-2 del bloque cargado = 87..662 del archivo
    th_tot = (theta[2:] - theta[:-2]) / (2.0 * dt)
    t_tot = (d["T"][2:] - d["T"][:-2]) / (2.0 * dt)
    core = slice(1, -1)
    fn_lw = d["LWD"][core] - d["LWU"][core]
    fn_sw = d["SWD"][core] - d["SWU"][core]
    dp = d["dp"][core][:, :NLW]
    t_lw = GRAVITY * (fn_lw[:, 1:] - fn_lw[:, :-1]) / (cp * dp)
    t_sw = GRAVITY * (fn_sw[:, 1:] - fn_sw[:, :-1]) / (cp * dp)
    ex = exner[core][:, :NLW]
    return {
        "theta": theta[core], "exner": exner[core],
        "th_tot": th_tot[:, :NLW], "th_lw": t_lw * ex, "th_sw": t_sw * ex,
        "th_res": th_tot[:, :NLW] - t_lw * ex,
        "t_tot": t_tot[:, :NLW], "t_lw": t_lw, "t_sw": t_sw, "t_res": t_tot[:, :NLW] - t_lw,
        "z": d["z_mass"][core][:, :NLW], "z_all": d["z_mass"][core], "dp": dp,
        "gsw": d["GSW"][core], "pblh": d["PBLH"][core], "times": d["Times"][1:-1],
    }


def interp_rows(z: np.ndarray, f: np.ndarray, z0: float) -> np.ndarray:
    return np.array([np.interp(z0, z[r], f[r]) for r in range(z.shape[0])])


def stability_rates(tend: dict, z_ref: float) -> dict:
    z = tend["z"]
    out = {}
    for key in ("th_tot", "th_lw", "th_sw", "th_res"):
        out[key] = interp_rows(z, tend[key], z_ref) - tend[key][:, 0]
    out["S"] = interp_rows(tend["z_all"], tend["theta"], z_ref) - tend["theta"][:, 0]
    return out


def verdict_q(q: float) -> str:
    if q >= 0.5:
        return "apoya"
    if q <= 0.1:
        return "no apoya"
    return "parcial"


def slope_0_1km(z: np.ndarray, f: np.ndarray) -> float:
    m = z <= 1000.0
    return float(np.polyfit(z[m], f[m], 1)[0])


# ---------------------------------------------------------------------------
# Atribucion offline
# ---------------------------------------------------------------------------

def offline_attribution(rec, base: dict, supp: dict, boundary: str, sol: float) -> dict:
    config = load_config("huygens_baseline", overrides={"radiation.surface_lw_boundary": boundary})
    cp = config.atmosphere.gas.specific_heat_cp
    g = config.planet.gravity
    core = slice(1, -1)
    n = base["T"][core].shape[0]
    nl = base["T"].shape[1]
    out = {k: np.empty((n, nl)) for k in ("H_bb", "H_sb", "H_ss", "H_ss_own")}

    def heating(state, scheme, p_int):
        down, up = scheme.longwave_fluxes(state)
        net = down - up
        return g * (net[1:] - net[:-1]) / (cp * (p_int[:-1] - p_int[1:]))

    for r in range(n):
        i = r + 1
        p_int, p_mass = base["p_int"][i], base["P"][i]
        state, scheme = rec.make_state(config, p_int, p_mass, base["T"][i], base["TSK"][i], base["PSFC"][i])
        out["H_bb"][r] = heating(state, scheme, p_int)
        state.skin_temperature = float(supp["TSK"][i])
        out["H_sb"][r] = heating(state, scheme, p_int)
        state.atmospheric_temperature = supp["T"][i].copy()
        out["H_ss"][r] = heating(state, scheme, p_int)
        s_state, s_scheme = rec.make_state(config, supp["p_int"][i], supp["P"][i], supp["T"][i], supp["TSK"][i],
                                           supp["PSFC"][i])
        out["H_ss_own"][r] = heating(s_state, s_scheme, supp["p_int"][i])
    emis = out["H_sb"] - out["H_bb"]
    prof = out["H_ss"] - out["H_sb"]
    total = out["H_ss"] - out["H_bb"]
    closure = float(np.max(np.abs(emis + prof - total)))
    return {"emis": emis, "prof": prof, "total": total, "closure_max_abs_K_s": closure,
            "total_own_grids": out["H_ss_own"] - out["H_bb"], "H_bb": out["H_bb"], "H_ss_own": out["H_ss_own"],
            "cp": cp}


# ---------------------------------------------------------------------------
# Principal
# ---------------------------------------------------------------------------

def sec5_values(runs: dict, tend: dict, dt: float) -> dict:
    """Three values of Section 5 (docstring, item 11)."""
    i_end = 662 - (FIRST - 1)  # registro 662 en el bloque cargado (FIRST - 1 .. LAST)
    d_tsk = float(runs["0.1flhc"]["TSK"][i_end] - runs["TMN91"]["TSK"][i_end])
    rs, rb = stability_rates(tend["0.1flhc"], 1000.0), stability_rates(tend["TMN91"], 1000.0)
    cum = np.cumsum(rs["th_tot"] - rb["th_tot"]) * dt  # curva "total" de la Figura 6a
    half = {}
    for sol_k, (a, bnd) in SOL_INDICES.items():
        seg = cum[a - FIRST:bnd - FIRST]
        half[sol_k] = float((seg.max() - seg.min()) / 2.0)
    config = load_config("huygens_baseline")
    state = initialize_state(config, build_grid(config))
    scheme = GrayRadiation(config)
    scheme.initialize(config, state)
    tau = scheme.infrared_optical_depth(state)  # interfaces, indice 0 = superficie
    z_int, _ = state.geopotential_heights()
    tau_up = tau[0] - tau  # profundidad optica acumulada desde la superficie
    k = int(np.argmax(tau_up >= 1.0))
    z1 = z_int[k - 1] + (1.0 - tau_up[k - 1]) / (tau_up[k] - tau_up[k - 1]) * (z_int[k] - z_int[k - 1])
    return {
        "dTSK_end_sol5_K": d_tsk, "dTSK_record": 662, "dTSK_time": runs["TMN91"]["Times"][i_end],
        "S1km_halfamp_per_sol_K": {str(s): v for s, v in half.items()},
        "S1km_halfamp_max_K": max(half.values()),
        "S1km_halfamp_median_K": float(np.median(list(half.values()))),
        "tau1_height_km": float(z1 / 1000.0),
        "tau1_bracketing_interfaces_m": [float(z_int[k - 1]), float(z_int[k])],
        "n_layers": int(tau.size - 1), "tau_ir_surface": float(tau[0]),
    }


def main() -> int:
    C.R2_OUT.mkdir(parents=True, exist_ok=True)
    C.QUICKLOOK.mkdir(parents=True, exist_ok=True)
    rec = C.load_r1_flux_module()
    p5 = p5_check(rec.reconstruct_interfaces)
    sol = float(load_config("huygens_baseline").planet.sol_length)
    dt = sol / 144.0
    runs = {name: load_run(path, rec.reconstruct_interfaces) for name, path in C.ARCHIVE_RUNS.items()}
    for name in runs:
        if runs[name]["Times"] != runs["TMN91"]["Times"]:
            raise SystemExit(f"Times de {name} no coinciden con TMN91")

    # 1. Emparejamiento --------------------------------------------------
    pairing_rows = []
    keys = sorted(set().union(*[set(r["attrs"]) for r in runs.values()]))
    s_at = runs["0.1flhc"]["attrs"]
    for k in keys:
        vals = {n: runs[n]["attrs"].get(k, "<ausente>") for n in runs}
        dif = {f"difiere_vs_{b}": str(vals["0.1flhc"]) != str(vals[b]) for b in ("TMN91", "baseline5sol")}
        if any(dif.values()):
            pairing_rows.append({"elemento": f"atributo global {k}", **{n: vals[n] for n in runs}, **dif})
    n_diff_attrs = len(pairing_rows)
    pairing_rows.insert(0, {"elemento": "atributos globales (conteo)",
                            **{n: f"{len(runs[n]['attrs'])} atributos" for n in runs},
                            "difiere_vs_TMN91": f"{n_diff_attrs} difieren",
                            "difiere_vs_baseline5sol": f"{n_diff_attrs} difieren"})
    pairing_rows.append({"elemento": "variable TMN (registro 0)", **{n: runs[n]["TMN0"] for n in runs},
                         **{f"difiere_vs_{b}": runs["0.1flhc"]["TMN0"] != runs[b]["TMN0"] for b in ("TMN91", "baseline5sol")}})
    for b in ("TMN91", "baseline5sol"):
        only_s = sorted(set(runs["0.1flhc"]["variables"]) - set(runs[b]["variables"]))
        only_b = sorted(set(runs[b]["variables"]) - set(runs["0.1flhc"]["variables"]))
        pairing_rows.append({"elemento": f"variables solo en 0.1flhc (vs {b})", "0.1flhc": " ".join(only_s),
                             b: " ".join(only_b) or "(ninguna)"})
    n_attrs = {n: len(r["attrs"]) for n, r in runs.items()}
    with (C.R2_OUT / "t3_run_pairing.csv").open("w", newline="", encoding="utf-8") as fh:
        fields = ["elemento", "0.1flhc", "TMN91", "baseline5sol", "difiere_vs_TMN91", "difiere_vs_baseline5sol"]
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in pairing_rows:
            writer.writerow({f: row.get(f, "") for f in fields})

    # 3-4. Tendencias ----------------------------------------------------
    tend = {name: tendencies(r, CP, dt) for name, r in runs.items()}
    day_mask = {name: t["gsw"] > 0.0 for name, t in tend.items()}
    z_ref_profile = tend["TMN91"]["z"].mean(axis=0)
    checks = {}
    # control con RTHRATLW (solo 0.1flhc): cociente RTHRATLW / (dT/dt)_LW = factor de Exner del archivo
    ex = runs["0.1flhc"]["extra"]
    if "RTHRATLW" in ex:
        # (control corregido tras la primera corrida: el ajuste en log fallaba
        #  por cocientes negativos; se usan solo puntos con cociente > 0)
        rth = ex["RTHRATLW"][1:-1][:, :NLW]
        tlw = tend["0.1flhc"]["t_lw"]
        ratio = rth / tlw
        pm = runs["0.1flhc"]["P"][1:-1][:, :NLW]
        big = np.abs(tlw) > 1e-2 * np.abs(tlw).max()
        sel = big & (ratio > 0)
        slope, icpt = np.polyfit(np.log(pm[sel]), np.log(ratio[sel]), 1)
        fit = np.exp(icpt + slope * np.log(pm))
        checks["RTHRATLW_vs_flux_divergence"] = {
            "note": "RTHRATLW esta en la theta propia del archivo; si ambos usan la misma divergencia, el cociente "
                    "RTHRATLW/(dT/dt)_LW es (p0/p)^kappa del archivo, funcion suave de p",
            "n_points_used": int(sel.sum()), "fraction_opposite_sign_among_large": float(np.mean(ratio[big] <= 0)),
            "fit_kappa_archive": float(-slope), "fit_p0_archive_Pa": float(np.exp(-icpt / slope)),
            "rms_rel_dev_of_ratio_from_fit_layers_0_10": float(np.sqrt(np.nanmean(
                np.where(sel[:, :11], ratio[:, :11] / fit[:, :11] - 1.0, np.nan) ** 2))),
            "median_ratio_by_layer_0_10": [float(np.median(ratio[:, k])) for k in range(11)],
            "reference_(1e5/p)^(296.8/1101)_layer0": float(np.median((1.0e5 / pm[:, 0]) ** (R_GAS / CP))),
        }
    if "ZNW" in ex:
        znw, ptop = ex["ZNW"][1:-1], ex["P_TOP"][1:-1]
        psfc = runs["0.1flhc"]["PSFC"][1:-1]
        p_eta = znw * (psfc - ptop)[:, None] + ptop[:, None]
        p_h = runs["0.1flhc"]["p_int"][1:-1]
        checks["hydrostatic_vs_ZNW_PTOP_0.1flhc"] = {
            "max_rel_diff": float(np.max(np.abs(p_h - p_eta) / p_eta)),
            "max_rel_diff_layer0_dp": float(np.max(np.abs((p_h[:, 0] - p_h[:, 1]) / (p_eta[:, 0] - p_eta[:, 1]) - 1.0))),
        }
    gsw_diff = float(np.max(np.abs(runs["0.1flhc"]["GSW"] - runs["TMN91"]["GSW"])))
    sw_prof_diff = float(np.max(np.abs(tend["0.1flhc"]["th_sw"] - tend["TMN91"]["th_sw"]))) * sol
    checks["GSW_max_abs_diff_0.1flhc_vs_TMN91_W_m2"] = gsw_diff
    checks["SW_theta_tendency_max_abs_diff_K_per_sol"] = sw_prof_diff

    # Perfiles de tendencia (medias de los Tsols 2-5; dia/noche)
    prof_rows = []
    for name, t in tend.items():
        for scope, m in (("todo", np.ones_like(day_mask[name])), ("dia", day_mask[name]), ("noche", ~day_mask[name])):
            m = m.astype(bool)
            for k in range(NLW):
                prof_rows.append({"run": name, "scope": scope, "layer": k, "z_mean_m": float(t["z"][:, k].mean()),
                                  **{f"{x}_K_per_sol": float(t[x][m, k].mean() * sol)
                                     for x in ("th_tot", "th_lw", "th_sw", "th_res")}})
    for pname, (s, b) in PAIRS.items():
        for scope in ("todo", "dia", "noche"):
            m = {"todo": np.ones(tend[b]["gsw"].size, bool), "dia": day_mask[b], "noche": ~day_mask[b]}[scope]
            for k in range(NLW):
                prof_rows.append({"run": f"delta {pname}", "scope": scope, "layer": k,
                                  "z_mean_m": float(tend[b]["z"][:, k].mean()),
                                  **{f"{x}_K_per_sol": float((tend[s][x][m, k] - tend[b][x][m, k]).mean() * sol)
                                     for x in ("th_tot", "th_lw", "th_sw", "th_res")}})
    with (C.R2_OUT / "t3_tendency_profiles.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(prof_rows[0]))
        writer.writeheader()
        writer.writerows(prof_rows)

    # 5. Atribucion offline ----------------------------------------------
    attribution = {}
    attr_rows = []
    for pname, (s, b) in PAIRS.items():
        attribution[pname] = {}
        for boundary in BOUNDARIES:
            res = offline_attribution(rec, runs[b], runs[s], boundary, sol)
            exner_b = tend[b]["exner"]
            res["emis_th"] = res["emis"] * exner_b
            res["prof_th"] = res["prof"] * exner_b
            res["total_th"] = res["total"] * exner_b
            attribution[pname][boundary] = res
            dm = day_mask[b]
            for scope, m in (("todo", np.ones_like(dm)), ("dia", dm), ("noche", ~dm)):
                m = m.astype(bool)
                for k in range(NLW):
                    arch = (tend[s]["t_lw"][m, k] - tend[b]["t_lw"][m, k]).mean() * sol
                    attr_rows.append({
                        "pair": pname, "lw_boundary": boundary, "scope": scope, "layer": k,
                        "z_mean_m": float(tend[b]["z"][:, k].mean()),
                        "dLW_emis_K_per_sol": float(res["emis"][m, k].mean() * sol),
                        "dLW_prof_K_per_sol": float(res["prof"][m, k].mean() * sol),
                        "dLW_total_offline_K_per_sol": float(res["total"][m, k].mean() * sol),
                        "closure_residual_K_per_sol": float((res["emis"][m, k] + res["prof"][m, k] - res["total"][m, k]).mean() * sol),
                        "dLW_total_offline_own_grids_K_per_sol": float(res["total_own_grids"][m, k].mean() * sol),
                        "dLW_archive_K_per_sol": float(arch),
                    })
    with (C.R2_OUT / "t3_attribution.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(attr_rows[0]))
        writer.writeheader()
        writer.writerows(attr_rows)

    # 6, 8. Estabilidad, Q_S y t_eros --------------------------------------
    stab_rows, stab = [], {}
    nine30 = {sol_k: next(i - FIRST for i in range(a, bnd) if runs["TMN91"]["Times"][i - FIRST + 1].endswith("_09:30:00"))
              for sol_k, (a, bnd) in SOL_INDICES.items()}
    for pname, (s, b) in PAIRS.items():
        stab[pname] = {}
        for z_ref in Z_REFS:
            rs, rb = stability_rates(tend[s], z_ref), stability_rates(tend[b], z_ref)
            d = {k: rs[k] - rb[k] for k in ("th_tot", "th_lw", "th_sw", "th_res")}
            ints = {k: float(np.sum(v) * dt) for k, v in d.items()}
            q = abs(ints["th_lw"]) / abs(ints["th_tot"]) if ints["th_tot"] != 0 else float("nan")
            # emision offline (full_absorption y kirchhoff) en la misma metrica
            off = {}
            for boundary in BOUNDARIES:
                a = attribution[pname][boundary]
                ze = tend[b]["z"]
                for comp in ("emis_th", "prof_th", "total_th"):
                    rate = interp_rows(ze, a[comp][:, :NLW], z_ref) - a[comp][:, 0]
                    off[f"int_dS_{comp}_{boundary}_K"] = float(np.sum(rate) * dt)
            mean_dlw = float(np.mean(d["th_lw"]))  # K/s
            s_base_0930 = {sk: float(rb["S"][nine30[sk]]) for sk in SOL_INDICES}
            t_eros = {sk: (v / (mean_dlw * sol) if mean_dlw != 0 else float("nan")) for sk, v in s_base_0930.items()}
            erodes = mean_dlw < 0
            mag = np.mean([abs(v) for v in t_eros.values()])
            if not erodes:
                eros_verdict = "no erosiona: <dS_LW> > 0 (el cambio LW estabiliza)"
            elif mag <= 1.0:
                eros_verdict = "puede actuar en escala diurna (|t_eros| <= 1 sol)"
            elif mag >= 10.0:
                eros_verdict = "no puede actuar en escala diurna (|t_eros| >= 10 sols)"
            else:
                eros_verdict = "intermedio (1 < |t_eros| < 10 sols)"
            # dS total observada entre corridas (cheque de la integral)
            ds_obs = float((rs["S"][-1] - rs["S"][0]) - (rb["S"][-1] - rb["S"][0]))
            row = {
                "pair": pname, "z_ref_m": z_ref,
                "int_dS_total_K": ints["th_tot"], "int_dS_LW_K": ints["th_lw"], "int_dS_SW_K": ints["th_sw"],
                "int_dS_res_K": ints["th_res"], "Q_S": q, "verdict": verdict_q(q),
                # descriptivo anadido tras la primera corrida (no altera el veredicto)
                "sign_LW_equals_sign_total": bool(np.sign(ints["th_lw"]) == np.sign(ints["th_tot"])),
                "check_dS_end_minus_start_K": ds_obs,
                "mean_dS_LW_K_per_sol": mean_dlw * sol,
                **{f"S_base_0930_sol{sk}_K": v for sk, v in s_base_0930.items()},
                **{f"t_eros_sol{sk}_sols": v for sk, v in t_eros.items()},
                "t_eros_verdict": eros_verdict,
                **off,
                "Q_S_offline_emis_full_absorption": abs(off["int_dS_emis_th_full_absorption_K"]) / abs(ints["th_tot"]),
                "Q_S_offline_emis_kirchhoff": abs(off["int_dS_emis_th_kirchhoff_K"]) / abs(ints["th_tot"]),
            }
            stab_rows.append(row)
            stab[pname][z_ref] = {"row": row, "series": d, "S_s": rs["S"], "S_b": rb["S"]}
    with (C.R2_OUT / "t3_stability_budget.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(stab_rows[0]))
        writer.writeheader()
        writer.writerows(stab_rows)

    # Sensibilidad a cp (+-5 %): solo cantidades del archivo
    cp_sens = {}
    for factor in (0.95, 1.0, 1.05):
        cpv = CP * factor
        tt = {n: tendencies(runs[n], cpv, dt) for n in ("0.1flhc", "TMN91")}
        cp_sens[f"cp={cpv:.2f}"] = {}
        for z_ref in Z_REFS:
            rs, rb = stability_rates(tt["0.1flhc"], z_ref), stability_rates(tt["TMN91"], z_ref)
            i_tot = float(np.sum(rs["th_tot"] - rb["th_tot"]) * dt)
            i_lw = float(np.sum(rs["th_lw"] - rb["th_lw"]) * dt)
            mean_dlw = float(np.mean(rs["th_lw"] - rb["th_lw"]))
            cp_sens[f"cp={cpv:.2f}"][str(int(z_ref))] = {
                "Q_S": abs(i_lw) / abs(i_tot), "verdict": verdict_q(abs(i_lw) / abs(i_tot)),
                "t_eros_sol3_sols": float(rb["S"][nine30[3]] / (mean_dlw * sol)),
            }

    # 7. Estructura vertical ----------------------------------------------
    vertical = {}
    for pname, (s, b) in PAIRS.items():
        vertical[pname] = {}
        zb = tend[b]["z"].mean(axis=0)
        for scope in ("todo", "dia", "noche"):
            m = {"todo": np.ones(tend[b]["gsw"].size, bool), "dia": day_mask[b], "noche": ~day_mask[b]}[scope]
            entry = {}
            for boundary in BOUNDARIES:
                a = attribution[pname][boundary]
                for comp in ("emis", "prof", "total"):
                    prof = a[comp][m, :NLW].mean(axis=0) * sol
                    sl = slope_0_1km(zb, prof)
                    entry[f"slope_dLW_{comp}_{boundary}_K_per_sol_per_km"] = sl * 1000.0
                    entry[f"sign_dLW_{comp}_{boundary}"] = "negativo (desestabiliza)" if sl < 0 else "positivo (estabiliza)"
            arch = (tend[s]["t_lw"][m] - tend[b]["t_lw"][m]).mean(axis=0) * sol
            sl = slope_0_1km(zb, arch)
            entry["slope_dLW_archive_K_per_sol_per_km"] = sl * 1000.0
            entry["sign_dLW_archive"] = "negativo (desestabiliza)" if sl < 0 else "positivo (estabiliza)"
            entry["n_levels_0_1km"] = int(np.sum(zb <= 1000.0))
            vertical[pname][scope] = entry

    # 9. PBLH -------------------------------------------------------------
    pblh = {}
    for name, t in tend.items():
        pblh[name] = {}
        for sk, (a, bnd) in SOL_INDICES.items():
            v = t["pblh"][a - FIRST:bnd - FIRST]
            pblh[name][sk] = {"max_m": float(v.max()), "mean_m": float(v.mean()),
                              "n_distinct_values": int(np.unique(np.round(v, 3)).size)}

    # 10. Energia por capa -------------------------------------------------
    energy = {}
    for pname, (s, b) in PAIRS.items():
        e = {}
        for x in ("t_lw", "t_tot", "t_res", "t_sw"):
            es = CP * tend[s]["dp"] / GRAVITY * tend[s][x]
            eb = CP * tend[b]["dp"] / GRAVITY * tend[b][x]
            e[x] = ((es - eb).mean(axis=0) * sol)
        zb = tend[b]["z"].mean(axis=0)
        energy[pname] = {
            "per_layer_J_m2_per_sol": {x: v.tolist() for x, v in e.items()},
            "z_mean_m": zb.tolist(),
            "column_0_1km_J_m2_per_sol": {x: float(v[zb <= 1000.0].sum()) for x, v in e.items()},
            "column_0_4km_J_m2_per_sol": {x: float(v[zb <= 4000.0].sum()) for x, v in e.items()},
            "column_all_0_54_J_m2_per_sol": {x: float(v.sum()) for x, v in e.items()},
        }

    # Resumen de (a): aumento del calentamiento LW de la atmosfera baja
    qa = {}
    for pname in PAIRS:
        qa[pname] = {}
        b = PAIRS[pname][1]
        zb = tend[b]["z"].mean(axis=0)
        for boundary in BOUNDARIES:
            a = attribution[pname][boundary]
            for scope in ("todo", "dia", "noche"):
                m = {"todo": np.ones(tend[b]["gsw"].size, bool), "dia": day_mask[b], "noche": ~day_mask[b]}[scope]
                emis = a["emis"][m, :NLW].mean(axis=0) * sol
                prof = a["prof"][m, :NLW].mean(axis=0) * sol
                low = zb <= 1000.0
                w = tend[b]["dp"][m][:, :NLW].mean(axis=0)
                qa[pname][f"{boundary}_{scope}"] = {
                    "layer0_dLW_emis_K_per_sol": float(emis[0]),
                    "layer0_dLW_prof_K_per_sol": float(prof[0]),
                    "mass_weighted_0_1km_dLW_emis_K_per_sol": float(np.sum(emis[low] * w[low]) / np.sum(w[low])),
                    "mass_weighted_0_1km_dLW_prof_K_per_sol": float(np.sum(prof[low] * w[low]) / np.sum(w[low])),
                    "closure_max_abs_K_per_sol": a["closure_max_abs_K_s"] * sol,
                }
        # comparacion del cambio total offline (mallas propias) con el archivado
        s = PAIRS[pname][0]
        arch = (tend[s]["t_lw"] - tend[b]["t_lw"]) * sol
        offl = attribution[pname]["full_absorption"]["total_own_grids"][:, :NLW] * sol
        qa[pname]["offline_vs_archive_dLW_rms_K_per_sol_layers_0_10"] = float(np.sqrt(np.mean((offl[:, :11] - arch[:, :11]) ** 2)))
        qa[pname]["archive_dLW_layer0_mean_K_per_sol"] = float(arch[:, 0].mean())
        qa[pname]["offline_dLW_layer0_mean_K_per_sol"] = float(offl[:, 0].mean())
        qa[pname]["delta_total_theta_tendency_layer0_K_per_sol"] = float((tend[s]["th_tot"][:, 0] - tend[b]["th_tot"][:, 0]).mean() * sol)
        qa[pname]["delta_residual_theta_tendency_layer0_K_per_sol"] = float((tend[s]["th_res"][:, 0] - tend[b]["th_res"][:, 0]).mean() * sol)

    verdicts = {p: [stab[p][z]["row"]["verdict"] for z in Z_REFS] for p in PAIRS}
    agg = {p: (v[0] if len(set(v)) == 1 else "depende de z_ref") for p, v in verdicts.items()}
    summary = {
        "task": "R2-T3", "guard": C.guard_info(),
        "constants": {"R": R_GAS, "cp": CP, "cp_source": "Hartwick et al. 2026, Tabla 2 (p. 5): 'Atmosphere (N 2 ) 1101.0'",
                      "p_ref": P_REF, "kappa": R_GAS / CP, "g": GRAVITY, "dt_record_s": dt, "sol_s": sol},
        "P5_interface_reconstruction": p5,
        "pairing": {"n_global_attributes": n_attrs, "rows_with_differences": len(pairing_rows)},
        "checks": checks,
        "question_a_lw_increase": qa,
        "question_b_vertical_structure": vertical,
        "question_c_erosion": {p: {str(int(z)): {k: stab[p][z]["row"][k] for k in stab[p][z]["row"]
                                                 if k.startswith(("t_eros", "S_base", "mean_dS_LW"))} for z in Z_REFS}
                               for p in PAIRS},
        "question_d_Q_S": {p: {str(int(z)): {"Q_S": stab[p][z]["row"]["Q_S"], "verdict": stab[p][z]["row"]["verdict"],
                                             "int_dS_total_K": stab[p][z]["row"]["int_dS_total_K"],
                                             "int_dS_LW_K": stab[p][z]["row"]["int_dS_LW_K"],
                                             "int_dS_SW_K": stab[p][z]["row"]["int_dS_SW_K"],
                                             "int_dS_res_K": stab[p][z]["row"]["int_dS_res_K"]}
                               for z in Z_REFS} for p in PAIRS},
        "verdict_by_pair": agg,
        "cp_sensitivity": cp_sens,
        "pblh_per_sol": pblh,
        "energy_budget": {p: {k: v for k, v in e.items() if k != "per_layer_J_m2_per_sol"} for p, e in energy.items()},
        "z_mean_lowest_levels_m": z_ref_profile[:12].tolist(),
    }
    # 11. Tres valores de la Seccion 5 (R4B) ---------------------------------
    sec5 = sec5_values(runs, tend, dt)
    summary["sec5"] = sec5
    C.dump_json(C.R2_OUT / "t3_summary.json", summary)
    C.dump_json(C.R2_OUT / "t3_energy_per_layer.json", energy)

    # numbers_r2.json
    src = C.rel(C.R2_OUT / "t3_summary.json")
    fn = "code/analysis/r2_t3_mayor3_budget.py::main"
    nums = []
    p = "primary_0.1flhc-TMN91"
    for z in Z_REFS:
        row = stab[p][z]["row"]
        zt = int(z)
        nums.append(C.entry(f"mayor3_Q_S_{zt}m", row["Q_S"], "-", "n/a (archivo)", "n/a", src, fn, 2,
                            "cociente de dos integrales del archivo; 2 cifras bastan frente a los cortes 0.1/0.5"))
        nums.append(C.entry(f"mayor3_int_dS_total_{zt}m", row["int_dS_total_K"], "K", "n/a (archivo)", "n/a", src, fn, 2, "idem"))
        nums.append(C.entry(f"mayor3_int_dS_LW_{zt}m", row["int_dS_LW_K"], "K", "n/a (archivo)", "n/a", src, fn, 2, "idem"))
        nums.append(C.entry(f"mayor3_t_eros_sol3_{zt}m", row["t_eros_sol3_sols"], "Titan sols", "n/a (archivo)", "n/a", src, fn, 2,
                            "cociente de lecturas del archivo; signo negativo = erosiona"))
    for boundary in BOUNDARIES:
        v = qa[p][f"{boundary}_todo"]
        nums.append(C.entry(f"mayor3_layer0_dLW_emis_{boundary}", v["layer0_dLW_emis_K_per_sol"], "K per Titan sol", boundary, "n/a (solo LW)",
                            src, fn, 2, "atribucion offline con la maquina LW de R1/T6"))
        nums.append(C.entry(f"mayor3_mass_weighted_0_1km_dLW_emis_{boundary}", v["mass_weighted_0_1km_dLW_emis_K_per_sol"],
                            "K per Titan sol", boundary, "n/a (solo LW)", src, fn, 2, "idem"))
        nums.append(C.entry(f"mayor3_slope_0_1km_dLW_emis_{boundary}",
                            vertical[p]["todo"][f"slope_dLW_emis_{boundary}_K_per_sol_per_km"], "K per Titan sol per km",
                            boundary, "n/a (solo LW)", src, fn, 2, "pendiente de minimos cuadrados 0-1 km"))
    nums.append(C.entry("mayor3_verdict_primary", agg[p], "texto", "n/a", "n/a", src, fn, None, "veredicto segun criterio declarado"))
    e = C.entry("sec5_dTSK_end_sol5", sec5["dTSK_end_sol5_K"], "K", "n/a (archivo)", "n/a", src, fn, 2,
                "archive reading")
    e["note"] = f"TSK(0.1flhc) - TSK(TMN91) at record {sec5['dTSK_record']} ({sec5['dTSK_time']}), last record of sol 5"
    nums.append(e)
    e = C.entry("sec5_S1km_oscillation_halfamp", sec5["S1km_halfamp_max_K"], "K", "n/a (archivo)", "n/a", src, fn, 1,
                "amplitude of an oscillation")
    per_sol = ", ".join(f"sol {s}: {v:.4g}" for s, v in sec5["S1km_halfamp_per_sol_K"].items())
    e["note"] = ("max over sols 2-5 of (max - min)/2 within each sol of the cumulative total stability difference "
                 f"at 1 km (Figure 6a); median over sols = {sec5['S1km_halfamp_median_K']:.4g} K; {per_sol} K")
    nums.append(e)
    e = C.entry("sec5_tau1_height", sec5["tau1_height_km"], "km", "n/a", "n/a", src, fn, 2,
                "grid geometry")
    e["note"] = ("height where the infrared optical depth accumulated from the surface reaches 1 "
                 f"(preset huygens_baseline, {sec5['n_layers']} layers, HASI initial state; interpolated in height)")
    nums.append(e)
    C.update_numbers_r2(nums, ("mayor3_", "sec5_"), fn)

    # Graficas de inspeccion
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    b, s = "TMN91", "0.1flhc"
    zb = tend[b]["z"].mean(axis=0)
    top = zb <= 5000.0
    fig, axes = plt.subplots(1, 3, figsize=(13, 5), sharey=True)
    for ax, scope in zip(axes, ("todo", "dia", "noche")):
        m = {"todo": np.ones(tend[b]["gsw"].size, bool), "dia": day_mask[b], "noche": ~day_mask[b]}[scope]
        for x, lab in (("th_tot", "total"), ("th_lw", "LW"), ("th_sw", "SW"), ("th_res", "residuo")):
            ax.plot(((tend[s][x][m] - tend[b][x][m]).mean(axis=0) * sol)[top], zb[top], label=lab)
        ax.axvline(0, color="k", lw=0.5)
        ax.set_title(f"supr - TMN91, {scope}", fontsize=9)
        ax.set_xlabel("d(theta)/dt difference [K/sol]")
    axes[0].set_ylabel("z [m]")
    axes[0].legend(fontsize=8)
    fig.suptitle("R2/T3 quicklook: tendency differences (not a final figure)", fontsize=9)
    fig.savefig(C.QUICKLOOK / "t3_tendency_differences.png", dpi=110)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    for ax, boundary in zip(axes, BOUNDARIES):
        a = attribution["primary_0.1flhc-TMN91"][boundary]
        for comp, lab in (("emis", "emision Ts"), ("prof", "perfil T"), ("total", "total offline")):
            for scope, ls in (("dia", "-"), ("noche", "--")):
                m = day_mask[b] if scope == "dia" else ~day_mask[b]
                ax.plot((a[comp][m, :NLW].mean(axis=0) * sol)[top], zb[top], ls, label=f"{lab} {scope}")
        arch = ((tend[s]["t_lw"] - tend[b]["t_lw"]).mean(axis=0) * sol)[top]
        ax.plot(arch, zb[top], "k:", label="archivo (todo)")
        ax.axvline(0, color="k", lw=0.5)
        ax.set_title(boundary, fontsize=9)
        ax.set_xlabel("dLW heating [K/sol]")
    axes[0].set_ylabel("z [m]")
    axes[0].legend(fontsize=7)
    fig.suptitle("R2/T3 quicklook: offline LW attribution (not a final figure)", fontsize=9)
    fig.savefig(C.QUICKLOOK / "t3_attribution.png", dpi=110)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(13, 4), sharex=True)
    tt = np.arange(FIRST, LAST) / 144.0
    for ax, z in zip(axes, Z_REFS):
        ser = stab["primary_0.1flhc-TMN91"][z]["series"]
        for k, lab in (("th_tot", "total"), ("th_lw", "LW"), ("th_res", "residuo")):
            ax.plot(tt, np.cumsum(ser[k]) * dt, label=lab)
        ax.set_title(f"cumulative dS, z_ref = {int(z)} m", fontsize=9)
        ax.set_xlabel("record / 144 (sols)")
    axes[0].set_ylabel("K")
    axes[0].legend(fontsize=8)
    fig.savefig(C.QUICKLOOK / "t3_cumulative_dS.png", dpi=110)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 3.5))
    for name in runs:
        ax.plot(tt, tend[name]["pblh"], label=name)
    ax.set_ylabel("PBLH [m]")
    ax.set_xlabel("record / 144 (sols)")
    ax.legend(fontsize=8)
    fig.savefig(C.QUICKLOOK / "t3_pblh.png", dpi=110)
    plt.close(fig)

    print(json.dumps({"P5": p5["pass"], "checks": checks, "verdict_by_pair": agg,
                      "Q_S": {p: {str(int(z)): round(stab[p][z]["row"]["Q_S"], 4) for z in Z_REFS} for p in PAIRS},
                      "question_a": qa["primary_0.1flhc-TMN91"],
                      "question_b": vertical["primary_0.1flhc-TMN91"]["todo"],
                      "cp_sens": cp_sens}, indent=1, ensure_ascii=False, default=float))
    for row in stab_rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
