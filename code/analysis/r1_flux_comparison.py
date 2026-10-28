#!/usr/bin/env python
"""
r1_flux_comparison.py -- stage r1, analysis T6.

Replaces the retired 0.9 % agreement (H-28): feeds the titan_pbl longwave
scheme with the grid, the temperature profile and the skin temperature of
each archived mtWRF record, and compares FLUXES (not heating rates) level by
level.

---------------------------------------------------------------------------
GEOMETRIA DEL ARCHIVO (determinada antes de comparar)
---------------------------------------------------------------------------
- FLUX_LW_D/U tienen 56 entradas (dimension bottom_top). Evidencia de que son
  INTERFACES 0..55 (la base de cada capa), sin la interfaz del tope:
  (1) FLUX_LW_D[:,0] = GLW exactamente en todas las corridas (R0/T3), es decir,
      el indice 0 es la superficie, no el primer nivel de masa;
  (2) FLUX_LW_D[:,55] es mayor que I_IR = 2 W m-2 y el perfil decrece
      monotonamente hacia ese valor, como corresponde a la interfaz 55 si la
      56 (tope) vale I_IR. Se registra ademas el I_IR implicito.
- El archivo no trae ZNW, P_TOP ni MU. Las presiones de interfaz se
  reconstruyen por relacion hidrostatica:
      p_int[0] = PSFC,  p_int[k+1] = p_int[k] * exp(-g dz_k / (R T_PHY[k])),
  con dz_k de (PH+PHB)/g en los niveles escalonados, g = 1.352 m s-2 y
  R = 296.8 J kg-1 K-1 (preset). Como control se compara con la reconstruccion
  alternativa de coordenada eta, p_int[k+1] = 2 p_masa[k] - p_int[k].

---------------------------------------------------------------------------
CONFIGURACION PRIMARIA (declarada antes de ejecutar)
---------------------------------------------------------------------------
- Esquema LW de titan_pbl con los parametros del preset huygens_baseline:
  tau_IR = 10, f = 0.15, I_IR = 2 W m-2, eps = 0.95, frontera full_absorption.
- Malla: interfaces reconstruidas; niveles de masa = P + PB del registro;
  T de capa = T_PHY; Ts = TSK; p_s de la Ec. A13 = PSFC del registro.
- Registros: todos los validos de los Tsols 2-5 de la corrida baseline TMN91
  (indices 87 a 662 del archivo, SOL_INDICES del script de figuras).

METRICAS: por interfaz k = 0..55 y por flujo (F_up, F_down): RMS relativo
sqrt(mean(((F_titan - F_arch)/F_arch)^2)) sobre los registros, maximo
absoluto, diferencia media; y la tasa de calentamiento LW de las capas 0..54,
calculada IGUAL para ambos: g (net[k+1] - net[k]) / (cp dp_k) * Tsol, con las
interfaces reconstruidas.

CRITERIO DECLARADO: "concordancia" si el RMS relativo de F_up y de F_down es
<= 1 % en TODAS las interfaces. Si no, se reporta la magnitud tal cual.

VARIANTES DIAGNOSTICAS (se reportan todas; ninguna sustituye a la primaria):
  (i) p_s fijo = 145000 Pa en la Ec. A13; (ii) I_IR = 0; (iii) frontera
  kirchhoff. Sirven solo para atribuir causas. Sin ajuste de parametros.

Salidas: build/r1/t6_flux_comparison.csv y build/r1/t6_flux_comparison.json.
"""

from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import netCDF4  # noqa: E402
import numpy as np  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT / "code" / "titan_pbl_subset") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "code" / "titan_pbl_subset"))

import titan_pbl  # noqa: E402

if (REPO_ROOT / "code" / "titan_pbl_subset") not in Path(titan_pbl.__file__).resolve().parents:
    raise ImportError(f"titan_pbl imported from {titan_pbl.__file__}, not from code/titan_pbl_subset.")

from titan_pbl.config import load_config  # noqa: E402
from titan_pbl.core.grid import AtmosphericGrid, ModelGrid, build_subsurface_grid  # noqa: E402
from titan_pbl.core.state import initialize_state  # noqa: E402
from titan_pbl.schemes.radiation import GrayRadiation  # noqa: E402

_RUN = Path(os.environ.get("TITAN_DATA_DIR", str(REPO_ROOT / "data"))) / "Dry Simulations" / "1D_Huygens_TMN91_5sol"
# Nombre publicado (con ':') o el que deja NTFS al descomprimir (con '_').
ARCHIVE = (_RUN / "wrfout_d01_2007-00001_09:30:00" if (_RUN / "wrfout_d01_2007-00001_09:30:00").exists()
           else _RUN / "wrfout_d01_2007-00001_09_30_00")
OUTPUT_DIR = REPO_ROOT / "build" / "r1"
FIRST, LAST = 87, 663  # Tsols 2-5 (SOL_INDICES 2..5 del script de figuras)
CRITERION = 0.01
FIXED_PS = 145000.0

VARIANTS = {
    "primaria": {},
    "i_ps_fijo_145000": {"fixed_ps": True},
    "ii_IIR_0": {"overrides": {"radiation.ir_top_flux": 0.0}},
    "iii_kirchhoff": {"overrides": {"radiation.surface_lw_boundary": "kirchhoff"}},
}


def reconstruct_interfaces(psfc, z_stag, t_layer, gravity, r_gas):
    p = np.empty(z_stag.size, dtype=np.float64)
    p[0] = psfc
    dz = np.diff(z_stag)
    for k in range(dz.size):
        p[k + 1] = p[k] * np.exp(-gravity * dz[k] / (r_gas * t_layer[k]))
    return p


def make_state(config, p_int, p_mass, t_layer, tsk, surface_pressure):
    eta_int = (p_int - p_int[-1]) / (p_int[0] - p_int[-1])
    eta_mid = (p_mass - p_int[-1]) / (p_int[0] - p_int[-1])
    atmosphere = AtmosphericGrid(
        n_layers=p_mass.size,
        eta_interface=eta_int,
        eta_midpoint=eta_mid,
        pressure_interface=p_int,
        pressure_midpoint=p_mass,
        surface_pressure=surface_pressure,
        top_pressure=float(p_int[-1]),
        spacing_type="archive",
        spacing_exponent=1.0,
    )
    grid = ModelGrid(atmosphere=atmosphere, subsurface=build_subsurface_grid(config),
                     metadata={"config_name": "archive_TMN91"})
    state = initialize_state(config, grid)
    state.atmospheric_temperature = np.asarray(t_layer, dtype=np.float64).copy()
    state.skin_temperature = float(tsk)
    scheme = GrayRadiation(config)
    scheme.initialize(config, state)
    return state, scheme


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    base = load_config("huygens_baseline")
    gravity = base.planet.gravity
    r_gas = base.atmosphere.gas.specific_gas_constant
    cp = base.atmosphere.gas.specific_heat_cp
    sol = base.planet.sol_length

    ds = netCDF4.Dataset(ARCHIVE)
    idx = np.arange(FIRST, LAST)
    D = np.asarray(ds["FLUX_LW_D"][FIRST:LAST, :, 0, 0], dtype=np.float64)
    U = np.asarray(ds["FLUX_LW_U"][FIRST:LAST, :, 0, 0], dtype=np.float64)
    GLW = np.asarray(ds["GLW"][FIRST:LAST, 0, 0], dtype=np.float64)
    P = np.asarray(ds["P"][FIRST:LAST, :, 0, 0] + ds["PB"][FIRST:LAST, :, 0, 0], dtype=np.float64)
    Z = np.asarray((ds["PH"][FIRST:LAST, :, 0, 0] + ds["PHB"][FIRST:LAST, :, 0, 0]) / gravity, dtype=np.float64)
    T = np.asarray(ds["T_PHY"][FIRST:LAST, :, 0, 0], dtype=np.float64)
    TSK = np.asarray(ds["TSK"][FIRST:LAST, 0, 0], dtype=np.float64)
    PSFC = np.asarray(ds["PSFC"][FIRST:LAST, 0, 0], dtype=np.float64)
    times = [str(netCDF4.chartostring(ds["Times"][i])) for i in (FIRST, LAST - 1)]
    ds.close()
    n_rec, n_lev = D.shape

    valid = ~((D[:, 0] == 0.0) & (U[:, 0] == 0.0))
    geometry = {
        "archive": str(ARCHIVE),
        "records": [int(FIRST), int(LAST - 1)],
        "time_strings_first_last": times,
        "n_records": int(n_rec),
        "n_valid": int(valid.sum()),
        "flux_entries": int(n_lev),
        "max_abs_FLUX_LW_D0_minus_GLW": float(np.max(np.abs(D[:, 0] - GLW))),
        "FLUX_LW_D55_range": [float(D[:, 55].min()), float(D[:, 55].max())],
        "FLUX_LW_D_monotonic_fraction": float(np.mean(np.all(np.diff(D, axis=1) < 0, axis=1))),
        "R_gas_used": r_gas,
        "gravity_used": gravity,
    }

    # Presiones de interfaz y control con la reconstruccion eta.
    p_int_all = np.array([reconstruct_interfaces(PSFC[r], Z[r], T[r], gravity, r_gas) for r in range(n_rec)])
    p_eta = np.empty_like(p_int_all)
    p_eta[:, 0] = PSFC
    for k in range(n_lev):
        p_eta[:, k + 1] = 2.0 * P[:, k] - p_eta[:, k]
    bracket_ok = np.all((p_int_all[:, :-1] > P) & (P > p_int_all[:, 1:]))
    geometry.update({
        "p_int_top_range_Pa": [float(p_int_all[:, -1].min()), float(p_int_all[:, -1].max())],
        "hydro_vs_eta_max_abs_diff_Pa": float(np.max(np.abs(p_int_all - p_eta))),
        "hydro_vs_eta_max_rel_diff": float(np.max(np.abs(p_int_all - p_eta) / p_eta)),
        "mass_levels_bracketed_by_interfaces": bool(bracket_ok),
        "dp_layer0_range_Pa": [float((p_int_all[:, 0] - p_int_all[:, 1]).min()), float((p_int_all[:, 0] - p_int_all[:, 1]).max())],
    })

    rows, summary = [], {}
    for name, spec in VARIANTS.items():
        config = load_config("huygens_baseline", overrides=spec.get("overrides") or None)
        Ft_d = np.empty((n_rec, n_lev + 1))
        Ft_u = np.empty((n_rec, n_lev + 1))
        implied_top = np.empty(n_rec)
        for r in range(n_rec):
            ps = FIXED_PS if spec.get("fixed_ps") else PSFC[r]
            state, scheme = make_state(config, p_int_all[r], P[r], T[r], TSK[r], ps)
            down, up = scheme.longwave_fluxes(state)
            Ft_d[r], Ft_u[r] = down, up
            # I_IR implicito del archivo: invierte la ultima capa con la tau del esquema.
            tau = scheme.infrared_optical_depth(state)
            t = np.exp(-(tau[55] - tau[56]))
            b = 5.670374419e-8 * T[r, 55] ** 4
            implied_top[r] = (D[r, 55] - b * (1.0 - t)) / t
        Ft_d, Ft_u = Ft_d[valid], Ft_u[valid]
        Da, Ua = D[valid], U[valid]
        rel_d = (Ft_d[:, :n_lev] - Da) / Da
        rel_u = (Ft_u[:, :n_lev] - Ua) / Ua
        rms_d = np.sqrt(np.mean(rel_d**2, axis=0))
        rms_u = np.sqrt(np.mean(rel_u**2, axis=0))
        for k in range(n_lev):
            rows.append({
                "variante": name, "interfaz": k,
                "rms_rel_F_down": rms_d[k], "rms_rel_F_up": rms_u[k],
                "max_abs_diff_F_down_W_m2": float(np.max(np.abs(Ft_d[:, k] - Da[:, k]))),
                "max_abs_diff_F_up_W_m2": float(np.max(np.abs(Ft_u[:, k] - Ua[:, k]))),
                "mean_diff_F_down_W_m2": float(np.mean(Ft_d[:, k] - Da[:, k])),
                "mean_diff_F_up_W_m2": float(np.mean(Ft_u[:, k] - Ua[:, k])),
                "mean_F_down_archive_W_m2": float(np.mean(Da[:, k])),
                "mean_F_up_archive_W_m2": float(np.mean(Ua[:, k])),
            })
        # Tasas de calentamiento LW de las capas 0..54, calculadas igual.
        p_int_v = p_int_all[valid]
        dp = p_int_v[:, :n_lev - 1] - p_int_v[:, 1:n_lev]
        net_a = Da - Ua
        net_t = Ft_d[:, :n_lev] - Ft_u[:, :n_lev]
        h_a = gravity * (net_a[:, 1:] - net_a[:, :-1]) / (cp * dp) * sol
        h_t = gravity * (net_t[:, 1:] - net_t[:, :-1]) / (cp * dp) * sol
        summary[name] = {
            "max_rms_rel_F_down": float(rms_d.max()),
            "interface_of_max_F_down": int(rms_d.argmax()),
            "max_rms_rel_F_up": float(rms_u.max()),
            "interface_of_max_F_up": int(rms_u.argmax()),
            "concordance_1pct_all_interfaces": bool(rms_d.max() <= CRITERION and rms_u.max() <= CRITERION),
            "n_interfaces_over_1pct_F_down": int(np.sum(rms_d > CRITERION)),
            "n_interfaces_over_1pct_F_up": int(np.sum(rms_u > CRITERION)),
            "surface_mean_diff_F_down_W_m2": float(np.mean(Ft_d[:, 0] - Da[:, 0])),
            "surface_mean_diff_F_up_W_m2": float(np.mean(Ft_u[:, 0] - Ua[:, 0])),
            "surface_rms_rel_F_down": float(rms_d[0]),
            "surface_rms_rel_F_up": float(rms_u[0]),
            "top_stored_interface": 55,
            "top_mean_diff_F_down_W_m2": float(np.mean(Ft_d[:, 55] - Da[:, 55])),
            "top_mean_diff_F_up_W_m2": float(np.mean(Ft_u[:, 55] - Ua[:, 55])),
            "top_rms_rel_F_down": float(rms_d[55]),
            "top_rms_rel_F_up": float(rms_u[55]),
            "implied_archive_I_IR_W_m2_mean": float(np.mean(implied_top[valid])),
            "implied_archive_I_IR_W_m2_range": [float(implied_top[valid].min()), float(implied_top[valid].max())],
            "heating_layer0_archive_K_per_sol_mean": float(np.mean(h_a[:, 0])),
            "heating_layer0_titan_K_per_sol_mean": float(np.mean(h_t[:, 0])),
            "heating_layer0_rms_rel": float(np.sqrt(np.mean(((h_t[:, 0] - h_a[:, 0]) / h_a[:, 0]) ** 2))),
            "heating_all_layers_rms_abs_K_per_sol": float(np.sqrt(np.mean((h_t - h_a) ** 2))),
            "heating_all_layers_max_abs_K_per_sol": float(np.max(np.abs(h_t - h_a))),
        }

    with (OUTPUT_DIR / "t6_flux_comparison.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "task": "R1-T6",
        "criterion": "concordancia si RMS relativo de F_up y F_down <= 1 % en todas las interfaces",
        "primary_configuration": "preset huygens_baseline (tau_IR=10, f=0.15, I_IR=2, eps=0.95, full_absorption); malla, T_PHY, TSK y PSFC del registro",
        "geometry": geometry,
        "summary": summary,
        "verdict_primary": "concordancia" if summary["primaria"]["concordance_1pct_all_interfaces"] else "sin concordancia segun el criterio",
    }
    (OUTPUT_DIR / "t6_flux_comparison.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")

    # Las cifras que sustituyen a la concordancia de 0.9 % se anaden a la
    # fuente unica de numeros (T5). Se reemplazan si ya estaban.
    numbers_path = REPO_ROOT / "build" / "source_a" / "numbers_r1.json"
    numbers = json.loads(numbers_path.read_text(encoding="utf-8"))
    numbers["numbers"] = [n for n in numbers["numbers"] if not n["key"].startswith("flux_comparison_")]
    p = summary["primaria"]
    src = "build/r1/t6_flux_comparison.json"
    fn = "code/analysis/r1_flux_comparison.py::main"
    why_pct = "RMS sobre 576 registros; basta una cifra para mostrar el orden frente al criterio de 1 %"
    for key, value, units, sig, why in (
        ("flux_comparison_n_records", int(n_rec), "records (Tsols 2-5, TMN91)", None, "conteo"),
        ("flux_comparison_max_rms_rel_F_down_pct", p["max_rms_rel_F_down"] * 100.0, "%", 1, why_pct),
        ("flux_comparison_max_rms_rel_F_up_pct", p["max_rms_rel_F_up"] * 100.0, "%", 1, why_pct),
        ("flux_comparison_concordance_1pct", p["concordance_1pct_all_interfaces"], "bool", None, "booleano"),
        ("flux_comparison_layer0_heating_archive", p["heating_layer0_archive_K_per_sol_mean"], "K per Titan sol", 2,
         "media de 576 registros; 2 cifras como las tasas de T5"),
        ("flux_comparison_layer0_heating_titan", p["heating_layer0_titan_K_per_sol_mean"], "K per Titan sol", 2, "idem"),
        ("flux_comparison_layer0_heating_rms_rel_pct", p["heating_layer0_rms_rel"] * 100.0, "%", 1, why_pct),
        ("flux_comparison_implied_archive_I_IR", p["implied_archive_I_IR_W_m2_mean"], "W m-2", 3,
         "comprobacion de geometria; 3 cifras bastan para mostrar 2.00"),
        ("flux_comparison_kirchhoff_surface_F_up_diff", summary["iii_kirchhoff"]["surface_mean_diff_F_up_W_m2"], "W m-2", 2,
         "variante diagnostica (iii)"),
    ):
        numbers["numbers"].append({
            "key": key, "value": value, "units": units, "lw_boundary": "full_absorption" if "kirchhoff" not in key else "kirchhoff",
            "sw_variant": "n/a (solo LW)", "source_file": src, "source_function": fn,
            "print_sig_figs": sig, "print_sig_figs_reason": why,
        })
    numbers_path.write_text(json.dumps(numbers, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=1, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
