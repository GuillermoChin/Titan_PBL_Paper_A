#!/usr/bin/env python
"""
r3_t3_flux_comparison_all.py -- stage r3, analysis T3.

Repeats the R1/T6 flux comparison on the
alternative baseline (1D_Huygens_baseline_5sol) and on the suppressed run
(1D_Huygens_0.1flhc_5sol), Tsols 2-5, with the primary configuration, the
same metrics and the same criterion. TMN91 is rerun as a reproduction
control (must match the values registered in R1/T6).

The R1 helpers ``reconstruct_interfaces`` and ``make_state`` are imported
unchanged from code/analysis/r1_flux_comparison.py.

CRITERIOS (copiados de R1/T6 y declarados antes de ejecutar):
- Configuracion primaria: preset huygens_baseline (tau_IR = 10, f = 0.15,
  I_IR = 2, eps = 0.95, full_absorption); malla con interfaces reconstruidas
  por hidrostatica (R = 296.8, g = 1.352), niveles de masa P + PB, T de capa
  T_PHY, Ts = TSK, p_s de la Ec. A13 = PSFC; registros 87-662.
- Metricas por interfaz k = 0..55: RMS relativo de F_up y F_down sobre los
  registros validos, maximo absoluto, diferencia media; tasa LW de la capa 0
  calculada igual para ambos.
- "Concordancia" si el RMS relativo de F_up y de F_down es <= 1 % en todas
  las interfaces.
- Control: TMN91 debe reproducir el maximo RMS de R1/T6 (diferencia relativa
  <= 1e-9); si no, se reporta.
- Diagnostico declarado (no sustituye a la primaria): 0.1flhc con interfaces
  de ZNW/P_TOP del propio archivo, p_int = ZNW (PSFC - P_TOP) + P_TOP
  (pregunta 10 de R2).

Salidas: build/r3/t3_flux_comparison_all.csv y .json, y
t3_flux_profiles.csv (perfiles medios para la figura de comparacion).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import r3_common as C  # noqa: E402

import netCDF4  # noqa: E402
import numpy as np  # noqa: E402

FIRST, LAST = 87, 663
CRITERION = 0.01


def compare(rec, config, ds, use_znw=False) -> tuple[dict, list[dict], list[dict]]:
    gravity = config.planet.gravity
    r_gas = config.atmosphere.gas.specific_gas_constant
    cp = config.atmosphere.gas.specific_heat_cp
    sol = config.planet.sol_length
    sl = slice(FIRST, LAST)
    D = np.asarray(ds["FLUX_LW_D"][sl, :, 0, 0], float)
    U = np.asarray(ds["FLUX_LW_U"][sl, :, 0, 0], float)
    P = np.asarray(ds["P"][sl, :, 0, 0] + ds["PB"][sl, :, 0, 0], float)
    Z = np.asarray((ds["PH"][sl, :, 0, 0] + ds["PHB"][sl, :, 0, 0]) / gravity, float)
    T = np.asarray(ds["T_PHY"][sl, :, 0, 0], float)
    TSK = np.asarray(ds["TSK"][sl, 0, 0], float)
    PSFC = np.asarray(ds["PSFC"][sl, 0, 0], float)
    n_rec, n_lev = D.shape
    valid = ~((D[:, 0] == 0.0) & (U[:, 0] == 0.0))
    if use_znw:
        znw = np.asarray(ds["ZNW"][sl, :], float)
        ptop = np.asarray(ds["P_TOP"][sl], float)
        p_int = znw * (PSFC - ptop)[:, None] + ptop[:, None]
    else:
        p_int = np.array([rec.reconstruct_interfaces(PSFC[r], Z[r], T[r], gravity, r_gas) for r in range(n_rec)])
    Fd = np.empty((n_rec, n_lev + 1))
    Fu = np.empty((n_rec, n_lev + 1))
    for r in range(n_rec):
        state, scheme = rec.make_state(config, p_int[r], P[r], T[r], TSK[r], PSFC[r])
        Fd[r], Fu[r] = scheme.longwave_fluxes(state)
    Fd, Fu, Da, Ua, pv = Fd[valid], Fu[valid], D[valid], U[valid], p_int[valid]
    rel_d = (Fd[:, :n_lev] - Da) / Da
    rel_u = (Fu[:, :n_lev] - Ua) / Ua
    rms_d = np.sqrt(np.mean(rel_d**2, axis=0))
    rms_u = np.sqrt(np.mean(rel_u**2, axis=0))
    rows, prof = [], []
    zint = (Z[valid] - Z[valid][:, :1]).mean(axis=0)
    for k in range(n_lev):
        rows.append({"interface": k, "rms_rel_F_down": rms_d[k], "rms_rel_F_up": rms_u[k],
                     "max_abs_diff_F_down_W_m2": float(np.max(np.abs(Fd[:, k] - Da[:, k]))),
                     "max_abs_diff_F_up_W_m2": float(np.max(np.abs(Fu[:, k] - Ua[:, k]))),
                     "mean_diff_F_down_W_m2": float(np.mean(Fd[:, k] - Da[:, k])),
                     "mean_diff_F_up_W_m2": float(np.mean(Fu[:, k] - Ua[:, k])),
                     "mean_F_down_archive_W_m2": float(np.mean(Da[:, k])),
                     "mean_F_up_archive_W_m2": float(np.mean(Ua[:, k]))})
        prof.append({"interface": k, "z_interface_m": float(zint[k]), "p_interface_Pa": float(pv[:, k].mean()),
                     "F_down_archive": float(Da[:, k].mean()), "F_down_titan_pbl": float(Fd[:, k].mean()),
                     "F_up_archive": float(Ua[:, k].mean()), "F_up_titan_pbl": float(Fu[:, k].mean()),
                     "mean_rel_diff_F_down": float(rel_d[:, k].mean()), "mean_rel_diff_F_up": float(rel_u[:, k].mean()),
                     "rms_rel_F_down": float(rms_d[k]), "rms_rel_F_up": float(rms_u[k])})
    dp = pv[:, :n_lev - 1] - pv[:, 1:n_lev]
    net_a, net_t = Da - Ua, Fd[:, :n_lev] - Fu[:, :n_lev]
    h_a = gravity * (net_a[:, 1:] - net_a[:, :-1]) / (cp * dp) * sol
    h_t = gravity * (net_t[:, 1:] - net_t[:, :-1]) / (cp * dp) * sol
    summary = {
        "n_records": int(n_rec), "n_valid": int(valid.sum()),
        "max_rms_rel_F_down": float(rms_d.max()), "interface_of_max_F_down": int(rms_d.argmax()),
        "max_rms_rel_F_up": float(rms_u.max()), "interface_of_max_F_up": int(rms_u.argmax()),
        "concordance_1pct_all_interfaces": bool(rms_d.max() <= CRITERION and rms_u.max() <= CRITERION),
        "n_interfaces_over_1pct_F_down": int(np.sum(rms_d > CRITERION)),
        "n_interfaces_over_1pct_F_up": int(np.sum(rms_u > CRITERION)),
        "surface_mean_diff_F_up_W_m2": float(np.mean(Fu[:, 0] - Ua[:, 0])),
        "surface_mean_diff_F_down_W_m2": float(np.mean(Fd[:, 0] - Da[:, 0])),
        "heating_layer0_archive_K_per_sol_mean": float(np.mean(h_a[:, 0])),
        "heating_layer0_titan_K_per_sol_mean": float(np.mean(h_t[:, 0])),
        "heating_layer0_rms_rel": float(np.sqrt(np.mean(((h_t[:, 0] - h_a[:, 0]) / h_a[:, 0]) ** 2))),
        "heating_all_layers_rms_abs_K_per_sol": float(np.sqrt(np.mean((h_t - h_a) ** 2))),
        "dp_layer0_mean_Pa": float(dp[:, 0].mean()),
    }
    return summary, rows, prof


def main() -> int:
    rec = C.C2.load_r1_flux_module()
    from titan_pbl.config import load_config

    config = load_config("huygens_baseline")
    r1 = __import__("json").loads((C.R1_OUT / "t6_flux_comparison.json").read_text(encoding="utf-8"))
    r1p = r1["summary"]["primaria"]
    rows_all, prof_all, out = [], [], {}
    plan = [("TMN91", False, "control de reproduccion de R1/T6"),
            ("baseline5sol", False, "primaria"), ("0.1flhc", False, "primaria"),
            ("0.1flhc", True, "diagnostico: interfaces ZNW/P_TOP")]
    for run, znw, role in plan:
        ds = netCDF4.Dataset(C.ARCHIVE_RUNS[run])
        summary, rows, prof = compare(rec, config, ds, use_znw=znw)
        ds.close()
        label = f"{run}{'_ZNW' if znw else ''}"
        summary["role"] = role
        out[label] = summary
        rows_all += [{"run": label, **r} for r in rows]
        prof_all += [{"run": label, **p} for p in prof]
        print(label, role, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in summary.items()
                            if k.startswith(("max_rms", "concordance", "n_interfaces", "heating_layer0", "surface"))})
    ctrl = out["TMN91"]
    reproduces = (abs(ctrl["max_rms_rel_F_down"] - r1p["max_rms_rel_F_down"]) <= 1e-9 * r1p["max_rms_rel_F_down"]
                  and abs(ctrl["max_rms_rel_F_up"] - r1p["max_rms_rel_F_up"]) <= 1e-9 * r1p["max_rms_rel_F_up"])
    verdicts = {k: ("concordancia" if v["concordance_1pct_all_interfaces"] else "sin concordancia segun el criterio")
                for k, v in out.items()}
    C.write_csv(C.R3_OUT / "t3_flux_comparison_all.csv", rows_all)
    C.write_csv(C.R3_OUT / "t3_flux_profiles.csv", prof_all)
    C.dump_json(C.R3_OUT / "t3_flux_comparison_all.json", {
        "task": "R3-T3", "guard": C.guard_info(), "head": C.git_head(),
        "criterion": "concordancia si RMS relativo de F_up y F_down <= 1 % en todas las interfaces",
        "primary_configuration": "igual a R1/T6", "r1_control_reproduced": reproduces,
        "r1_registered": {"max_rms_rel_F_down": r1p["max_rms_rel_F_down"], "max_rms_rel_F_up": r1p["max_rms_rel_F_up"]},
        "summary": out, "verdicts": verdicts,
    })
    src = C.rel(C.R3_OUT / "t3_flux_comparison_all.json")
    fn = "code/analysis/r3_t3_flux_comparison_all.py::main"
    ents = []
    for label, tag in (("baseline5sol", "baseline5sol"), ("0.1flhc", "0p1flhc")):
        s = out[label]
        ents += [C.r3_entry(f"flux_comparison_{tag}_max_rms_rel_F_down_pct", s["max_rms_rel_F_down"] * 100, "%", src, fn,
                            boundary="full_absorption", variant="n/a (solo LW)"),
                 C.r3_entry(f"flux_comparison_{tag}_max_rms_rel_F_up_pct", s["max_rms_rel_F_up"] * 100, "%", src, fn,
                            boundary="full_absorption", variant="n/a (solo LW)"),
                 C.r3_entry(f"flux_comparison_{tag}_concordance_1pct", s["concordance_1pct_all_interfaces"], "bool", src, fn,
                            boundary="full_absorption", variant="n/a (solo LW)"),
                 C.r3_entry(f"flux_comparison_{tag}_n_records", s["n_valid"], "records (Tsols 2-5)", src, fn)]
    C.stage_numbers(ents, ("flux_comparison_baseline5sol", "flux_comparison_0p1flhc"), fn)
    print("control R1 reproducido:", reproduces, "| veredictos:", verdicts)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
