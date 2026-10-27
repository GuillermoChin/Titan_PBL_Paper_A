#!/usr/bin/env python
"""
r2_t3b_column_energy.py -- stage r2, analysis T3.

Column energy closure of each archived run and the surface heat flux implied
by the archived exchange coefficient. Added at closure (T10), after the T3
energy budget showed that the 0-4 km column of the suppressed run gains
~3 MJ m-2 per sol more than the baseline, far more than any difference in the
archived surface fluxes. Descriptive: it does not change the declared T3
verdicts.

METODO (declarado antes de ejecutar este script):
- Columna = capas 0..k_top con interfaces reconstruidas (R1/T6); contenido
  de entalpia E = sum cp T dp / g, cp = 1101 (Hartwick, Tabla 2).
- Cambio de E en los Tsols 2-5 = sum_t sum_k cp dp/g (dT/dt) dt, dT/dt por
  diferencia centrada (dt = 9570 s).
- Entradas archivadas: HFX en la base; convergencia LW y SW de la columna
  (F_net[k_top] - F_net[0], F_net = D - U, + hacia abajo).
- Residuo_HFX = cambio de E - (HFX + LW + SW).
- Flujo implicito por el coeficiente de intercambio: H_AKHS = rho cp AKHS
  (theta_s - theta_1), rho = p_1/(R T_1), theta_s = TSK, theta_1 = T_1
  (PSFC/p_1)^kappa, R = 296.8; AKHS = "SFC EXCH COEFF FOR HEAT" [m s-1].
  Residuo_AKHS = cambio de E - (H_AKHS + LW + SW).
- Criterio descriptivo (declarado): se considera que H_AKHS explica la
  energia de la columna si |Residuo_AKHS| <= 0.1 |Residuo_HFX| y
  |Residuo_AKHS| <= 0.1 * integral de H_AKHS.
- Columnas: 0-4 km (k_top = nivel mas cercano a 4000 m) y completa (capas
  0-54). Por sol y en los Tsols 2-5.

Salidas: build/r2/t3b_column_energy.csv y .json.
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

import r2_common as C  # noqa: E402

GRAVITY, R_GAS, CP = 1.352, 296.8, 1101.0
SOL_INDICES = {2: (87, 231), 3: (231, 375), 4: (375, 519), 5: (519, 663)}
MJ = 1.0e6


def main() -> int:
    rec = C.load_r1_flux_module()
    from titan_pbl.config import load_config

    sol = float(load_config("huygens_baseline").planet.sol_length)
    dt = sol / 144.0
    kappa = R_GAS / CP
    rows, out = [], {}
    for run, path in C.ARCHIVE_RUNS.items():
        ds = netCDF4.Dataset(path)
        sl = slice(86, 664)
        g2 = lambda v: np.asarray(ds[v][sl, :, 0, 0], dtype=np.float64)  # noqa: E731
        g1 = lambda v: np.asarray(ds[v][sl, 0, 0], dtype=np.float64)  # noqa: E731
        T, Z = g2("T_PHY"), (g2("PH") + g2("PHB")) / GRAVITY
        P = g2("P") + g2("PB")
        PS, TSK, HFX, AKHS = g1("PSFC"), g1("TSK"), g1("HFX"), g1("AKHS")
        lw = g2("FLUX_LW_D") - g2("FLUX_LW_U")
        sw = g2("FLUX_SW_D") - g2("FLUX_SW_U")
        flhc = g1("FLHC") if "FLHC" in ds.variables else None
        desc = {"AKHS": getattr(ds["AKHS"], "description", ""), "HFX": getattr(ds["HFX"], "description", "")}
        if flhc is not None:
            desc["FLHC"] = getattr(ds["FLHC"], "description", "")
        ds.close()
        pint = np.array([rec.reconstruct_interfaces(PS[r], Z[r], T[r], GRAVITY, R_GAS) for r in range(PS.size)])
        dp = pint[:, :-1] - pint[:, 1:]
        zm = (0.5 * (Z[:, :-1] + Z[:, 1:]) - Z[:, :1]).mean(axis=0)
        core = slice(1, -1)
        rho = P[:, 0] / (R_GAS * T[:, 0])
        dtheta = TSK - T[:, 0] * (PS / P[:, 0]) ** kappa
        h_akhs = (rho * CP * AKHS * dtheta)[core]
        hfx = HFX[core]
        dT = (T[2:] - T[:-2]) / (2.0 * dt)
        info = {"descriptions": desc,
                "median_HFX_over_H_AKHS_where_H_AKHS_gt_0.05": float(np.median((hfx / h_akhs)[h_akhs > 0.05])),
                "mean_HFX_W_m2": float(hfx.mean()), "mean_H_AKHS_W_m2": float(h_akhs.mean()),
                "max_HFX_W_m2": float(hfx.max()), "max_H_AKHS_W_m2": float(h_akhs.max())}
        if flhc is not None:
            info["median_FLHC_over_rho_cp_AKHS"] = float(np.median(flhc[core] / (rho * CP * AKHS)[core]))
        for label, ktop in (("0-4 km", int(np.argmin(np.abs(zm - 4000.0)))), ("columna completa (capas 0-54)", 55)):
            dE_t = (CP * dp[core] / GRAVITY * dT)[:, :ktop].sum(axis=1)
            lw_c = (lw[:, ktop] - lw[:, 0])[core]
            sw_c = (sw[:, ktop] - sw[:, 0])[core]
            periods = {f"Tsol {s}": slice(a - 87, b - 87) for s, (a, b) in SOL_INDICES.items()}
            periods["Tsols 2-5"] = slice(0, 576)
            for pname, ps in periods.items():
                e = dE_t[ps].sum() * dt / MJ
                ih, ia = hfx[ps].sum() * dt / MJ, h_akhs[ps].sum() * dt / MJ
                il, isw = lw_c[ps].sum() * dt / MJ, sw_c[ps].sum() * dt / MJ
                res_h, res_a = e - (ih + il + isw), e - (ia + il + isw)
                explains = bool(abs(res_a) <= 0.1 * abs(res_h) and abs(res_a) <= 0.1 * abs(ia))
                rows.append({"run": run, "column": label, "z_top_m": float(zm[ktop - 1]) if ktop < 55 else float(zm[54]),
                             "period": pname, "dE_column_MJ_m2": e, "HFX_MJ_m2": ih, "H_AKHS_MJ_m2": ia,
                             "LW_convergence_MJ_m2": il, "SW_convergence_MJ_m2": isw,
                             "residual_with_HFX_MJ_m2": res_h, "residual_with_H_AKHS_MJ_m2": res_a,
                             "H_AKHS_explains_column": explains})
        out[run] = info
    with (C.R2_OUT / "t3b_column_energy.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    C.dump_json(C.R2_OUT / "t3b_column_energy.json", {
        "task": "R2-T3 complemento (cierre)", "guard": C.guard_info(), "runs": out,
        "totals_Tsols_2_5": [r for r in rows if r["period"] == "Tsols 2-5"],
        "note": "descriptivo; anadido en el cierre tras ver el presupuesto de energia de T3",
    })
    fn = "code/analysis/r2_t3b_column_energy.py::main"
    src = C.rel(C.R2_OUT / "t3b_column_energy.json")
    nums = []
    for run in C.ARCHIVE_RUNS:
        tag = run.replace(".", "p")
        nums.append(C.entry(f"colenergy_HFX_over_H_AKHS_{tag}", out[run]["median_HFX_over_H_AKHS_where_H_AKHS_gt_0.05"], "-",
                            "n/a", "n/a", src, fn, 2, "mediana de cocientes de lecturas del archivo"))
        tot = [r for r in rows if r["run"] == run and r["period"] == "Tsols 2-5" and r["column"].startswith("columna")][0]
        nums.append(C.entry(f"colenergy_residual_HFX_full_column_{tag}", tot["residual_with_HFX_MJ_m2"], "MJ m-2 (Tsols 2-5)",
                            "n/a", "n/a", src, fn, 2, "cierre de columna con HFX"))
        nums.append(C.entry(f"colenergy_residual_H_AKHS_full_column_{tag}", tot["residual_with_H_AKHS_MJ_m2"], "MJ m-2 (Tsols 2-5)",
                            "n/a", "n/a", src, fn, 2, "cierre de columna con el flujo de AKHS"))
    C.update_numbers_r2(nums, "colenergy_", fn)
    for run, info in out.items():
        print(run, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in info.items() if k != "descriptions"})
    for r in rows:
        if r["period"] == "Tsols 2-5":
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
