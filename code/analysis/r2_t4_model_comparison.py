#!/usr/bin/env python
"""
r2_t4_model_comparison.py -- stage r2, analysis T4.

Diagnostic comparison table: TAM (Han & Lora 2025), Huygens cases 1-4,
against the mtWRF archive (Hartwick et al. 2026), baseline and suppressed
runs. No causal claims.

---------------------------------------------------------------------------
COLUMNAS Y DEFINICIONES (declaradas antes de ejecutar)
---------------------------------------------------------------------------
Por caso o corrida:
- amplitud diurna de Ts: pico a pico (max - min) dentro de cada dia, media
  de los dias; tambien la semiamplitud (= pico a pico / 2). mtWRF: TSK de los
  144 registros instantaneos de cada Tsol 2-5. TAM: tsurf de los registros
  cuyo centro cae en cada uno de los 30 dias centrales de T2 (registros que
  son medias de 1/24 de dia: la amplitud queda ligeramente subestimada).
- Ts - Ta, maximo y medio, con la altura de Ta: mtWRF TSK - T_PHY[0], Ta en
  el nivel de masa mas bajo; TAM tsurf - temp[nivel mas bajo], Ta en zg.
- H integrado por sol/dia: de T2 (HFX; sens), media de los dias.
- LW neto en superficie (+ hacia abajo): mtWRF XINET; TAM lwdn_sfc - lwup_sfc
  (diagnosticos archivados) y, aparte, lwdn - 0.95 sigma Ts^4 (lectura del
  balance de T1, descriptiva).
- termino R_O1: de T2.
- clase de frontera: T1 (TAM) y R0 (mtWRF: O1).
- altura del nivel mas bajo y numero de niveles por debajo de 4 km (alturas
  medias en el tiempo).
- PBL maxima por dia (descriptivo): mtWRF PBLH; TAM z_pbl.
Casos TAM descritos con sus atributos (los archivos solo traen 'filename' y
'title') y con la cita literal de Han & Lora (2025), Tabla 1 y p. 4.

Salida: build/r2/t4_model_comparison.csv (+ .json con las citas).
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

GRAVITY = 1.352
SOL_INDICES = {2: (87, 231), 3: (231, 375), 4: (375, 519), 5: (519, 663)}
TAM_HUYGENS = ("case1_huygens", "case2_huygens", "case3_huygens", "case4_huygens")
CENTRAL = range(10, 40)
SIGMA_FMS = 5.6734e-8

# Citas literales de Han & Lora (2025, PSJ 6, 222), transcritas de la fuente (uso interno).
HAN_LORA_CASES = {
    "case1_huygens": "Table 1 (p. 3): \"Case 1 3 × 10 6 0.01 173\" (rho c_p [J m-3 K-1], k [W m-1 K-1], I [TIU]); "
                     "Fig. 1 (p. 4): \"Case 1 (I = 173 TIU)\"",
    "case2_huygens": "Table 1 (p. 3): \"Case 2 3 × 10 6 0.1 548\"; Fig. 1 (p. 4): \"Case 2 (I = 548 TIU)\"",
    "case3_huygens": "Table 1 (p. 3): \"Case 3 3 × 10 6 1.0 1732\"; Fig. 1 (p. 4): \"Case 3 (I = 1732 TIU)\"",
    "case4_huygens": "Table 1 (p. 3): \"Case 4 f A spatially varying thermal inertia map is adopted.\"; p. 4: "
                     "\"In Case 4, the plains terrain type (I = 601 TIU) is assumed at the Huygens landing site\"",
}
HAN_LORA_COMMON = ("p. 3: \"A shortwave albedo of 0.2 and surface thermal emissivity of 0.95 are adopted.\"; "
                   "p. 3: \"50 atmospheric vertical levels and 15 subsurface vertical levels, using a 600 s time step\"")


def archive_rows(t2: dict) -> list[dict]:
    rows = []
    labels = {"TMN91": "mtWRF baseline TMN91 (principal, H-23)", "baseline5sol": "mtWRF baseline_5sol (alternativa)",
              "0.1flhc": "mtWRF 0.1flhc (flujo sensible suprimido)"}
    for run, path in C.ARCHIVE_RUNS.items():
        ds = netCDF4.Dataset(path)
        tsk = np.asarray(ds["TSK"][:, 0, 0], dtype=np.float64)
        ta = np.asarray(ds["T_PHY"][:, 0, 0, 0], dtype=np.float64)
        xinet = np.asarray(ds["XINET"][:, 0, 0], dtype=np.float64)
        pblh = np.asarray(ds["PBLH"][:, 0, 0], dtype=np.float64)
        z = np.asarray((ds["PH"][:, :, 0, 0] + ds["PHB"][:, :, 0, 0]) / GRAVITY, dtype=np.float64)
        ds.close()
        sl = slice(SOL_INDICES[2][0], SOL_INDICES[5][1])
        zm = (0.5 * (z[:, :-1] + z[:, 1:]) - z[:, :1])[sl].mean(axis=0)
        p2p = [float(np.ptp(tsk[i:j])) for i, j in SOL_INDICES.values()]
        contrast = tsk[sl] - ta[sl]
        ps = t2["archive"][run]["per_sol"]
        h = np.mean([ps[str(s)]["HFX"] if str(s) in ps else ps[s]["HFX"] for s in SOL_INDICES])
        r = np.mean([ps[str(s)]["R_O1"] if str(s) in ps else ps[s]["R_O1"] for s in SOL_INDICES])
        rows.append({
            "modelo": "mtWRF (Hartwick et al. 2026)", "caso": labels[run],
            "periodo": "Tsols 2-5",
            "Ts_amplitud_pico_a_pico_K": float(np.mean(p2p)), "Ts_semiamplitud_K": float(np.mean(p2p) / 2.0),
            "Ts_menos_Ta_max_K": float(contrast.max()), "Ts_menos_Ta_medio_K": float(contrast.mean()),
            "altura_Ta_m": float(zm[0]),
            "H_integrado_MJ_m2_por_dia": float(h),
            "LW_neto_sfc_W_m2_archivado": float(xinet[sl].mean()),
            "LW_neto_sfc_W_m2_lectura_balance": "",
            "R_O1_MJ_m2_por_dia": float(r), "R_O1_estado": "(1-0.95) GLW; clase O1 (R0)",
            "clase_frontera": "O1 (R0/T3: FLUX_LW_U[:,0] = 0.95 sigma TSK^4; XINET = GLW - 0.95 sigma TSK^4)",
            "altura_nivel_mas_bajo_m": float(zm[0]),
            "niveles_bajo_4km": int(np.sum(zm <= 4000.0)),
            "PBL_max_por_dia_m": float(np.mean([pblh[i:j].max() for i, j in SOL_INDICES.values()])),
            "descripcion": "corrida seca 1D del archivo Zenodo 17246644",
        })
    return rows


def tam_rows(t1: dict, t2: dict) -> list[dict]:
    rows = []
    period = t1["day_length_summary"]["method_A_mean_days"]
    for case in TAM_HUYGENS:
        ds = netCDF4.Dataset(C.tam_path(case))
        g = {v: np.asarray(ds[v][:], dtype=np.float64) for v in ("tsurf", "temp", "zg", "lwup_sfc", "lwdn_sfc", "z_pbl", "time")}
        attrs = {k: ds.getncattr(k) for k in ds.ncattrs()}
        ds.close()
        t = g["time"]
        dt = float(np.median(np.diff(t)))
        t0 = t[0] - dt / 2.0
        day = np.floor((t - t0) / period).astype(int)
        p2p, pbl = [], []
        for k in CENTRAL:
            m = day == k
            p2p.append(float(np.ptp(g["tsurf"][m])))
            pbl.append(float(g["z_pbl"][m].max()))
        mk = np.isin(day, list(CENTRAL))
        contrast = g["tsurf"][mk] - g["temp"][mk, -1]
        zmean = g["zg"][mk].mean(axis=0)
        f1 = t1["files"][case]
        s2 = t2["tam"][case]["30 centrales (k=10..39)"]
        r_o1 = s2["R_O1_equiv"]
        rows.append({
            "modelo": "TAM (Han & Lora 2025)", "caso": case,
            "periodo": "30 dias centrales de 50 (k = 10..39)",
            "Ts_amplitud_pico_a_pico_K": float(np.mean(p2p)), "Ts_semiamplitud_K": float(np.mean(p2p) / 2.0),
            "Ts_menos_Ta_max_K": float(contrast.max()), "Ts_menos_Ta_medio_K": float(contrast.mean()),
            "altura_Ta_m": float(zmean[-1]),
            "H_integrado_MJ_m2_por_dia": s2["sens"]["mean"],
            "LW_neto_sfc_W_m2_archivado": float(np.mean(g["lwdn_sfc"][mk] - g["lwup_sfc"][mk])),
            "LW_neto_sfc_W_m2_lectura_balance": float(np.mean(g["lwdn_sfc"][mk] - 0.95 * SIGMA_FMS * g["tsurf"][mk] ** 4)),
            "R_O1_MJ_m2_por_dia": "No determinable" if r_o1 is None else r_o1["mean"],
            "R_O1_estado": (t2["tam"][case]["R_O1_status"]
                            + f"; condicional eps=0.95 (balance): {s2['R_O1_equiv_if_O_eps095_balance']['mean']:.3f}"),
            "clase_frontera": (f"T1: {f1['class_all_sigmas']} (criterio sobre lwup_sfc; eps_O = "
                               f"{f1['per_sigma']['FMS_5.6734e-8']['median_eps_O']:.4f}); balance: O1 con eps = "
                               f"{f1['balance_test_explicit_emission']['FMS_5.6734e-8']['eps_best_fit_full_lwdn_absorption']:.4f} (descriptivo)"),
            "altura_nivel_mas_bajo_m": float(zmean[-1]),
            "niveles_bajo_4km": int(np.sum(zmean <= 4000.0)),
            "PBL_max_por_dia_m": float(np.mean(pbl)),
            "descripcion": f"atributos: {attrs}; Han & Lora: {HAN_LORA_CASES[case]}",
        })
    return rows


def main() -> int:
    t1 = json.loads((C.R2_OUT / "t1_tam_boundary.json").read_text(encoding="utf-8"))
    t2 = json.loads((C.R2_OUT / "t2_energy_partition.json").read_text(encoding="utf-8"))
    rows = tam_rows(t1, t2) + archive_rows(t2)
    with (C.R2_OUT / "t4_model_comparison.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    C.dump_json(C.R2_OUT / "t4_model_comparison.json", {
        "task": "R2-T4", "guard": C.guard_info(), "rows": rows,
        "han_lora_quotes": {"cases": HAN_LORA_CASES, "common": HAN_LORA_COMMON,
                            "source": "Han and Lora (2025), Planet. Sci. J. 6, 222 (transcription)"},
        "note": "solo diagnostico, sin afirmaciones causales",
    })
    for r in rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k not in ("descripcion",)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
