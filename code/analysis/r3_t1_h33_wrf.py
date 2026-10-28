#!/usr/bin/env python
"""
r3_t1_h33_wrf.py -- stage r3, analysis T1.

Verification of H-33 against the public WRF architecture and the archive.

6.1  Base version: literal quote of Hartwick (R2/T7) plus the TITLE attribute
     of the wrfout files.
6.2  Public WRF source (read-only clone outside the repo, path given with
     --wrf-files): literal transcription with file:line of (a) FLHC and HFX
     in MYJSFC and in the slab/thermal-diffusion LSM (SF_SURFACE_PHYSICS=1),
     (b) the lower boundary of MYJ PBL (BL_PBL_PHYSICS=2), (c) AKHS and THZ0
     in MYJSFC (SF_SFCLAY_PHYSICS=2).
6.3  Quantitative test in the archive: r = HFX / F_MYJ for daytime records.

CRITERIOS (declarados antes de ejecutar):
- 6.2 "Consistente": en la arquitectura revisada, la frontera inferior del
  MYJ no pasa por FLHC/HFX, de modo que escalar FLHC no altera el flujo que
  entra a la atmosfera. "Inconsistente": el MYJ usa HFX o FLHC. "No
  determinable": no se localiza el acoplamiento, o difiere entre versiones.
- 6.3 "Escalado confirmado en la salida": mediana(r) = 1.00 +/- 0.02 en las
  dos baselines y mediana(r) = 0.100 +/- 0.005 en 0.1flhc, con RIC <= 0.02 en
  cada corrida. "Parcial": se cumple el patron pero no las tolerancias. "No
  determinable": no hay variables suficientes.
- Filtro declarado: registros diurnos con |F_MYJ| > 0.05 W m-2.

DEFINICIONES (declaradas antes de ejecutar):
- Registro diurno: GSW > 0 (onda corta neta en superficie del propio archivo).
- F_MYJ (formula exacta del MYJ, VDIFH/frontera inferior): F = CP * rho_1 *
  AKHS * (THZ0 - theta_1), con THZ0 = THSK = TSK*(1e5/PSFC)**CAPA en puntos de
  tierra (MYJPBL sobrescribe THZ0 con THSK si SEAMASK = XLAND-1 = 0),
  theta_1 = T[0] + T0 (theta del propio modelo), rho_1 = p_1/(R_D*T_PHY[0])
  (seco: QVAPOR ausente o nulo). THZ0 no esta archivado: se reconstruye con la
  regla del codigo y se marca como tal. La tierra se infiere de LU_INDEX !=
  ISWATER (XLAND no esta archivado).
- Constantes del codigo: P0 = 1e5 Pa (literal en MYJPBL), CAPA = R_D/CP. CP =
  1101 (Hartwick, Tabla 2). T0, CAPA y P0 se obtienen del archivo por ajuste de
  theta = T + T0 = T_PHY*(P0/p)**CAPA (registro de la convencion del modelo);
  R_D = CAPA*CP.
- En paralelo se reporta la definicion de R2 (H_AKHS con R = 296.8 y theta
  referida a PSFC).
- Integrales por sol (Tsols 2-5, tres corridas): HFX, F_MYJ y residuo de
  energia de la columna (metodo de R2/T3b, columna completa y 0-4 km) con cada
  flujo.

Salidas: build/r3/t1_wrf_coupling.md, t1_h33_ratio.csv, t1_h33_summary.json,
quicklook/t1_h33_inspection.png.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import r3_common as C  # noqa: E402  (guard de R1)

import netCDF4  # noqa: E402
import numpy as np  # noqa: E402

TAGS = ["V3.6", "V3.0", "V3.9.1.1", "v4.8.0"]
TAG_ROLE = {"V3.6": "declarada por el atributo TITLE del archivo",
            "V3.0": "v3 mas antiguo disponible",
            "V3.9.1.1": "v3 tardio (ultimo v3)",
            "v4.8.0": "v4 mas reciente"}
MJ = 1.0e6

# Patrones de transcripcion: (item, archivo, regex). Las lineas se copian literalmente.
PATTERNS = [
    ("a", "module_sf_myjsfc.F", r"^\s*RLOW=PLOW/\(R_D\*TLOW\)"),
    ("a", "module_sf_myjsfc.F", r"^\s*HSFLX=AKHS\*\(THLOW-THZ0\)"),
    ("a", "module_sf_myjsfc.F", r"^\s*HFX=-RLOW\*CP\*HSFLX"),
    ("a", "module_sf_myjsfc.F", r"^\s*FLHC=RLOW\*CP\*AKHS"),
    ("a", "module_surface_driver.F", r"CALL SLAB\(t_phy,qv_curr,p_phy,flhc,flqc"),
    ("a", "module_sf_slab.F", r"^\s*THCON=\("),
    ("a", "module_sf_slab.F", r"^\s*THX\(I\)=SCR3\(I\)\*THCON"),
    ("a", "module_sf_slab.F", r"^\s*THG=TSK\(I\)\*\(100\./PSFC\(I\)\)\*\*ROVCP"),
    ("a", "module_sf_slab.F", r"^\s*HFX\(I\)=FLHC\(I\)\*\(THG-THX\(I\)\)"),
    ("a", "module_sf_slab.F", r"^\s*HFXT=FLHC\(I\)\*\(THG-THX\(I\)\)"),
    ("a", "module_sf_slab.F", r"^\s*QS\(I\)=HFX\(I\)\+QFX\(I\)\*XLV"),
    ("a", "module_sf_slab.F", r"^\s*DTHGDT\(I\)=\(RNET\(I\)-QS\(I\)\)/CAPG\(I\)"),
    ("b", "module_bl_myjpbl.F", r"REAL,PARAMETER :: CAPA=R_D/CP"),
    ("b", "module_bl_myjpbl.F", r"^\s*THSK\(I,J\)=TSK\(I,J\)\*\(1\.E5/PSFC\)\*\*CAPA"),
    ("b", "module_bl_myjpbl.F", r"^\s*RHOK\(K,I\)=PMID\(I,KFLIP,J\)/\(R_D\*T\(I,KFLIP,J\)"),
    ("b", "module_bl_myjpbl.F", r"^\s*SEAMASK=XLAND\(I,J\)-1\."),
    ("b", "module_bl_myjpbl.F", r"^\s*THZ0\(I,J\)=\(1\.-SEAMASK\)\*THSK\(I,J\)\+SEAMASK\*THZ0\(I,J\)"),
    ("b", "module_bl_myjpbl.F", r"^\s*AKHS_DENS=AKHS\(I,J\)\*RHOK\(KTE\+1-LLOW,I\)"),
    ("b", "module_bl_myjpbl.F", r"^\s*SZ0\(1\)=THZ0\(I,J\)"),
    ("b", "module_bl_myjpbl.F", r"^\s*CALL VDIFH\("),
    ("b", "module_bl_myjpbl.F", r"^\s*RKSS=RKHS\*CLOW\(M\)"),
    ("b", "module_bl_myjpbl.F", r"^\s*SPECIES\(M,LMH\)=\(DTOZS\*RKSS\*SZ0\(M\)\+RSSB\)/CMSB"),
    ("b", "module_bl_myjpbl.F", r"^\s*THE\(LMH\)=\(DTOZS\*RKHS\*THZ0\+RSTB\)/CMTB"),
    ("b", "module_pbl_driver.F", r"CALL myjpbl\("),
    ("b", "module_pbl_driver.F", r"TSK=tsk,QSFC=qsfc,CHKLOWQ=chklowq,THZ0=thz0"),
    ("b", "module_pbl_driver.F", r",AKHS=akhs,AKMS=akms,ELFLX=lh"),
    ("c", "module_sf_myjsfc.F", r"^\s*THSK\(I,J\)=TSK\(I,J\)/\(PSFC"),
    ("c", "module_sf_myjsfc.F", r"^\s*THZ0=THS\s*$"),
    ("c", "module_sf_myjsfc.F", r"^\s*AKHS=MAX\(USTARK/SIMH,CXCHL\)"),
    ("c", "module_sf_myjsfc.F", r"^\s*AKHS=MAX\(USTARK/SIMH,CXCHS\)"),
    ("c", "module_sf_myjsfc.F", r"^\s*CHS=AKHS"),
]

ITEM_TITLE = {
    "a": "(a) FLHC en la capa superficial MYJ y HFX a partir de FLHC en el esquema de superficie 1 (slab / difusion termica)",
    "b": "(b) Frontera inferior del MYJ PBL (BL_PBL_PHYSICS = 2)",
    "c": "(c) AKHS y THZ0 en la capa superficial MYJ (SF_SFCLAY_PHYSICS = 2), puntos de tierra",
}


def tag_commit(repo: Path, tag: str) -> str:
    try:
        return subprocess.run(["git", "-C", str(repo), "rev-parse", f"{tag}^{{}}"], capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:
        return "desconocido"


def myjpbl_call_block(lines: list[str]) -> tuple[int, list[str]]:
    """Return (first line number, lines) of the ``CALL myjpbl(`` block inside CASE (MYJPBLSCHEME)."""
    in_case = False
    for i, line in enumerate(lines):
        if "CASE (MYJPBLSCHEME)" in line:
            in_case = True
        if in_case and "CALL myjpbl(" in line:
            block = []
            for j in range(i, min(i + 40, len(lines))):
                block.append(lines[j])
                if re.match(r"^\s*\)\s*$", lines[j]) or re.match(r"^\s*&?\s*\)\s*$", lines[j]):
                    break
            return i + 1, block
    return -1, []


def analyse_wrf(wrf_files: Path, wrf_repo: Path | None) -> dict:
    result = {}
    for tag in TAGS:
        d = wrf_files / tag
        files = {p.name: p.read_text(encoding="latin-1").splitlines() for p in d.glob("*.F")}
        hits = []
        for item, fname, rx in PATTERNS:
            for n, line in enumerate(files.get(fname, []), start=1):
                if re.search(rx, line):
                    hits.append({"item": item, "file": f"phys/{fname}", "line": n, "text": line.rstrip()})
        pbl = files["module_bl_myjpbl.F"]
        refs_hfx = [(n, l.rstrip()) for n, l in enumerate(pbl, 1)
                    if re.search(r"\b(HFX|FLHC)\b", l, re.I) and not l.lstrip().startswith("!")]
        call_line, call_block = myjpbl_call_block(files["module_pbl_driver.F"])
        call_has_hfx = any(re.search(r"\b(hfx|flhc)\b", l, re.I) for l in call_block)
        call_has_thz0 = any("THZ0=thz0" in l for l in call_block)
        call_has_akhs = any("AKHS=akhs" in l for l in call_block)
        bottom_bc = [h for h in hits if h["item"] == "b" and ("SPECIES(M,LMH)" in h["text"] or "THE(LMH)" in h["text"])]
        thz0_land = [h for h in hits if h["item"] == "b" and "THZ0(I,J)=(1.-SEAMASK)" in h["text"]]
        located = bool(bottom_bc) and bool(thz0_land) and call_line > 0
        consistent = located and not refs_hfx and not call_has_hfx and call_has_thz0 and call_has_akhs
        result[tag] = {
            "role": TAG_ROLE[tag],
            "commit": tag_commit(wrf_repo, tag) if wrf_repo else "desconocido",
            "hits": hits,
            "myjpbl_uncommented_refs_to_HFX_or_FLHC": refs_hfx,
            "pbl_driver_call_line": call_line,
            "pbl_driver_call_block": [l.rstrip() for l in call_block],
            "call_passes_hfx_or_flhc": call_has_hfx,
            "call_passes_THZ0_and_AKHS": call_has_thz0 and call_has_akhs,
            "coupling_located": located,
            "verdict": ("Consistente" if consistent else ("Inconsistente" if (refs_hfx or call_has_hfx) else
                                                          "No determinable")),
        }
    verdicts = {r["verdict"] for r in result.values()}
    overall = verdicts.pop() if len(verdicts) == 1 else "No determinable"
    return {"per_version": result, "overall": overall}


def fit_theta_convention(ds) -> dict:
    """theta = T + T0 = T_PHY*(P0/p)**kappa over Tsols 2-5 (model convention)."""
    T = np.asarray(ds["T"][C.CORE, :, 0, 0], float)
    TP = np.asarray(ds["T_PHY"][C.CORE, :, 0, 0], float)
    P = np.asarray(ds["P"][C.CORE, :, 0, 0] + ds["PB"][C.CORE, :, 0, 0], float)
    x = np.log(P).ravel()
    A = np.vstack([np.ones_like(x), -x]).T
    best = None
    for t0 in np.arange(50.0, 400.0, 0.25):
        y = np.log((T + t0) / TP).ravel()
        coef = np.linalg.lstsq(A, y, rcond=None)[0]
        rms = float(np.sqrt(np.mean((A @ coef - y) ** 2)))
        if best is None or rms < best[0]:
            best = (rms, t0, coef)
    rms, t0, coef = best
    kappa = float(coef[1])
    return {"T0_K": float(t0), "kappa": kappa, "P0_Pa": float(np.exp(coef[0] / kappa)), "rms_log": rms}


def load_run(path: Path) -> dict:
    ds = netCDF4.Dataset(path)
    g2 = lambda v: np.asarray(ds[v][C.CORE, :, 0, 0], dtype=np.float64)  # noqa: E731
    g1 = lambda v: np.asarray(ds[v][C.CORE, 0, 0], dtype=np.float64)  # noqa: E731
    d = {
        "T": g2("T"), "T_PHY": g2("T_PHY"), "P": g2("P") + g2("PB"), "Z": (g2("PH") + g2("PHB")) / C.GRAVITY,
        "PSFC": g1("PSFC"), "TSK": g1("TSK"), "HFX": g1("HFX"), "AKHS": g1("AKHS"), "GSW": g1("GSW"),
        "LU_INDEX": g1("LU_INDEX"),
        "LW": g2("FLUX_LW_D") - g2("FLUX_LW_U"), "SW": g2("FLUX_SW_D") - g2("FLUX_SW_U"),
        "FLHC": g1("FLHC") if "FLHC" in ds.variables else None,
        "QVAPOR_max": float(np.max(np.abs(ds["QVAPOR"][:]))) if "QVAPOR" in ds.variables else None,
        "vars": {v: (v in ds.variables) for v in ("AKHS", "FLHC", "THZ0", "TSK", "T_PHY", "PSFC", "RHO",
                                                  "ALT", "QVAPOR", "XLAND", "LU_INDEX")},
        "attrs": {a: getattr(ds, a) for a in ("TITLE", "BL_PBL_PHYSICS", "SF_SFCLAY_PHYSICS",
                                              "SF_SURFACE_PHYSICS", "ISFFLX", "ISWATER")},
        "theta_fit": fit_theta_convention(ds),
    }
    ds.close()
    return d


WRF_TRANSCRIPTION = C.INPUTS / "wrf_myj_transcription.json"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--wrf-files", default=os.environ.get("WRF_SRC"),
                    help="dir with <tag>/module_*.F extracted from the WRF clone (default: $WRF_SRC)")
    ap.add_argument("--wrf-repo", default=None, help="WRF clone (to record tag commits)")
    args = ap.parse_args()
    if args.wrf_files and Path(args.wrf_files).is_dir():
        wrf = analyse_wrf(Path(args.wrf_files), Path(args.wrf_repo) if args.wrf_repo else None)
        wrf["source"] = "recomputed from the WRF source files"
    else:
        # Sin el codigo fuente de WRF: transcripcion registrada (archivo, linea y texto de cada coincidencia).
        wrf = json.loads(WRF_TRANSCRIPTION.read_text(encoding="utf-8"))
        wrf["source"] = "recorded transcription (inputs/wrf_myj_transcription.json); set WRF_SRC to recompute"

    rec = C.C2.load_r1_flux_module()
    dt = C.SOL_S / 144.0
    runs = {}
    for run, path in C.ARCHIVE_RUNS.items():
        runs[run] = load_run(path)

    # Constantes del modelo a partir del archivo (deben coincidir entre corridas).
    fits = {r: d["theta_fit"] for r, d in runs.items()}
    T0 = float(np.round(np.mean([f["T0_K"] for f in fits.values()]), 2))
    KAPPA = float(np.mean([f["kappa"] for f in fits.values()]))
    P0 = 1.0e5
    CP = C.CP_HARTWICK
    R_D = KAPPA * CP

    rows, summary, sol_rows = [], {}, []
    for run, d in runs.items():
        land = bool(np.all(d["LU_INDEX"] != d["attrs"]["ISWATER"]))
        p1, T1 = d["P"][:, 0], d["T_PHY"][:, 0]
        theta1 = d["T"][:, 0] + T0
        thz0 = d["TSK"] * (P0 / d["PSFC"]) ** KAPPA          # THSK = THZ0 en tierra (reconstruido)
        rho1 = p1 / (R_D * T1)
        f_myj = CP * rho1 * d["AKHS"] * (thz0 - theta1)
        # Definicion de R2 (paralela)
        kap_r2 = C.R_GAS_R1 / CP
        f_r2 = (p1 / (C.R_GAS_R1 * T1)) * CP * d["AKHS"] * (d["TSK"] - T1 * (d["PSFC"] / p1) ** kap_r2)
        hfx = d["HFX"]
        day = d["GSW"] > 0.0
        sel = day & (np.abs(f_myj) > 0.05)
        sel_r2 = day & (np.abs(f_r2) > 0.05)
        r = np.where(sel, hfx / np.where(f_myj == 0, np.nan, f_myj), np.nan)
        r2 = np.where(sel_r2, hfx / np.where(f_r2 == 0, np.nan, f_r2), np.nan)
        core = slice(1, -1)
        for k in range(1, len(hfx) - 1):
            rows.append({"run": run, "record": k + C.CORE.start, "GSW_W_m2": d["GSW"][k], "daytime": bool(day[k]),
                         "TSK_K": d["TSK"][k], "theta1_K": theta1[k], "THZ0_reconstructed_K": thz0[k],
                         "rho1_kg_m3": rho1[k], "AKHS_m_s": d["AKHS"][k], "HFX_W_m2": hfx[k],
                         "F_MYJ_W_m2": f_myj[k], "F_R2_def_W_m2": f_r2[k],
                         "r_HFX_over_F_MYJ": r[k] if sel[k] else None,
                         "r_HFX_over_F_R2def": r2[k] if sel_r2[k] else None,
                         "FLHC": d["FLHC"][k] if d["FLHC"] is not None else None})
        rr = r[core][np.isfinite(r[core])]
        rr2 = r2[core][np.isfinite(r2[core])]
        q25, q50, q75 = np.percentile(rr, [25, 50, 75])
        info = {
            "n_records_core": int(len(hfx) - 2), "n_daytime": int(day[core].sum()), "n_used": int(rr.size),
            "median_r": float(q50), "q25_r": float(q25), "q75_r": float(q75), "IQR_r": float(q75 - q25),
            "median_r_R2_definition": float(np.median(rr2)),
            "IQR_r_R2_definition": float(np.subtract(*np.percentile(rr2, [75, 25]))),
            "land_point": land, "QVAPOR_max": d["QVAPOR_max"], "variables_present": d["vars"],
            "attributes": d["attrs"], "theta_convention_fit": d["theta_fit"],
            "mean_HFX_daytime_W_m2": float(hfx[core][day[core]].mean()),
            "mean_F_MYJ_daytime_W_m2": float(f_myj[core][day[core]].mean()),
            "max_HFX_W_m2": float(hfx[core].max()), "max_F_MYJ_W_m2": float(f_myj[core].max()),
        }
        if d["FLHC"] is not None:
            ok = d["AKHS"][core] > 0
            info["median_FLHC_over_rho1_CP_AKHS"] = float(np.median(d["FLHC"][core][ok] /
                                                                    (rho1 * CP * d["AKHS"])[core][ok]))
            info["IQR_FLHC_over_rho1_CP_AKHS"] = float(np.subtract(*np.percentile(
                d["FLHC"][core][ok] / (rho1 * CP * d["AKHS"])[core][ok], [75, 25])))
        summary[run] = info

        # Energia de columna (metodo de R2/T3b; interfaces de R1/T6 con R = 296.8 por comparabilidad).
        pint = np.array([rec.reconstruct_interfaces(d["PSFC"][k], d["Z"][k], d["T_PHY"][k], C.GRAVITY, C.R_GAS_R1)
                         for k in range(d["PSFC"].size)])
        dp = pint[:, :-1] - pint[:, 1:]
        zm = (0.5 * (d["Z"][:, :-1] + d["Z"][:, 1:]) - d["Z"][:, :1]).mean(axis=0)
        dT = (d["T_PHY"][2:] - d["T_PHY"][:-2]) / (2.0 * dt)
        for label, ktop in (("0-4 km", int(np.argmin(np.abs(zm - 4000.0)))), ("full column (layers 0-54)", 55)):
            dE = (C.CP_HARTWICK * dp[core] / C.GRAVITY * dT)[:, :ktop].sum(axis=1)
            lw_c = (d["LW"][:, ktop] - d["LW"][:, 0])[core]
            sw_c = (d["SW"][:, ktop] - d["SW"][:, 0])[core]
            periods = {f"Tsol {s}": slice(a - 87, b - 87) for s, (a, b) in C.SOL_INDICES.items()}
            periods["Tsols 2-5"] = slice(0, 576)
            for pname, ps in periods.items():
                e = dE[ps].sum() * dt / MJ
                ih, im = hfx[core][ps].sum() * dt / MJ, f_myj[core][ps].sum() * dt / MJ
                ir2 = f_r2[core][ps].sum() * dt / MJ
                il, isw = lw_c[ps].sum() * dt / MJ, sw_c[ps].sum() * dt / MJ
                sol_rows.append({"run": run, "column": label, "period": pname,
                                 "dE_column_MJ_m2": e, "HFX_MJ_m2": ih, "F_MYJ_MJ_m2": im, "F_R2def_MJ_m2": ir2,
                                 "LW_convergence_MJ_m2": il, "SW_convergence_MJ_m2": isw,
                                 "residual_with_HFX_MJ_m2": e - (ih + il + isw),
                                 "residual_with_F_MYJ_MJ_m2": e - (im + il + isw)})

    # Veredicto 6.3
    base = [summary["TMN91"], summary["baseline5sol"]]
    sup = summary["0.1flhc"]
    tol_ok = (all(abs(b["median_r"] - 1.0) <= 0.02 for b in base) and abs(sup["median_r"] - 0.100) <= 0.005
              and all(s["IQR_r"] <= 0.02 for s in summary.values()))
    pattern = all(b["median_r"] > 0.5 for b in base) and sup["median_r"] < 0.2
    verdict63 = "Escalado confirmado en la salida" if tol_ok else ("Parcial" if pattern else "No determinable")

    C.write_csv(C.R3_OUT / "t1_h33_ratio.csv", rows)
    C.write_csv(C.R3_OUT / "t1_h33_sol_integrals.csv", sol_rows)
    out = {
        "task": "R3-T1", "guard": C.guard_info(), "head": C.git_head(),
        "version_6_1": {
            "hartwick_quote": ("The mtWRF model was developed in 2020 based on the Advance Research WRF (ARW) "
                               "dynamical core."),
            "hartwick_quote_page": "p. 3 (Sec. 2.1)",
            "hartwick_declares_version": False,
            "archive_TITLE": {r: d["attrs"]["TITLE"] for r, d in runs.items()},
            "versions_reviewed": TAGS,
        },
        "wrf_6_2": {"overall_verdict": wrf["overall"],
                    "per_version": {t: {k: v for k, v in r.items() if k != "hits"} for t, r in wrf["per_version"].items()}},
        "model_constants_from_archive": {"T0_K": T0, "kappa": KAPPA, "P0_Pa": P0, "CP": CP, "R_D": R_D,
                                         "fits": fits},
        "ratio_6_3": {"verdict": verdict63, "criterion": "baselines 1.00+/-0.02; 0.1flhc 0.100+/-0.005; IQR<=0.02",
                      "runs": summary},
        "sol_integrals": sol_rows,
    }
    C.dump_json(C.R3_OUT / "t1_h33_summary.json", out)
    write_markdown(wrf, out, runs)
    plot(rows, sol_rows)

    tags = {"TMN91": "TMN91", "baseline5sol": "baseline5sol", "0.1flhc": "0p1flhc"}
    src = C.rel(C.R3_OUT / "t1_h33_summary.json")
    fn = "code/analysis/r3_t1_h33_wrf.py::main"
    ents = [C.r3_entry("h33_wrf_coupling_verdict", wrf["overall"], "-", src, fn,
                       note="veredicto 6.2 (arquitectura WRF publica, 4 versiones)"),
            C.r3_entry("h33_ratio_verdict", verdict63, "-", src, fn, note="veredicto 6.3"),
            C.r3_entry("mtwrf_R_D_inferred", R_D, "J kg-1 K-1", src, fn,
                       note="R_D = kappa*CP, kappa ajustado de theta = T + T0 en el archivo"),
            C.r3_entry("mtwrf_theta_T0", T0, "K", src, fn, note="theta = T + T0 en el wrfout")]
    for run, s in summary.items():
        t = tags[run]
        ents.append(C.r3_entry(f"h33_median_r_{t}", s["median_r"], "-", src, fn,
                               uncertainty=s["IQR_r"], note="mediana de HFX/F_MYJ (diurno, |F|>0.05); uncertainty = RIC"))
        tot = [x for x in sol_rows if x["run"] == run and x["period"] == "Tsols 2-5"
               and x["column"].startswith("full")][0]
        ents.append(C.r3_entry(f"h33_colresid_HFX_{t}", tot["residual_with_HFX_MJ_m2"], "MJ m-2 (Tsols 2-5)", src, fn))
        ents.append(C.r3_entry(f"h33_colresid_FMYJ_{t}", tot["residual_with_F_MYJ_MJ_m2"], "MJ m-2 (Tsols 2-5)", src, fn))
        ents.append(C.r3_entry(f"h33_int_HFX_{t}", tot["HFX_MJ_m2"], "MJ m-2 (Tsols 2-5)", src, fn))
        ents.append(C.r3_entry(f"h33_int_FMYJ_{t}", tot["F_MYJ_MJ_m2"], "MJ m-2 (Tsols 2-5)", src, fn))
    C.stage_numbers(ents, ("h33_", "mtwrf_"), fn)

    print("6.2:", wrf["overall"], {t: r["verdict"] for t, r in wrf["per_version"].items()})
    print("constantes:", out["model_constants_from_archive"]["T0_K"], KAPPA, R_D)
    for run, s in summary.items():
        print(run, {k: (round(v, 4) if isinstance(v, float) else v) for k, v in s.items()
                    if k.startswith(("median", "IQR", "q", "n_", "mean", "max"))})
    print("6.3:", verdict63)
    for x in sol_rows:
        if x["period"] == "Tsols 2-5":
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in x.items()})
    return 0


def write_markdown(wrf: dict, out: dict, runs: dict) -> None:
    L = ["# R3/T1 — Acoplamiento superficie–MYJ en el código público de WRF (H-33)", "",
         "*Generado por `code/analysis/r3_t1_h33_wrf.py`. Transcripciones literales del clon público "
         "`https://github.com/wrf-model/WRF` (solo lectura, fuera del repo). Líneas comentadas excluidas de las "
         "comprobaciones automáticas.*", "",
         "## 6.1 Versión base", "",
         f"- Hartwick et al. (2026), {out['version_6_1']['hartwick_quote_page']}: "
         f"\"{out['version_6_1']['hartwick_quote']}\" → **versión no declarada**.",
         "- Atributo global `TITLE` de los tres wrfout: " +
         "; ".join(f"`{r}`: \"{t.strip()}\"" for r, t in out["version_6_1"]["archive_TITLE"].items()) +
         " → se usa **V3.6** como versión base, más V3.0, V3.9.1.1 y v4.8.0 (procedimiento 6.2).",
         "- Atributos de física del archivo: " +
         ", ".join(f"{k} = {v}" for k, v in runs["TMN91"]["attrs"].items() if k not in ("TITLE",)) + ".",
         "- Advertencia: mtWRF es una derivación modificada (radiación gris, subsuelo propio). El código de mtWRF "
         "no es público; lo que sigue verifica la **arquitectura estándar** de WRF, no el código de mtWRF.", ""]
    L.append("## 6.2 Transcripciones por versión")
    for tag, r in wrf["per_version"].items():
        L += ["", f"### WRF {tag} ({r['role']}; commit `{r['commit']}`)", ""]
        for item in ("a", "b", "c"):
            L += [f"**{ITEM_TITLE[item]}**", "", "```fortran"]
            for h in r["hits"]:
                if h["item"] == item:
                    L.append(f"{h['file']}:{h['line']}: {h['text'].strip()}")
            L += ["```", ""]
        L += [f"Llamada a `myjpbl` en `phys/module_pbl_driver.F:{r['pbl_driver_call_line']}` (bloque completo):", "",
              "```fortran"] + r["pbl_driver_call_block"] + ["```", "",
              f"- Referencias no comentadas a `HFX`/`FLHC` en `module_bl_myjpbl.F`: "
              f"{len(r['myjpbl_uncommented_refs_to_HFX_or_FLHC'])}",
              f"- La llamada pasa `hfx` o `flhc`: {'sí' if r['call_passes_hfx_or_flhc'] else 'no'}; "
              f"pasa `THZ0` y `AKHS`: {'sí' if r['call_passes_THZ0_and_AKHS'] else 'no'}",
              f"- **Veredicto {tag}: {r['verdict']}**"]
    mc = out["model_constants_from_archive"]
    L += ["", "## Lectura de la cadena (común a las cuatro versiones)", "",
          "1. La capa superficial MYJ calcula `AKHS` (similitud) y, en tierra, `THZ0 = THS` (θ de la piel); "
          "de ahí `FLHC = RLOW*CP*AKHS` y un `HFX` provisional.",
          "2. El esquema de superficie 1 recibe `FLHC`, recalcula `HFX = FLHC*(THG-THX)` y lo usa en el balance "
          "de la piel (`QS = HFX + QFX*XLV`; `DTHGDT = (RNET-QS)/CAPG ...`). Ese `HFX` es el que se archiva.",
          "3. El MYJ PBL **no recibe** `HFX` ni `FLHC`: impone su frontera inferior con "
          "`AKHS_DENS = AKHS*ρ₁` y `THZ0`, que en tierra sobrescribe con `THSK = TSK*(1e5/PSFC)**CAPA`, en la "
          "solución implícita de `VDIFH`. Flujo implícito: F = CP·ρ₁·AKHS·(THSK − θ₁).",
          "4. Por tanto, un factor aplicado a `FLHC` (o a `HFX`) **dentro del esquema de superficie** reduce el "
          "`HFX` archivado y el enfriamiento de la piel, pero no el flujo que el MYJ inyecta a la atmósfera.", "",
          f"**Veredicto 6.2 (criterio declarado): {wrf['overall']}.**", "",
          "## Constantes del modelo inferidas del archivo (para 6.3)", "",
          f"- θ = T + T0 con T0 = {mc['T0_K']:.2f} K; θ = T_PHY·(P0/p)^κ con P0 = {mc['P0_Pa']:.0f} Pa y "
          f"κ = {mc['kappa']:.6f} (ajuste, RMS en log ≤ "
          f"{max(f['rms_log'] for f in mc['fits'].values()):.1e}).",
          f"- Con CP = {mc['CP']:.1f} (Hartwick, Tabla 2): R_D = κ·CP = {mc['R_D']:.2f} J kg⁻¹ K⁻¹ "
          "(R1/R2 usaban 296.8).", ""]
    ratio = out["ratio_6_3"]
    L += ["## 6.3 Cociente r = HFX / F_MYJ (registros diurnos, |F_MYJ| > 0.05 W m⁻²)", "",
          "| Corrida | n | mediana(r) | RIC | mediana con la definición de R2 | FLHC/(ρ₁·CP·AKHS) |", "|---|---|---|---|---|---|"]
    for run, s in ratio["runs"].items():
        fl = s.get("median_FLHC_over_rho1_CP_AKHS")
        L.append(f"| {run} | {s['n_used']} | {s['median_r']:.4f} | {s['IQR_r']:.4f} | "
                 f"{s['median_r_R2_definition']:.4f} | {('%.4f' % fl) if fl is not None else 'no archivado'} |")
    L += ["", f"**Veredicto 6.3 (criterio declarado): {ratio['verdict']}.**", "",
          "Lectura (descriptiva): en `0.1flhc` el `FLHC` archivado es exactamente ρ₁·CP·AKHS, es decir, **no está "
          "escalado**; el factor 0.1 aparece solo en el `HFX` archivado. Con la cadena de 6.2, eso corresponde a "
          "un factor aplicado dentro del esquema de superficie (a la copia de `FLHC` que usa, o a su `HFX`); con "
          "las salidas no se puede distinguir entre ambas. El MYJ recibe `AKHS` sin escalar.", ""]
    (C.R3_OUT / "t1_wrf_coupling.md").write_text("\n".join(L) + "\n", encoding="utf-8")


def plot(rows: list[dict], sol_rows: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    C.QUICKLOOK.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    for run, col in (("TMN91", "C0"), ("baseline5sol", "C2"), ("0.1flhc", "C3")):
        rr = [r for r in rows if r["run"] == run]
        t = 2.0 + (np.array([r["record"] for r in rr]) - 87) / 144.0  # Tsol 2 empieza en el registro 87
        ax[0].plot(t, [r["HFX_W_m2"] for r in rr], color=col, lw=0.8, label=f"HFX {run}")
        ax[0].plot(t, [r["F_MYJ_W_m2"] for r in rr], color=col, lw=0.8, ls="--", label=f"F_MYJ {run}")
        ax[1].plot(t, [r["r_HFX_over_F_MYJ"] if r["r_HFX_over_F_MYJ"] is not None else np.nan for r in rr],
                   ".", ms=2, color=col, label=run)
    ax[0].set_ylabel("W m-2"); ax[0].legend(fontsize=7, ncol=3)
    ax[1].set_ylabel("r = HFX / F_MYJ"); ax[1].axhline(1, color="k", lw=0.5); ax[1].axhline(0.1, color="k", lw=0.5)
    ax[1].legend(fontsize=7)
    for run, col in (("TMN91", "C0"), ("baseline5sol", "C2"), ("0.1flhc", "C3")):
        s = [x for x in sol_rows if x["run"] == run and x["column"].startswith("full") and x["period"].startswith("Tsol ")]
        sols = [int(x["period"].split()[1]) for x in s]
        ax[2].plot(np.array(sols) + 1, np.cumsum([x["residual_with_HFX_MJ_m2"] for x in s]), "o-", color=col, label=f"{run} HFX")
        ax[2].plot(np.array(sols) + 1, np.cumsum([x["residual_with_F_MYJ_MJ_m2"] for x in s]), "s--", color=col, label=f"{run} F_MYJ")
    ax[2].set_ylabel("cumulative column residual [MJ m-2]"); ax[2].set_xlabel("Tsol (residual accumulated at the end of each sol)")
    ax[2].legend(fontsize=7, ncol=3)
    fig.tight_layout()
    fig.savefig(C.QUICKLOOK / "t1_h33_inspection.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
