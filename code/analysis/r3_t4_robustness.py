#!/usr/bin/env python
"""
r3_t4_robustness.py -- stage r3, analysis T4.

HASI uncertainty and clarification of the most favourable corner. Uses the
R2/T5 helpers (make_ctx, dT_req, mixing_time, OAT, PRESET) imported unchanged
from code/analysis/r2_t5_sw_dtreq.py, which in turn use
code/analysis/source_a_core.py unchanged.

CRITERIOS Y METODO (declarados antes de ejecutar):
1. Para las dos fronteras y la variante principal (solo LW), Delta T_req(4) y
   Delta T_req(5) con:
   - Delta_theta = 0.652 K (nominal);
   - Delta_theta - sigma y Delta_theta + sigma, sigma = sqrt(2) * 0.25 K
     (errores independientes, D20);
   - error sistematico (sin cambio respecto del nominal).
   Si Delta_theta - sigma <= 0: "sin estratificacion resoluble", tal cual.
2. La esquina: recalcular la esquina mas favorable por separado para las
   variantes solo LW y mediodia, y declarar cual corresponde al 1.76 K de R2.
   Regla de esquina (la de R2/T5.4): para cada factor del barrido OAT, el
   valor que minimiza Delta T_req(4); empate -> el primero de la lista. Se
   aplica con kirchhoff (como R2) y, como complemento declarado, con
   full_absorption.
3. Combinacion mas desfavorable para la conclusion: esquina de la variante
   principal con Delta_theta - sigma; Delta T_req(4) frente al contraste
   maximo verificado (0.69 K, archivo; clave
   archive_contrast_max_sols2to5_suppressed de numbers_r1).
Complemento declarado: la sigma de R2 (0.25 K * theta/T en cada extremo,
0.37 K) como sensibilidad, en filas aparte.

Salidas: build/r3/t4_robustness.csv y t4_robustness.json.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import r3_common as C  # noqa: E402

import numpy as np  # noqa: E402

SIGMA_T = 0.25
SIGMA_IND = float(np.sqrt(2.0) * SIGMA_T)
BOUNDARIES = ("full_absorption", "kirchhoff")
VARIANT = "lw"


def load_r2_t5():
    spec = importlib.util.spec_from_file_location("r2_t5_sw_dtreq", C.AUDIT / "r2_t5_sw_dtreq.py")
    m = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(m)
    return m


def overrides_for(boundary: str, values: dict, s0_preset: float) -> dict:
    ov = {"radiation.surface_lw_boundary": boundary}
    for k, val in values.items():
        if k == "S0_factor":
            ov["radiation.solar_constant"] = s0_preset * val
        else:
            ov[k] = val
    return ov


def corner(t5, mod, boundary: str, variant: str, s0_preset: float) -> dict:
    best, ties, sweep = {}, {}, []
    for factor, values in t5.OAT.items():
        scores = {}
        for val in values:
            vals = dict(t5.PRESET)
            vals[factor] = val
            ctx = t5.make_ctx(mod, overrides_for(boundary, vals, s0_preset))
            scores[val] = t5.dT_req(mod, ctx, 4, variant)
            sweep.append({"factor": factor, "value": val, "dT_req_4sol_K": scores[val]})
        finite = {k: s for k, s in scores.items() if np.isfinite(s)}
        best[factor] = min(finite, key=finite.get)
        ties[factor] = bool(np.ptp(list(finite.values())) == 0.0)
    ctx = t5.make_ctx(mod, overrides_for(boundary, best, s0_preset))
    thr = mod.compute_thresholds_variant(ctx, variant)
    return {"boundary": boundary, "sw_variant": variant, "chosen_values": best, "factor_without_effect_tie": ties,
            "heating_threshold_layer0_K": thr["heating"], "destabilisation_threshold_0_4km_K": thr["destabilisation"],
            "dT_req_4sol_K": t5.dT_req(mod, ctx, 4, variant), "dT_req_5sol_K": t5.dT_req(mod, ctx, 5, variant),
            "ctx_overrides": overrides_for(boundary, best, s0_preset), "sweep": sweep}


def main() -> int:
    t5 = load_r2_t5()
    mod = t5.C.load_figure_module()
    dtheta = float(mod.DELTA_THETA_4KM)
    r2_unc = __import__("json").loads((C.R2_OUT / "t5_hasi_uncertainty.json").read_text(encoding="utf-8"))
    sigma_r2 = float(r2_unc["independent"]["sigma_delta_theta_K"])
    r2_corner = __import__("json").loads((C.R2_OUT / "t5_corner.json").read_text(encoding="utf-8"))
    contrast_max = float(t5.C.number_r1("archive_contrast_max_sols2to5_suppressed"))
    base_cfg = t5.make_ctx(mod, {"radiation.surface_lw_boundary": "kirchhoff"})["config"]
    s0_preset = float(base_cfg.radiation.solar_constant)

    cases = [("nominal", dtheta, "Delta_theta nominal"),
             ("minus_sigma_independent", dtheta - SIGMA_IND, "Delta_theta - sqrt(2)*0.25 K"),
             ("plus_sigma_independent", dtheta + SIGMA_IND, "Delta_theta + sqrt(2)*0.25 K"),
             ("systematic", dtheta, "error sistematico: se cancela en Delta_theta"),
             ("minus_sigma_R2_sensitivity", dtheta - sigma_r2, "sensibilidad: sigma de R2 (theta/T)"),
             ("plus_sigma_R2_sensitivity", dtheta + sigma_r2, "sensibilidad: sigma de R2 (theta/T)")]

    def req_row(boundary, config_label, ctx, case, dth, note, variant=VARIANT):
        if dth <= 0:
            return {"lw_boundary": boundary, "sw_variant": variant, "configuration": config_label, "case": case,
                    "delta_theta_K": dth, "dT_req_4sol_K": None, "dT_req_5sol_K": None,
                    "status": "sin estratificacion resoluble", "note": note}
        r4 = t5.dT_req(mod, ctx, 4, variant, dth)
        r5 = t5.dT_req(mod, ctx, 5, variant, dth)
        return {"lw_boundary": boundary, "sw_variant": variant, "configuration": config_label, "case": case,
                "delta_theta_K": dth, "dT_req_4sol_K": r4, "dT_req_5sol_K": r5,
                "contrast_max_verified_K": contrast_max, "ratio_contrast_over_dT_req4": contrast_max / r4,
                "contrast_reaches_dT_req4": bool(contrast_max >= r4), "status": "ok", "note": note}

    rows = []
    for b in BOUNDARIES:
        ctx = t5.make_ctx(mod, {"radiation.surface_lw_boundary": b})
        for case, dth, note in cases:
            rows.append(req_row(b, "preset", ctx, case, dth, note))

    corners = {}
    for b in ("kirchhoff", "full_absorption"):
        for v in ("lw", "noon"):
            corners[f"{b}|{v}"] = corner(t5, mod, b, v, s0_preset)
            c = corners[f"{b}|{v}"]
            print("esquina", b, v, c["chosen_values"], round(c["dT_req_4sol_K"], 4), round(c["dT_req_5sol_K"], 4))
    r2_val = float(r2_corner["corner_results"]["dT_req_4sol_K"])
    match = {k: bool(abs(c["dT_req_4sol_K"] - r2_val) <= 1e-6) for k, c in corners.items()}
    corresponds = [k for k, m in match.items() if m]

    # Peor caso: esquina de la variante principal con Delta_theta - sigma (y la de mediodia como complemento).
    worst = {}
    for key in ("kirchhoff|lw", "full_absorption|lw", "kirchhoff|noon"):
        c = corners[key]
        b, v = key.split("|")
        ctx = t5.make_ctx(mod, c["ctx_overrides"])
        for case, dth, note in cases[:4]:
            rows.append(req_row(b, f"corner_{v}", ctx, case, dth, note, v))
        wr = req_row(b, f"corner_{v}", ctx, "minus_sigma_independent", dtheta - SIGMA_IND, "peor caso", v)
        worst[key] = wr
    C.write_csv(C.R3_OUT / "t4_robustness.csv", rows)
    out = {
        "task": "R3-T4", "guard": C.guard_info(), "head": C.git_head(),
        "delta_theta_nominal_K": dtheta, "sigma_independent_K": SIGMA_IND, "sigma_R2_K": sigma_r2,
        "contrast_max_verified_K": contrast_max,
        "corners": {k: {kk: vv for kk, vv in c.items() if kk not in ("sweep", "ctx_overrides")} for k, c in corners.items()},
        "corner_sweeps": {k: c["sweep"] for k, c in corners.items()},
        "r2_corner_dT_req4_K": r2_val, "r2_corner_matches": match,
        "r2_corner_corresponds_to": corresponds,
        "r2_corner_code_evidence": ("r2_t5_sw_dtreq.py::main, evaluate(best) usa variant=headline; headline = "
                                    f"'{r2_corner['sw_variant']}' en t5_corner.json"),
        "worst_case": worst,
        "rows": rows,
    }
    C.dump_json(C.R3_OUT / "t4_robustness.json", out)

    src = C.rel(C.R3_OUT / "t4_robustness.json")
    fn = "code/analysis/r3_t4_robustness.py::main"
    ents = [C.r3_entry("hasi_sigma_dtheta_independent_brief", SIGMA_IND, "K", src, fn, note="sqrt(2)*0.25 K (D20)")]
    for r in rows:
        if r["status"] != "ok" or r["case"] == "systematic":
            continue
        tag = f"{r['lw_boundary']}_{r['sw_variant']}_{r['configuration']}_{r['case']}"
        ents.append(C.r3_entry(f"robust_dTreq4_{tag}", r["dT_req_4sol_K"], "K", src, fn,
                               boundary=r["lw_boundary"], variant=r["sw_variant"]))
        ents.append(C.r3_entry(f"robust_dTreq5_{tag}", r["dT_req_5sol_K"], "K", src, fn,
                               boundary=r["lw_boundary"], variant=r["sw_variant"]))
    for k, c in corners.items():
        b, v = k.split("|")
        ents.append(C.r3_entry(f"corner_{b}_{v}_dTreq4", c["dT_req_4sol_K"], "K", src, fn, boundary=b, variant=v))
        ents.append(C.r3_entry(f"corner_{b}_{v}_dTreq5", c["dT_req_5sol_K"], "K", src, fn, boundary=b, variant=v))
    w = worst["kirchhoff|lw"]
    ents.append(C.r3_entry("worstcase_dTreq4", w["dT_req_4sol_K"], "K", src, fn, boundary="kirchhoff", variant="lw",
                           note="esquina lw con Delta_theta - sqrt(2)*0.25 K"))
    ents.append(C.r3_entry("worstcase_ratio_contrast_over_dTreq4", w["ratio_contrast_over_dT_req4"], "-", src, fn,
                           boundary="kirchhoff", variant="lw"))
    C.stage_numbers(ents, ("hasi_sigma_dtheta_independent_brief", "robust_", "corner_", "worstcase_"), fn)

    for r in rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if k not in ("note",)})
    print("R2 1.76 corresponde a:", corresponds)
    print("peor caso:", {k: (round(v["dT_req_4sol_K"], 4), round(v["ratio_contrast_over_dT_req4"], 4)) for k, v in worst.items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
