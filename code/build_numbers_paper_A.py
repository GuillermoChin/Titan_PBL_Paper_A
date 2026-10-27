#!/usr/bin/env python
"""
build_numbers_paper_A.py -- single source of the numbers of the revised Paper A.

Consolidates the staged numbers of the three analysis stages

  build/source_a/numbers_r1.json   (framework diagnostics and archive readings)
  build/r2/numbers_r2.json         (archive and TAM analyses)
  build/r3/numbers_r3_staging.json (column closure, robustness, timescales)

into numbers/numbers_paper_A.json, proposes the printed form of every number
with a one-line rationale, and writes numbers/numbers_review.md (readable
table) and build/numbers/literature_contrasts.csv (Table 7).

Rules
- Each entry: key, value, uncertainty, units, lw_boundary, sw_variant, print,
  sig_figs_rationale, source_file, source_function, stage, status,
  superseded_by, mandatory, section, note.
- Rounding: no printed number carries more precision than its uncertainty or
  its least precise input justifies. "print" is the proposed form; the form
  printed in the manuscript is the author's decision (numbers_review.md).
- When two stages give a value to the same quantity (same key, or a pair in
  EQUIVALENCES), the values must agree to 1e-9 (relative); otherwise a
  documented decision states which one is valid and the other is kept with
  status "superseded". An unresolved conflict stops the build.
- Text fields are written in English (TRANSLATIONS); the analysis code keeps
  its comments and internal labels in Spanish.

Format of "print": the figure only (no units), "a × 10^n" for |x| >= 1e5,
booleans and texts as they are, "never" for infinite times.
"""

from __future__ import annotations

import json
import math
import os
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "code" / "analysis"))
import r3_common as C  # noqa: E402  (guard: titan_pbl del subconjunto)

# PAPER_A_NUMBERS_OUT_DIR redirige todas las salidas (lo usa audit_paper_A_numbers.py para no tocar numbers/).
_ALT = os.environ.get("PAPER_A_NUMBERS_OUT_DIR")
OUT_JSON = Path(_ALT) / "numbers_paper_A.json" if _ALT else C.NUMBERS_REV
OUT_REVIEW = Path(_ALT) / "numbers_review.md" if _ALT else C.NUMBERS_DIR / "numbers_review.md"
OUT_CONTRASTS = Path(_ALT) / "literature_contrasts.csv" if _ALT else C.DATA_REV / "literature_contrasts.csv"
OUT_CONSOLIDATION = Path(_ALT) / "consolidation.json" if _ALT else C.DATA_REV / "consolidation.json"

# --------------------------------------------------------------------------
# Equivalencias (misma cantidad) y su resolucion
# --------------------------------------------------------------------------
EQUIVALENCES = [
    # (clave A, clave B, resolucion si difieren)
    ("tau_rad_eq6", "tau_rad_r2_eq6_s", None),
    ("tau_rad_eq6_sols", "tau_rad_r2_eq6_sols", None),
    ("tau_rad_eq6_years", "tau_rad_r2_eq6_years",
     ("D21", "taurad_tau_lin_titan_years",
      "the stage-r1 value used the year of the preset (no source); the Titan year is taken from Strobel et al. (2009, p. 242)")),
    ("tau_rad_eq6", "taurad_tau_lin_s", None),
    ("tau_rad_eq6_sols", "taurad_tau_lin_Tsol", None),
    ("tau_rad_eq6_years", "taurad_tau_lin_titan_years",
     ("D21", "taurad_tau_lin_titan_years", "Titan year with a source (29.5 Earth years) instead of the preset")),
    ("tau_rad_r2_stone_form_s", "taurad_tau_bulk_s", None),
    ("hasi_unc_sigma_dtheta_independent", "hasi_sigma_dtheta_independent_brief",
     ("D20", "hasi_sigma_dtheta_independent_brief",
      "sigma = sqrt(2) x 0.25 K (independent errors); the theta/T-weighted value is kept as a sensitivity")),
    ("corner_dTreq_4sol", "corner_kirchhoff_lw_dTreq4", None),
    ("corner_dTreq_5sol", "corner_kirchhoff_lw_dTreq5", None),
    ("dTreq_4sol_kirchhoff_lw", "robust_dTreq4_kirchhoff_lw_preset_nominal", None),
    ("dTreq_5sol_kirchhoff_lw", "robust_dTreq5_kirchhoff_lw_preset_nominal", None),
    ("dTreq_4sol_full_absorption_lw", "robust_dTreq4_full_absorption_lw_preset_nominal", None),
    ("dTreq_5sol_full_absorption_lw", "robust_dTreq5_full_absorption_lw_preset_nominal", None),
    ("colenergy_residual_HFX_full_column_TMN91", "h33_colresid_HFX_TMN91", None),
    ("colenergy_residual_HFX_full_column_baseline5sol", "h33_colresid_HFX_baseline5sol", None),
    ("colenergy_residual_HFX_full_column_0p1flhc", "h33_colresid_HFX_0p1flhc", None),
]
for run in ("TMN91", "baseline5sol", "0p1flhc"):
    EQUIVALENCES.append((f"colenergy_HFX_over_H_AKHS_{run}", f"h33_median_r_{run}",
                         ("exact MYJ flux", f"h33_median_r_{run}",
                          "the exact MYJ lower-boundary flux (R_D = 290, theta with P0 = 1e5) replaces the stage-r2 approximation")))
    EQUIVALENCES.append((f"colenergy_residual_H_AKHS_full_column_{run}", f"h33_colresid_FMYJ_{run}",
                         ("exact MYJ flux", f"h33_colresid_FMYJ_{run}", "exact MYJ flux instead of the stage-r2 approximation")))

# Claves que usa el manuscrito revisado.
MANDATORY = {
    "dTreq_4sol_kirchhoff_lw", "dTreq_4sol_full_absorption_lw", "archive_contrast_max_sols2to5_suppressed",
    "archive_contrast_max_sols2to5_baseline", "taurad_tau_bulk_titan_years", "taurad_tau_lin_titan_years",
    "worstcase_dTreq4", "h33_median_r_0p1flhc", "flux_comparison_max_rms_rel_F_down_pct",
    "mayor3_Q_S_4000m", "delta_theta_0_4km",
}

# Secciones del manuscrito revisado (numeracion final).
SECTION_RULES = [
    (r"^(analytic_threshold|Ta_layer0)", "§2.3 The heating threshold"),
    (r"^(tau_rad|taurad|titan_year)", "§2.5 The radiative timescale"),
    (r"^(z_mid|thickness|d_tau)", "§3.1 Method"),
    (r"^(profile_|layer0_heating)", "§3.2 The heating rate profile"),
    (r"^threshold_", "§3.3 Threshold contrasts"),
    (r"^(delta_theta|mixing_time|destabilisation_rate|dTreq|swmatrix|headline|realdiurnal_)",
     "§3.4 The mixing timescale and the required contrast"),
    (r"^(hasi_|robust_|corner_|worstcase)", "§3.5 Sensitivity and uncertainty"),
    (r"^flux_comparison", "§3.6 Comparison with archived output"),
    (r"^archive_contrast", "§4.1 Verified contrasts"),
    (r"^archive_", "§5.1 The experiment and the proposed mechanism"),
    (r"^mayor3_", "§5.2–5.3 Radiative heating and static stability"),
    (r"^(colenergy_|h33_|mtwrf_)", "§5.4 Column energy closure"),
    (r"^(energy_archive|energy_tam)", "§5.5 Implications for comparing mtWRF and TAM"),
    (r"^sweep_", "Other (parameter sweep of the original submission)"),
]
SECTION_ORDER = [s for _, s in SECTION_RULES if not s.startswith("Other")] + [
    "Other (parameter sweep of the original submission)", "Other"]


def section_of(key: str) -> str:
    for rx, sec in SECTION_RULES:
        if re.match(rx, key):
            return sec
    return "Other"


# Textos de las etapas (en espanol en el codigo de analisis) -> ingles en la fuente unica.
TRANSLATIONS = [
    ("Titan sols (inf = nunca)", "Titan sols (inf = never)"),
    ("n/a (archivo)", "n/a (archive)"),
    ("ciclo diurno real", "real diurnal cycle"),
    ("n/a (solo LW)", "n/a (longwave only)"),
    ("insumos v0 Sec. 2.4; D21", "inputs of Section 2.4 of the original submission; two explicit definitions (bulk, linearized)"),
    ("mediana de HFX/F_MYJ (diurno, |F|>0.05); uncertainty = RIC",
     "median of HFX/F_MYJ over daytime records with |F_MYJ| > 0.05 W m-2; uncertainty = interquartile range"),
    ("sqrt(2)*0.25 K (D20)", "sqrt(2) x 0.25 K (independent errors at the surface and at 4 km)"),
    ("esquina lw con Delta_theta - sqrt(2)*0.25 K", "longwave-only corner with Delta_theta - sqrt(2) x 0.25 K"),
    ("veredicto 6.2 (arquitectura WRF publica, 4 versiones)", "verdict on the public WRF architecture (four versions)"),
    ("veredicto 6.3", "verdict on the archived ratio HFX/F_MYJ"),
    ("R_D = kappa*CP, kappa ajustado de theta = T + T0 en el archivo",
     "R_D = kappa x CP, kappa fitted from theta = T + T0 in the archive"),
    ("theta = T + T0 en el wrfout", "theta = T + T0 in the wrfout files"),
    ("derivacion: con F_up = eps sigma Ts^4 + (1-eps) F_down y F_down ~ sigma Ta^4, dE/dt ~ dtau eps sigma (Ts^4 - Ta^4)",
     "derivation: with F_up = eps sigma Ts^4 + (1-eps) F_down and F_down ~ sigma Ta^4, dE/dt ~ dtau eps sigma (Ts^4 - Ta^4)"),
    ("(D12: argmin dT_req(4) con kirchhoff)", "(argmin of dT_req(4) with the Kirchhoff-consistent condition)"),
    ("(Ec. 6 de v0: cp p / (4 g sigma T^3), insumos de v0 2.4)",
     "(Eq. 6 of the original submission: cp p / (4 g sigma T^3), inputs of its Section 2.4)"),
    ("(Tsol de v0 Sec. 3.1)", "(Titan sol of Section 3.1 of the original submission)"),
    ("(cp p / (g sigma T^3), forma de Stone 1975 / Tokano 1999)",
     "(cp p / (g sigma T^3), form of Stone 1975 / Tokano et al. 1999)"),
    ("Ec. 6: cp p / (4 g sigma T^3), T = 93 K", "Eq. 6 of the original submission: cp p / (4 g sigma T^3), T = 93 K"),
    ("Ec. 6 / planet.sol_length", "Eq. 6 of the original submission / planet.sol_length"),
    ("Ec. 6 / planet.orbital_period", "Eq. 6 of the original submission / planet.orbital_period"),
    ("depende de z_ref", "depends on z_ref"),
    ("Escalado confirmado en la salida", "Scaling confirmed in the output"),
    ("Consistente", "Consistent"),
    ("texto", "text"),
]


def en(s):
    if not isinstance(s, str):
        return s
    for a, b in TRANSLATIONS:
        if s == a or (len(a) > 12 and a in s):
            s = s.replace(a, b)
    return s


# --------------------------------------------------------------------------
# Regla de redondeo
# --------------------------------------------------------------------------
def fmt_sig(x: float, sig: int) -> str:
    if x == 0:
        return "0"
    if abs(x) >= 1e5:
        e = int(math.floor(math.log10(abs(x))))
        m = round(x / 10**e, sig - 1)
        if abs(m) >= 10:
            m, e = m / 10, e + 1
        return f"{m:.{sig - 1}f} × 10^{e}"
    d = sig - 1 - int(math.floor(math.log10(abs(x))))
    v = round(x, d)
    if d <= 0:
        return f"{int(round(v)):d}"
    return f"{v:.{d}f}"


def fmt_dec(x: float, decimals: int) -> str:
    return f"{round(x, decimals):.{decimals}f}"


def propose(e: dict) -> tuple[str, str]:
    k, v, units = e["key"], e["value"], str(e.get("units", ""))
    if isinstance(v, bool):
        return str(v).lower(), "boolean: not numerical"
    if isinstance(v, str):
        return v, "text: not numerical"
    if v is None:
        if "nunca" in units or "inf" in units:
            return "never", "infinite time: the rate is not positive"
        return "n/a", "no value (not determinable at its stage)"
    if isinstance(v, int) and not isinstance(v, bool) and ("record" in units):
        return str(v), "exact count"
    x = float(v)
    if re.match(r"^(tau_rad|taurad)", k):
        return fmt_sig(x, 2), "2 significant figures: T = 93 K (2 figures) is the least precise input of the timescale"
    if k == "titan_year_earth_years":
        return "29.5", "as written in the source (Strobel et al. 2009, p. 242)"
    if k == "Ta_layer0":
        return fmt_dec(x, 1), "0.1 K: HASI uncertainty of ±0.25 K (Fulchignoni et al. 2005)"
    if k == "delta_theta_0_4km":
        return fmt_sig(x, 2), "2 significant figures: sigma = 0.35 K fixes the hundredths"
    if k.startswith("hasi_"):
        return fmt_sig(x, 2), "uncertainty: 2 significant figures"
    if re.match(r"^(archive_contrast)", k) or k.startswith("worstcase_ratio"):
        return fmt_sig(x, 2), "2 significant figures, as the thresholds they are compared with"
    if re.match(r"^(threshold_|dTreq|robust_|corner_|worstcase_dTreq|analytic_threshold)", k):
        return fmt_sig(x, 2), ("2 significant figures: eps is given with 2 figures and a change of 0.05 moves "
                               "~0.3 K; with independent HASI errors the uncertainty is ~±2 K (reported as a range)")
    if k.startswith(("mixing_time", "swmatrix_mix")):
        return fmt_sig(x, 2), "2 significant figures: rate with 2 figures; independent HASI errors give ±54 % (range)"
    if "per Titan sol" in units or "K per" in units:
        return fmt_sig(x, 2), "2 significant figures: rates; eps with 2 figures"
    if k.startswith("flux_comparison") and units == "%":
        return fmt_sig(x, 1), "1 significant figure suffices for the order against the 1 % criterion"
    if k == "flux_comparison_implied_archive_I_IR":
        return fmt_sig(x, 3), "3 significant figures: geometry check (2.00)"
    if k.startswith("flux_comparison"):
        return fmt_sig(x, 2), "2 significant figures"
    if k.startswith("h33_median_r"):
        return fmt_sig(x, 3), "3 significant figures: IQR <= 4e-4 and criterion tolerance 0.005"
    if k.startswith(("h33_colresid", "h33_int", "colenergy_residual")):
        return fmt_dec(x, 1), "0.1 MJ m-2: the closure of the baselines leaves ~0.2 MJ m-2 of numerical residual"
    if k.startswith("colenergy_HFX_over"):
        return fmt_sig(x, 2), "stage-r2 definition (parallel); 2 significant figures"
    if units == "%":
        return fmt_sig(x, 2), "percentage of archive readings: 2 significant figures"
    if "W m-2" in units:
        return fmt_sig(x, 3), "archive reading; 3 significant figures as in the original submission"
    if "MJ" in units:
        return fmt_sig(x, 2), "integrals per sol: 2 significant figures"
    if k.startswith("mayor3_"):
        return fmt_sig(x, 2), "2 significant figures: cp ±5 % moves Q_S(4 km) by ±0.005"
    if k == "mtwrf_R_D_inferred":
        return fmt_sig(x, 3), "3 significant figures: exact fit (290.00)"
    if k == "mtwrf_theta_T0":
        return fmt_sig(x, 3), "model constant"
    if k in ("z_mid_layer0", "thickness_layer0", "d_tau_layer0"):
        return fmt_sig(x, 2), "grid geometry: 2 significant figures suffice"
    return fmt_sig(x, 2), "default rule: 2 significant figures"


# --------------------------------------------------------------------------
def load_stage(path: Path, stage: str) -> list[dict]:
    out = []
    for n in C.numbers(path):
        e = {"key": n["key"], "value": n["value"], "uncertainty": n.get("uncertainty"), "units": n["units"],
             "lw_boundary": n.get("lw_boundary"), "sw_variant": n.get("sw_variant"),
             "source_file": n["source_file"], "source_function": n["source_function"], "stage": stage}
        if n.get("note"):
            e["note"] = n["note"]
        out.append(e)
    return out


def consolidate() -> tuple[list[dict], dict]:
    entries = (load_stage(C.NUMBERS_R1, "r1") + load_stage(C.NUMBERS_R2, "r2") +
               load_stage(C.NUMBERS_R3, "r3"))
    by_key: dict[str, list[dict]] = {}
    for e in entries:
        by_key.setdefault(e["key"], []).append(e)
    same_key = {k: v for k, v in by_key.items() if len(v) > 1}
    checks, unresolved = [], []
    for k, v in same_key.items():
        vals = {json.dumps(x["value"]) for x in v}
        status = "agree" if len(vals) == 1 else "unresolved conflict"
        checks.append({"a": k, "b": k, "stages": [x["stage"] for x in v], "status": status})
        if status != "agree":
            unresolved.append(k)
    idx = {e["key"]: e for e in entries}
    for a, b, res in EQUIVALENCES:
        ea, eb = idx.get(a), idx.get(b)
        if ea is None or eb is None:
            checks.append({"a": a, "b": b, "status": "missing key", "present": [ea is not None, eb is not None]})
            continue
        va, vb = ea["value"], eb["value"]
        same = (va == vb) or (isinstance(va, (int, float)) and isinstance(vb, (int, float))
                              and abs(va - vb) <= 1e-9 * max(abs(va), abs(vb), 1e-300))
        row = {"a": a, "b": b, "stage_a": ea["stage"], "stage_b": eb["stage"], "value_a": va, "value_b": vb}
        if same:
            row["status"] = "agree"
        elif res is not None:
            decision, winner, why = res
            loser = a if winner == b else b
            row.update({"status": "resolved by decision", "decision": decision, "valid_key": winner, "why": why})
            idx[loser]["status"] = "superseded"
            idx[loser].setdefault("superseded_by", winner)
        else:
            row["status"] = "unresolved conflict"
            unresolved.append(f"{a} <> {b}")
        checks.append(row)
    for e in entries:
        e.setdefault("status", "valid")
        e["mandatory"] = e["key"] in MANDATORY
        e["print"], e["sig_figs_rationale"] = propose(e)
        if e["status"] == "superseded":
            e["sig_figs_rationale"] += f" | SUPERSEDED by {e['superseded_by']}"
        e["section"] = section_of(e["key"])
    # Incertidumbres conocidas de las entradas de la etapa r1 (declaradas)
    sig = idx["hasi_sigma_dtheta_independent_brief"]["value"]
    idx["delta_theta_0_4km"]["uncertainty"] = sig
    idx["delta_theta_0_4km"]["note"] = "uncertainty = sqrt(2) x 0.25 K (independent errors, upper bound); systematic: 0"
    idx["Ta_layer0"]["uncertainty"] = 0.25
    idx["Ta_layer0"]["note"] = "uncertainty = ±0.25 K (Fulchignoni et al. 2005, p. 787)"
    # Textos en ingles
    for e in entries:
        for f in ("value", "units", "lw_boundary", "sw_variant", "source_file", "source_function", "note"):
            if f in e:
                e[f] = en(e[f])
        if isinstance(e["value"], str):
            e["print"] = e["value"]
    return entries, {"checks": checks, "unresolved": unresolved, "stop": bool(unresolved), "n_entries": len(entries)}


# --------------------------------------------------------------------------
def write_review(entries: list[dict]) -> None:
    L = ["# Paper A (revised) — single source of numbers, readable table", "",
         "*Generated by `code/build_numbers_paper_A.py` from `numbers/numbers_paper_A.json`. The column `print` is "
         "the form proposed by the rounding rule; the form printed in the manuscript is the author's decision.*", "",
         "Rounding rule: no printed number carries more precision than its uncertainty or its least precise input "
         "justifies. Keys with status *superseded* are kept for traceability and are not used in the text.", ""]
    for sec in SECTION_ORDER:
        rows = [e for e in entries if e["section"] == sec]
        if not rows:
            continue
        L += [f"## {sec}", "", "| key | stage | full value | ± | units | print | rationale | status |",
              "|---|---|---|---|---|---|---|---|"]
        for e in rows:
            v = e["value"]
            vs = repr(v) if isinstance(v, float) else str(v)
            u = "" if e.get("uncertainty") is None else f"{e['uncertainty']:.3g}"
            st = e["status"] + (" · mandatory" if e["mandatory"] else "")
            L.append(f"| `{e['key']}` | {e['stage']} | {vs} | {u} | {e['units']} | {e['print']} | "
                     f"{e['sig_figs_rationale']} | {st} |")
        L.append("")
    OUT_REVIEW.parent.mkdir(parents=True, exist_ok=True)
    OUT_REVIEW.write_text("\n".join(L) + "\n", encoding="utf-8", newline="\n")


def write_contrasts() -> list[dict]:
    """Table 7: surface temperature quantities of the literature and of the two archives."""
    r2 = C.read_csv(C.R2_OUT / "t6_literature_contrasts.csv")
    tam = C.read_csv(C.R2_OUT / "t4_model_comparison.csv")
    rows = []
    cite = {
        "v0_generic_low": ("original submission, Section 3.3: 'of order 0.1 to 2 K'", "not found in the sources"),
        "v0_generic_high": ("original submission, Section 3.3: 'of order 0.1 to 2 K'", "not found in the sources"),
        "hanlora_amp_low": ("Han & Lora 2025, p. 5", "verified"),
        "hanlora_amp_high": ("Han & Lora 2025, p. 5", "verified"),
        "hartwick_gt1K": ("Hartwick et al. 2026, p. 10: 'increasing surface skin temperature by greater than one degree kelvin'", "verified"),
        "archive_baseline": ("Archive of Hartwick et al. 2026 (Zenodo 17246644), TMN91: max TSK − T_PHY[0], sols 2–5", "archive reading"),
        "archive_suppressed": ("Archive of Hartwick et al. 2026, 0.1flhc: max TSK − T_PHY[0], sols 2–5", "archive reading"),
        "archive_dTs": ("Archive of Hartwick et al. 2026: δTs between runs", "archive reading"),
        "charnay_theta_10m": ("Charnay & Lebonnois 2012, p. 106: 'decreased from 93.65 K at the ground to 93.48 K in the first 10 m'", "verified"),
        "charnay_Ts_amp": ("Charnay & Lebonnois 2012, p. 106: 'ground temperature varies by 0.4 K over one Titan day'", "verified"),
        "tokano2006_qual": ("Tokano et al. 2006: 'marginally warmer than the air'", "verified (qualitative)"),
    }
    magnitude = {
        "v0_generic_low": "surface-to-air contrast (stated in the original submission, no source)",
        "v0_generic_high": "surface-to-air contrast (stated in the original submission, no source)",
        "hanlora_amp_low": "amplitude of Ts (peak to peak)", "hanlora_amp_high": "amplitude of Ts (peak to peak)",
        "hartwick_gt1K": "other: increase of Ts between runs", "archive_dTs": "other: increase of Ts between runs",
        "archive_baseline": "surface-to-air contrast", "archive_suppressed": "surface-to-air contrast",
        "charnay_theta_10m": "surface-to-air contrast (θ of the ground − θ at 10 m, HASI)",
        "charnay_Ts_amp": "amplitude of Ts (model)", "tokano2006_qual": "other: qualitative",
    }
    is_contrast = {"archive_baseline": "yes", "archive_suppressed": "yes", "charnay_theta_10m": "yes"}
    for r in r2:
        i = r["id"]
        c, status = cite.get(i, ("", ""))
        rows.append({"source": r["estudio"], "id": i, "quantity": magnitude.get(i, r["naturaleza"]),
                     "value_K": r["contraste_K"], "Ta_height": r["altura_Ta"], "citation": c,
                     "citation_status": status, "surface_to_air_contrast": is_contrast.get(i, "no"),
                     "note": ("withdrawn: no source" if i.startswith("v0_generic") else
                              "Ta at the lowest mass level of the archive (2.7 m)" if i.startswith("archive_") and i != "archive_dTs"
                              else "classified as a contrast: θ of the air at the ground ≈ Ts" if i == "charnay_theta_10m"
                              else "")})
    for t in tam:
        if not t["modelo"].startswith("TAM"):
            continue
        rows.append({"source": f"TAM (Han & Lora 2025), {t['caso']}", "id": f"tam_{t['caso']}",
                     "quantity": "surface-to-air contrast (daily maximum, FMS means over 0.667 day)",
                     "value_K": t["Ts_menos_Ta_max_K"], "Ta_height": f"{float(t['altura_Ta_m']):.1f} m (lowest TAM level)",
                     "citation": "Zenodo 16791001, atmos_diurnal (tsurf − temp at the lowest level)",
                     "citation_status": "archive reading", "surface_to_air_contrast": "yes",
                     "note": "FMS time means; not in the original submission"})
    # Estimaciones bulk de Han & Lora (2025, p. 17), citas literales.
    hl_bulk = "T(z=0) − Tatm, bulk estimate (sensible heat flux 0.5–1.0 W m⁻² / a₀ ≈ 0.5 W m⁻² K⁻¹), Dragonfly site"
    hl_cite = 'Han & Lora 2025, p. 17: "this results in T(z = 0) − Tatm ≈ 1.0–2.0 K"'
    for bid, fuente, magnitud, valor, cita in (
            ("hanlora_bulk_dragonfly_low", "Han & Lora 2025", hl_bulk, "1.0", hl_cite),
            ("hanlora_bulk_dragonfly_high", "Han & Lora 2025", hl_bulk, "2.0", hl_cite),
            ("lorenz2021_bulk", "Lorenz 2021 (quoted by Han & Lora)",
             "T(z=0) − Tatm, simplified energy balance", "1.2",
             'Han & Lora 2025, p. 17: "T(z = 0) − Tatm ≈ 1.2 K (R. D. Lorenz 2021)"')):
        rows.append({"source": fuente, "id": bid, "quantity": magnitud, "value_K": valor,
                     "Ta_height": "not stated", "citation": cita, "citation_status": "verified",
                     "surface_to_air_contrast": "bulk",
                     "note": "bulk estimate; not a reading of a model or of an archive"})
    C.write_csv(OUT_CONTRASTS, rows)
    return rows


def main() -> int:
    C.DATA_REV.mkdir(parents=True, exist_ok=True)
    entries, cons = consolidate()
    report = {"head": C.git_head(), **cons}
    if cons["stop"]:
        C.dump_json(OUT_CONSOLIDATION, report)
        print("unresolved conflicts:", cons["unresolved"])
        return 3
    data = {"description": ("Single source of numbers for the revised version of 'Radiative Timescale Constraints on "
                            "Boundary Layer Evolution in Models of Titan' (Planetary and Space Science, PSS-D-26-00233)"),
            "generated_by": "code/build_numbers_paper_A.py",
            "sol_length_s": C.SOL_S,
            "print_rule": "proposed by the rounding rule of the generator; the printed form in the manuscript is the author's decision",
            "numbers": [{k: e.get(k) for k in ("key", "value", "uncertainty", "units", "lw_boundary", "sw_variant",
                                                "print", "sig_figs_rationale", "source_file", "source_function",
                                                "stage", "status", "superseded_by", "mandatory", "section", "note")}
                        for e in entries]}
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(C.clean_nan(data), indent=2, ensure_ascii=False, default=C.C2._json_default) + "\n",
                        encoding="utf-8", newline="\n")
    write_review(entries)
    contrasts = write_contrasts()
    report.update({"n_superseded": sum(e["status"] == "superseded" for e in entries),
                   "n_mandatory": sum(e["mandatory"] for e in entries),
                   "by_stage": {p: sum(e["stage"] == p for e in entries) for p in ("r1", "r2", "r3")},
                   "contrasts_rows": len(contrasts)})
    C.dump_json(OUT_CONSOLIDATION, report)
    print({k: v for k, v in report.items() if k not in ("checks",)})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
