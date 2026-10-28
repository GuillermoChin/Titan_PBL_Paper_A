#!/usr/bin/env python
"""
make_tables.py -- Tables 1-8 of the revised Paper A and the auxiliary tables.

Every table is written in long format (one row per cell):

    row, column, value, print, source

where "value" is the full value, "print" the printed form, and "source" the
key of numbers/numbers_paper_A.json or the build/ file and field the value
comes from. No number is typed here: values come from the single source,
from the analysis outputs under build/, or from the configuration preset.

Auxiliary tables (tables/aux_*.csv) are copies of the analysis outputs that
support numbers quoted in the text without a key of their own (sensitivity
study, threshold verification, heating profiles, per-sol archive readings,
longwave attribution, stability budget, energy partition, column closure per
sol, flux comparison per interface).

Usage:
    python code/make_tables.py
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
sys.path.insert(0, str(ROOT / "code" / "analysis"))
import r3_common as C  # noqa: E402  (guard: titan_pbl del subconjunto)

from build_numbers_paper_A import fmt_sig  # noqa: E402
from titan_pbl.config import load_config  # noqa: E402

TABLES = ROOT / "tables"
BUILD = ROOT / "build"
NUM = {e["key"]: e for e in json.loads(C.NUMBERS_REV.read_text(encoding="utf-8"))["numbers"]}
BOUNDARIES = (("kirchhoff", "Kirchhoff-consistent"), ("full_absorption", "Full absorption"))
SW = (("lw", "none"), ("noon", "noon"), ("diurnal", "diurnal mean"))
RUNS = (("TMN91", "Baseline (1D_Huygens_TMN91_5sol)"), ("baseline5sol", "Baseline (1D_Huygens_baseline_5sol)"),
        ("0p1flhc", "Suppressed flux (1D_Huygens_0.1flhc_5sol)"))


def key(row: str, col: str, k: str) -> dict:
    e = NUM[k]
    if e["status"] != "valid":
        raise SystemExit(f"clave superseded en una tabla: {k}")
    return {"row": row, "column": col, "value": e["value"], "print": e["print"], "source": f"numbers:{k}"}


def val(row: str, col: str, v, src: str, pr: str | None = None) -> dict:
    if pr is None:
        pr = fmt_sig(float(v), 2) if isinstance(v, (int, float)) and not isinstance(v, bool) else str(v)
    return {"row": row, "column": col, "value": v, "print": pr, "source": src}


def write(name: str, rows: list[dict]) -> None:
    with (TABLES / name).open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["row", "column", "value", "print", "source"], lineterminator="\n")
        w.writeheader()
        for r in rows:
            v = r["value"]
            w.writerow({**r, "value": repr(v) if isinstance(v, float) else v})


def table01() -> list[dict]:
    src = BUILD / "r3" / "t5_tau_rad_table.csv"
    rows = []
    tr = {"calculo propio": "own calculation", "dado en anos de Titan": "given in Titan years",
          "anos terrestres / 29.5 (Strobel 2009, p. 242)": "Earth years / 29.5 (Strobel et al. 2009, p. 242)",
          "sin valor numerico": "no numerical value"}
    v0 = (("v0 Sec. 2.4", "original submission, Sec. 2.4"), ("v0 Eq. 6", "original submission, Eq. 6"))
    for r in C.read_csv(src):
        label = r["source"] + (f" ({r['page']})" if r["page"] not in ("", "-") else "")
        definition = r["definition"] + "; " + r["region_or_inputs"]
        for a, b in v0:
            label, definition = label.replace(a, b), definition.replace(a, b)
        # value_as_written ya trae la unidad en los calculos propios ("6.472291e+08 s").
        value, units = r["value_as_written"], r["units_as_written"]
        if units not in ("", "-") and not value.endswith(" " + units):
            value = f"{value} {units}"
        rows.append(val(label, "definition or region", definition,
                        f"build/r3/t5_tau_rad_table.csv:definition,region_or_inputs", pr=""))
        rows.append(val(label, "value as reported", value,
                        "build/r3/t5_tau_rad_table.csv:value_as_written,units_as_written", pr=""))
        if r["value_titan_years"]:
            rows.append(val(label, "Titan years", float(r["value_titan_years"]),
                            "build/r3/t5_tau_rad_table.csv:value_titan_years"))
        rows.append(val(label, "conversion", tr.get(r["conversion_method"], r["conversion_method"]),
                        "build/r3/t5_tau_rad_table.csv:conversion_method", pr=""))
    for r in rows:
        if r["print"] == "":
            r["print"] = r["value"]
    return rows


def table02() -> list[dict]:
    oat = C.read_csv(BUILD / "r2" / "t5_sensitivity_oat.csv")
    tau = sorted({float(r["value"]) for r in oat if r["factor"] == "radiation.optical_depth_ir"})
    f = sorted({float(r["value"]) for r in oat if r["factor"] == "radiation.ir_weighting_coefficient"})
    d_heat = NUM["threshold_heating_layer0_full_absorption_lw"]["value"] - NUM["threshold_heating_layer0_kirchhoff_lw"]["value"]
    u_hasi = NUM["Ta_layer0"]["uncertainty"]
    t = [
        ("Gray, single-band longwave transfer", "Neglects spectral windows and band structure; direction not determined here",
         f"Effective opacity varied through τ_IR ({tau[0]:g}–{tau[-1]:g}) and f ({f[0]:.2f}, {f[-1]:.2f}), Section 3.5"),
        ("Two-stream integration without a diffusivity factor", "Equivalent to a rescaling of the infrared optical depth",
         "Covered by the range of τ_IR"),
        ("No longwave scattering", "As in Chatain et al. (2022)", "Not assessed"),
        ("Frozen atmospheric state, no time integration",
         "Ignores that the lower layers warm and emit more as they destabilize, which slows further destabilization; mixing times are lower bounds",
         "Conservative"),
        ("Contrast held constant through the sol", "Real contrasts peak for a few hours; destabilization is overestimated",
         "Conservative; the real diurnal cycle is evaluated in Section 3.4"),
        ("Shortwave treatment", "The longwave-only treatment gives the smallest required contrast",
         "Longwave only in the text; noon and diurnal mean in Table 4 and Figure S1"),
        ("Surface boundary condition", f"Full absorption raises all thresholds by about {fmt_sig(d_heat, 2)} K",
         "Both reported; Kirchhoff-consistent values in the text"),
        ("HASI profile as initial state", f"Temperature uncertainty of ±{u_hasi:g} K per point", "Propagated to Δθ, Section 3.5"),
        ("Mixing criterion", "Radiative erosion of Δθ over 0–4 km; other mixing processes not represented", "Definition, Eq. (10)"),
    ]
    src = ("build/r2/t5_sensitivity_oat.csv (sweep design); numbers:threshold_heating_layer0_*_lw (difference); "
           "numbers:Ta_layer0 (uncertainty)")
    rows = []
    for a, eff, tr in t:
        rows += [val(a, "effect on the results", eff, src, pr=eff), val(a, "treatment", tr, src, pr=tr)]
    return rows


def table03() -> list[dict]:
    cfg = load_config("huygens_baseline")
    r, a, p = cfg.radiation, cfg.atmosphere, cfg.planet
    items = [
        ("Solar constant", "S_⊙", r.solar_constant, "W m-2", "radiation.solar_constant", "Hartwick et al. (2026)"),
        ("Downward IR at model top", "I_IR", r.ir_top_flux, "W m-2", "radiation.ir_top_flux", "Hartwick et al. (2026)"),
        ("Surface albedo", "α", r.surface_albedo, "-", "radiation.surface_albedo", "Hartwick et al. (2026)"),
        ("Surface emissivity", "ε", r.surface_emissivity, "-", "radiation.surface_emissivity", "Hartwick et al. (2026)"),
        ("Visible optical depth", "τ_vis", r.optical_depth_visible, "-", "radiation.optical_depth_visible", "Chatain et al. (2022)"),
        ("Infrared optical depth", "τ_IR", r.optical_depth_ir, "-", "radiation.optical_depth_ir", "Hartwick et al. (2026)"),
        ("Visible scaling exponent", "γ", r.visible_scaling_exponent, "-", "radiation.visible_scaling_exponent", "Schneider et al. (2012)"),
        ("IR weighting coefficient", "f", r.ir_weighting_coefficient, "-", "radiation.ir_weighting_coefficient", "Chatain et al. (2022)"),
        ("Single scattering albedo", "ω_0", r.single_scattering_albedo, "-", "radiation.single_scattering_albedo", "Chatain et al. (2022)"),
        ("Asymmetry factor", "g_a", r.asymmetry_factor, "-", "radiation.asymmetry_factor", "Chatain et al. (2022)"),
        ("Atmospheric layers", "N", a.vertical_grid.n_layers, "-", "atmosphere.vertical_grid.n_layers", "Hartwick et al. (2026)"),
        ("Surface pressure", "p_s", a.vertical_grid.surface_pressure, "Pa", "atmosphere.vertical_grid.surface_pressure", "Hartwick et al. (2026)"),
        ("Model top pressure", "p_t", a.vertical_grid.top_pressure, "Pa", "atmosphere.vertical_grid.top_pressure", "Hartwick et al. (2026)"),
        ("Specific heat", "c_p", a.gas.specific_heat_cp, "J kg-1 K-1", "atmosphere.gas.specific_heat_cp", "Hartwick et al. (2026)"),
        ("Surface gravity", "g", p.gravity, "m s-2", "planet.gravity", "standard"),
        ("Length of the Titan sol", "-", p.sol_length, "s", "planet.sol_length", "archive of Hartwick et al. (2026)"),
    ]
    def pressure(x: float) -> str:
        e = len(str(int(x))) - 1
        return f"{x / 10**e:.2f} × 10^{e}"

    rows = []
    for name, sym, v, units, path, source in items:
        pr = (pressure(float(v)) if path.endswith("pressure") else f"{int(v):,}".replace(",", " ") if path.endswith("sol_length")
              else f"{v:g}" if isinstance(v, float) else str(v))
        rows += [val(name, "symbol", sym, f"preset huygens_baseline:{path}", pr=sym),
                 val(name, "value", v, f"preset huygens_baseline:{path}", pr=pr),
                 val(name, "units", units, f"preset huygens_baseline:{path}", pr=units),
                 val(name, "source", source, f"preset huygens_baseline:{path}", pr=source)]
    rows += [val("Surface longwave boundary condition", "value", "Kirchhoff-consistent; full absorption for comparison",
                 "radiation.surface_lw_boundary = kirchhoff | full_absorption (Section 2.2)",
                 pr="Kirchhoff-consistent; full absorption for comparison")]
    return rows


def table04() -> list[dict]:
    rows = []
    for b, bl in BOUNDARIES:
        for v, vl in SW:
            lab = f"{bl}, shortwave {vl}"
            rows += [key(lab, "Heating threshold (K)", f"threshold_heating_layer0_{b}_{v}"),
                     key(lab, "Destabilization threshold (K)", f"threshold_destabilisation_0_4km_{b}_{v}"),
                     key(lab, "ΔT_req(4) (K)", f"dTreq_4sol_{b}_{v}"),
                     key(lab, "ΔT_req(5) (K)", f"dTreq_5sol_{b}_{v}"),
                     key(lab, "t_mix at 3 K (sols)", f"mixing_time_3K_{b}_{v}")]
    return rows


def table05() -> list[dict]:
    dth, sig = NUM["delta_theta_0_4km"]["value"], NUM["hasi_sigma_dtheta_independent_brief"]["value"]
    cases = (("Nominal (Table 3)", "preset", "nominal", dth), ("Δθ − σ_Δθ", "preset", "minus_sigma_independent", dth - sig),
             ("Δθ + σ_Δθ", "preset", "plus_sigma_independent", dth + sig), ("Most favorable corner", "corner_lw", "nominal", dth),
             ("Corner with Δθ − σ_Δθ", "corner_lw", "minus_sigma_independent", dth - sig))
    rows = []
    for lab, cfg, case, d in cases:
        rows.append(val(lab, "Δθ (K)", d, "numbers:delta_theta_0_4km ± numbers:hasi_sigma_dtheta_independent_brief"))
        for n in (4, 5):
            for b, bl in BOUNDARIES:
                rows.append(key(lab, f"ΔT_req({n}), {bl} (K)", f"robust_dTreq{n}_{b}_lw_{cfg}_{case}"))
    return rows


def table06() -> list[dict]:
    summ = json.loads((BUILD / "r3" / "t3_flux_comparison_all.json").read_text(encoding="utf-8"))["summary"]
    prefix = {"TMN91": "flux_comparison_", "baseline5sol": "flux_comparison_baseline5sol_", "0p1flhc": "flux_comparison_0p1flhc_"}
    rows = []
    for run, lab in RUNS:
        p = prefix[run]
        s = summ["0.1flhc" if run == "0p1flhc" else run]
        rows += [key(lab, "Records", f"{p}n_records"),
                 key(lab, "Max. relative RMS, F↓ (%)", f"{p}max_rms_rel_F_down_pct"),
                 key(lab, "Max. relative RMS, F↑ (%)", f"{p}max_rms_rel_F_up_pct"),
                 val(lab, "Lowest-layer rate, archive (K sol-1)", s["heating_layer0_archive_K_per_sol_mean"],
                     "build/r3/t3_flux_comparison_all.json:summary.heating_layer0_archive_K_per_sol_mean"),
                 val(lab, "Lowest-layer rate, this work (K sol-1)", s["heating_layer0_titan_K_per_sol_mean"],
                     "build/r3/t3_flux_comparison_all.json:summary.heating_layer0_titan_K_per_sol_mean")]
    return rows


def table08() -> list[dict]:
    rows = []
    for run, lab in RUNS:
        rows += [key(lab, "Median HFX/F_MYJ", f"h33_median_r_{run}"),
                 key(lab, "∫HFX dt (MJ m-2)", f"h33_int_HFX_{run}"),
                 key(lab, "∫F_MYJ dt (MJ m-2)", f"h33_int_FMYJ_{run}"),
                 key(lab, "Residual with HFX (MJ m-2)", f"h33_colresid_HFX_{run}"),
                 key(lab, "Residual with F_MYJ (MJ m-2)", f"h33_colresid_FMYJ_{run}")]
    return rows


def table07() -> None:
    """Table 7: the literature and archive quantities of the single source, without the
    values of the original submission that were withdrawn (ids v0_*)."""
    with (C.DATA_REV / "literature_contrasts.csv").open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        fields, rows = reader.fieldnames, [r for r in reader if not r["id"].startswith("v0_")]
    with (TABLES / "table07_literature_contrasts.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


AUX = [
    ("aux_oat_sensitivity.csv", BUILD / "r2" / "t5_sensitivity_oat.csv"),
    ("aux_threshold_verification_384.csv", BUILD / "r1" / "t4_thresholds.csv"),
    ("aux_thresholds_sw_variants.csv", BUILD / "source_a" / "thresholds_sw_variants.csv"),
    ("aux_mixing_times_sw_variants.csv", BUILD / "source_a" / "mixing_times_sw_variants.csv"),
    ("aux_heating_profiles.csv", BUILD / "source_a" / "fig2_profiles_sw_variants.csv"),
    ("aux_archive_per_sol.csv", BUILD / "source_a" / "kirchhoff" / "table4_contrasts.csv"),
    ("aux_lw_attribution.csv", BUILD / "r2" / "t3_attribution.csv"),
    ("aux_stability_budget.csv", BUILD / "r2" / "t3_stability_budget.csv"),
    ("aux_energy_partition.csv", BUILD / "r2" / "t2_energy_partition.csv"),
    ("aux_column_closure_per_sol.csv", BUILD / "r3" / "t1_h33_sol_integrals.csv"),
    ("aux_flux_comparison_interfaces.csv", BUILD / "r3" / "t3_flux_comparison_all.csv"),
]
# Celdas de texto de las salidas de analisis, traducidas por valor exacto de la celda.
AUX_TRANSLATE = {
    "barrido declarado": "declared sweep",
    "complemento descriptivo (no entra en la esquina)": "descriptive complement (not used for the corner)",
    "30 centrales (k=10..39)": "30 central days (k=10..39)",
    "50 dias": "50 days",
    "- (condicional)": "- (conditional)",
    "- (condicional, lectura del balance)": "- (conditional, read from the balance)",
    "MJ m-2 per Titan day (media de los dias)": "MJ m-2 per Titan day (mean over the days)",
    "No determinable": "not determinable",
    "dia": "day",
    "noche": "night",
    "todo": "all",
    "apoya": "supports",
    "no apoya": "does not support",
    "parcial": "partial",
    "intermedio (1 < |t_eros| < 10 sols)": "intermediate (1 < |t_eros| < 10 sols)",
    "no erosiona: <dS_LW> > 0 (el cambio LW estabiliza)": "no erosion: <dS_LW> > 0 (the longwave change stabilizes)",
    "no puede actuar en escala diurna (|t_eros| >= 10 sols)": "cannot act on the diurnal scale (|t_eros| >= 10 sols)",
}


def write_aux(name: str, src: Path) -> None:
    with src.open(encoding="utf-8", newline="") as fh:
        rows = list(csv.reader(fh))
    with (TABLES / name).open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerows([[AUX_TRANSLATE.get(c, c) for c in r] for r in rows])


def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    write("table01_radiative_timescales.csv", table01())
    write("table02_assumptions.csv", table02())
    write("table03_configuration.csv", table03())
    write("table04_thresholds_required_contrasts.csv", table04())
    write("table05_hasi_uncertainty_corner.csv", table05())
    write("table06_flux_comparison.csv", table06())
    table07()
    write("table08_column_energy_closure.csv", table08())
    for name, src in AUX:
        write_aux(name, src)
    print("tablas:", sorted(p.name for p in TABLES.glob("*.csv")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
