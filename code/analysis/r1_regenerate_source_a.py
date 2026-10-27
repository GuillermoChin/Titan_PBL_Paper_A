#!/usr/bin/env python
"""
r1_regenerate_source_a.py -- stage r1, analysis T5.3-T5.7.

Regenerates the Source A quantities of Paper A for both surface longwave
boundaries and three shortwave variants, recomputes the archive contrasts
(decision D4), runs the declared consistency checks and writes the single
source of numbers ``build/source_a/numbers_r1.json``.

Todas las cantidades salen de funciones de ``code/analysis/source_a_core.py``
(las existentes y las nuevas de R1), importado como modulo. Las tablas de
siempre (umbrales, mezcla, Tabla 2, Tabla 4) las produce el propio script de
figuras con ``--lw-boundary <b> --outdir build/source_a/<b>`` (T5.1-T5.2);
este script las lee de ahi.

Variantes de onda corta:
  (a) "lw"       solo onda larga;
  (b) "noon"     onda larga + onda corta al mediodia del Tsol 1 (lo publicado);
  (c) "diurnal"  onda larga + promedio de 48 evaluaciones equiespaciadas de la
                 onda corta en el Tsol 1, sobre el estado fijo.

COMPROBACIONES DE CONSISTENCIA (criterio declarado):
  1. Con kirchhoff, columna isoterma (93.65 K, I_IR = sigma T^4, S0 = 0): el
     calentamiento es nulo en todas las capas, |.| <= 1e-12 K s-1.
  2. HASI al mediodia: el desplazamiento LW de la capa 0 entre K y O1 es
     g (1-eps) F_down(0) (1 - e^-dtau0) / (cp dp0) * sol, con tolerancia
     relativa 1e-6.

Salidas:
  build/source_a/thresholds_sw_variants.csv
  build/source_a/mixing_times_sw_variants.csv
  build/source_a/fig2_profiles_sw_variants.csv
  build/source_a/archive_contrasts_D4.csv
  build/source_a/consistency_checks.json
  build/source_a/numbers_r1.json
  build/r1/quicklook/*.png
"""

from __future__ import annotations

import csv
import importlib.util
import json
import os
import sys
import tempfile
from pathlib import Path

sys.dont_write_bytecode = True

import numpy as np  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT / "code" / "titan_pbl_subset") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "code" / "titan_pbl_subset"))

import titan_pbl  # noqa: E402

if (REPO_ROOT / "code" / "titan_pbl_subset") not in Path(titan_pbl.__file__).resolve().parents:
    raise ImportError(f"titan_pbl imported from {titan_pbl.__file__}, not from code/titan_pbl_subset.")

from titan_pbl.config import load_config  # noqa: E402
from titan_pbl.config.physical_constants import STEFAN_BOLTZMANN  # noqa: E402
from titan_pbl.core import build_grid, initialize_state  # noqa: E402
from titan_pbl.schemes.radiation import GrayRadiation  # noqa: E402

OUT = REPO_ROOT / "build" / "source_a"
QUICKLOOK = REPO_ROOT / "build" / "r1" / "quicklook"
BOUNDARIES = ("full_absorption", "kirchhoff")
CONTRASTS_FIG2 = (0.0, 1.0, 2.0, 3.0, 4.0)
PROFILE_ALTITUDES = (5.4, 171.0, 541.0, 1382.0, 1652.0)
T_ISO = 93.65


def load_figure_module():
    os.environ["PAPER_A_OUT"] = tempfile.mkdtemp(prefix="r1_t5_")
    spec = importlib.util.spec_from_file_location(
        "make_paper_A_figures_r1", REPO_ROOT / "code" / "analysis" / "source_a_core.py"
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def rel(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT)).replace("\\", "/")


def entry(key, value, units, boundary, variant, source_file, source_function, sig, why):
    return {
        "key": key,
        "value": None if value is None or (isinstance(value, float) and not np.isfinite(value)) else value,
        "units": units,
        "lw_boundary": boundary,
        "sw_variant": variant,
        "source_file": source_file,
        "source_function": source_function,
        "print_sig_figs": sig,
        "print_sig_figs_reason": why,
    }


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    QUICKLOOK.mkdir(parents=True, exist_ok=True)
    mod = load_figure_module()
    sol = float(mod.SOL_LENGTH)

    thresholds_rows, mixing_rows, profile_rows = [], [], []
    ctxs = {}
    results: dict = {}
    for boundary in BOUNDARIES:
        mod.LW_BOUNDARY_OVERRIDES.clear()
        mod.LW_BOUNDARY_OVERRIDES["radiation.surface_lw_boundary"] = boundary
        ctx = mod.load_framework()
        ctxs[boundary] = ctx
        results[boundary] = {}
        for variant in mod.SW_VARIANTS:
            thr = mod.compute_thresholds_variant(ctx, variant)
            mix = mod.mixing_times_variant(ctx, thr, variant)
            layer0 = float(mod.heating_profile_variant(ctx, 0.0, variant)[0])
            results[boundary][variant] = {"thresholds": thr, "mixing": mix, "layer0_dT0": layer0}
            thresholds_rows.append({
                "lw_boundary": boundary, "sw_variant": variant,
                "analytic_eq4_K": thr["analytic"],
                "heating_threshold_layer0_K": thr["heating"],
                "destabilisation_threshold_0_4km_K": thr["destabilisation"],
                "layer0_heating_at_dT0_K_per_sol": layer0,
            })
            for row in mix:
                mixing_rows.append({"lw_boundary": boundary, "sw_variant": variant, **row})
            for contrast in CONTRASTS_FIG2:
                rates = mod.heating_profile_variant(ctx, contrast, variant)
                for k, z in enumerate(ctx["altitude"]):
                    if z <= 4000.0:
                        profile_rows.append({
                            "lw_boundary": boundary, "sw_variant": variant, "contrast_K": contrast,
                            "layer": k, "altitude_m": float(z), "heating_K_per_sol": float(rates[k]),
                        })
        # Variante publicada: confirmar que coincide con compute_thresholds.
        legacy = mod.compute_thresholds(ctx)
        results[boundary]["legacy_compute_thresholds"] = legacy

    def write_csv(path, rows):
        with path.open("w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    write_csv(OUT / "thresholds_sw_variants.csv", thresholds_rows)
    write_csv(OUT / "mixing_times_sw_variants.csv", mixing_rows)
    write_csv(OUT / "fig2_profiles_sw_variants.csv", profile_rows)

    # ---- D4: contrastes del archivo ------------------------------------
    contrasts = mod.archive_contrast_maxima()
    rows = []
    for label, info in contrasts.items():
        for s, v in info["per_sol_max_K"].items():
            rows.append({"run": label, "sol": s, "max_contrast_K": v})
        rows.append({"run": label, "sol": "2-5", "max_contrast_K": info["max_sols_2_to_5_K"]})
    write_csv(OUT / "archive_contrasts_D4.csv", rows)

    # ---- Comprobaciones de consistencia ---------------------------------
    cfg = load_config("huygens_baseline", overrides={
        "radiation.surface_lw_boundary": "kirchhoff",
        "atmosphere.initialization.profile_source": "isothermal",
        "atmosphere.initialization.isothermal_temperature": T_ISO,
        "radiation.ir_top_flux": STEFAN_BOLTZMANN * T_ISO**4,
        "radiation.solar_constant": 0.0,
    })
    st = initialize_state(cfg, build_grid(cfg))
    sc = GrayRadiation(cfg)
    sc.initialize(cfg, st)
    iso_max = float(np.max(np.abs(sc.heating_rate(st))))

    ctx_a, ctx_k = ctxs["full_absorption"], ctxs["kirchhoff"]
    lw_a = mod.heating_components(ctx_a, 0.0)["lw"][0]
    lw_k = mod.heating_components(ctx_k, 0.0)["lw"][0]
    state, scheme, config = ctx_a["state"], ctx_a["scheme"], ctx_a["config"]
    down, _ = scheme.longwave_fluxes(state)
    tau = scheme.infrared_optical_depth(state)
    d_tau0 = float(tau[0] - tau[1])
    dp0 = float(state.grid.atmosphere.delta_pressure[0])
    eps = float(config.radiation.surface_emissivity)
    expected = (
        config.planet.gravity * (1.0 - eps) * float(down[0]) * (1.0 - np.exp(-d_tau0))
        / (config.atmosphere.gas.specific_heat_cp * dp0) * sol
    )
    shift = float(lw_k - lw_a)
    checks = {
        "isothermal_kirchhoff_max_abs_heating_K_s": iso_max,
        "isothermal_kirchhoff_pass": bool(iso_max <= 1.0e-12),
        "layer0_lw_shift_K_per_sol": shift,
        "layer0_lw_shift_expected_K_per_sol": expected,
        "layer0_lw_shift_relative_error": abs(shift - expected) / abs(expected),
        "layer0_lw_shift_pass": bool(abs(shift - expected) <= 1.0e-6 * abs(expected)),
        "noon_variant_equals_legacy": {
            b: {
                "heating": results[b]["noon"]["thresholds"]["heating"] == results[b]["legacy_compute_thresholds"]["heating"],
                "destabilisation": results[b]["noon"]["thresholds"]["destabilisation"]
                == results[b]["legacy_compute_thresholds"]["destabilisation"],
            }
            for b in BOUNDARIES
        },
    }
    (OUT / "consistency_checks.json").write_text(json.dumps(checks, indent=2), encoding="utf-8")

    # ---- Fuente unica de numeros ----------------------------------------
    fig = "code/analysis/source_a_core.py"
    script = rel(Path(__file__))
    nums = []
    sig_thr = (2, "eps se da con dos cifras (0.95); en T4 un cambio de 0.05 en eps mueve el umbral del orden de 0.5 K, de modo que la tercera cifra no es significativa")
    sig_rate = (2, "tasas derivadas de flujos con parametros de 1-2 cifras (eps, f, I_IR); el valor relativo entre variantes es lo informativo")
    sig_mix = (2, "cociente de dTheta (0.652 K, 3 cifras) entre una tasa de 2 cifras")
    for b in BOUNDARIES:
        tdir = OUT / b
        for v in mod.SW_VARIANTS:
            r = results[b][v]
            nums.append(entry(f"threshold_heating_layer0_{b}_{v}", r["thresholds"]["heating"], "K", b, v,
                              rel(OUT / "thresholds_sw_variants.csv"), f"{fig}::compute_thresholds_variant", *sig_thr))
            nums.append(entry(f"threshold_destabilisation_0_4km_{b}_{v}", r["thresholds"]["destabilisation"], "K", b, v,
                              rel(OUT / "thresholds_sw_variants.csv"), f"{fig}::compute_thresholds_variant", *sig_thr))
            nums.append(entry(f"layer0_heating_dT0_{b}_{v}", r["layer0_dT0"], "K per Titan sol", b, v,
                              rel(OUT / "thresholds_sw_variants.csv"), f"{fig}::heating_profile_variant", *sig_rate))
            for row in r["mixing"][1:]:
                c = row["contrast_K"]
                nums.append(entry(f"destabilisation_rate_{c:g}K_{b}_{v}", row["destabilisation_K_per_sol"], "K per Titan sol",
                                  b, v, rel(OUT / "mixing_times_sw_variants.csv"), f"{fig}::mixing_times_variant", *sig_rate))
                nums.append(entry(f"mixing_time_{c:g}K_{b}_{v}", row["mixing_time_sols"], "Titan sols", b, v,
                                  rel(OUT / "mixing_times_sw_variants.csv"), f"{fig}::mixing_times_variant", *sig_mix))
        # Perfil de la Fig. 2 a dT = 0 (variante publicada) en las alturas que cita el texto.
        ctx = ctxs[b]
        prof = mod.heating_profile(ctx, 0.0)
        for z in PROFILE_ALTITUDES:
            k = int(np.argmin(np.abs(ctx["altitude"] - z)))
            nums.append(entry(f"profile_noon_dT0_{z:g}m_{b}", float(prof[k]), "K per Titan sol", b, "noon",
                              rel(OUT / "fig2_profiles_sw_variants.csv"), f"{fig}::heating_profile", *sig_rate))
        # Tabla 2 (barrido) regenerada por el script con --lw-boundary.
        import pandas as pd

        t2 = pd.read_csv(tdir / "table2_parameter_sweep.csv")
        ref = t2.loc[np.isclose(t2["I_IR_W_m2"], 2.0) & np.isclose(t2["tau_IR"], 10.0) & np.isclose(t2["f"], 0.15)].iloc[0]
        nums.append(entry(f"sweep_max_heating_below_2p5km_{b}", float(t2["max_heating_K_per_sol"].max()), "K per Titan sol",
                          b, "noon", rel(tdir / "table2_parameter_sweep.csv"), f"{fig}::table2_parameter_sweep", *sig_rate))
        nums.append(entry(f"sweep_reference_layer0_{b}", float(ref["surface_heating_K_per_sol"]), "K per Titan sol",
                          b, "noon", rel(tdir / "table2_parameter_sweep.csv"), f"{fig}::table2_parameter_sweep", *sig_rate))
        nums.append(entry(f"sweep_all_gradients_negative_{b}", bool((t2["gradient_K_per_sol"] < 0).all()), "bool",
                          b, "noon", rel(tdir / "table2_parameter_sweep.csv"), f"{fig}::table2_parameter_sweep", None,
                          "booleano"))
        nums.append(entry(f"sweep_any_layer_heats_below_2p5km_{b}", bool((t2["max_heating_K_per_sol"] > 0).any()), "bool",
                          b, "noon", rel(tdir / "table2_parameter_sweep.csv"), f"{fig}::table2_parameter_sweep", None,
                          "booleano"))

    # Umbral analitico de la Ec. 4 (O1) y su valor con Kirchhoff (cero por derivacion).
    ctx = ctxs["full_absorption"]
    eps = float(ctx["config"].radiation.surface_emissivity)
    nums.append(entry("analytic_threshold_eq4_full_absorption", ctx["air_temperature"] * (eps ** -0.25 - 1.0), "K",
                      "full_absorption", "n/a", fig, f"{fig}::compute_thresholds (Ta*(eps^-1/4 - 1))", 2,
                      "Ta (2 decimales de HASI interpolado) y eps de 2 cifras"))
    nums.append(entry("analytic_threshold_kirchhoff", 0.0, "K", "kirchhoff", "n/a", "manuscript, Section 2.3 (Eqs. 5-6)",
                      "derivacion: con F_up = eps sigma Ts^4 + (1-eps) F_down y F_down ~ sigma Ta^4, dE/dt ~ dtau eps sigma (Ts^4 - Ta^4)",
                      None, "cero exacto por derivacion"))

    # Geometria y estado inicial (independientes de la frontera).
    state = ctx["state"]
    z_int, z_mid = state.geopotential_heights()
    tau = ctx["scheme"].infrared_optical_depth(state)
    theta = state.potential_temperature()
    top = int(np.argmin(np.abs(z_mid - 4000.0)))
    geometry = [
        ("Ta_layer0", ctx["air_temperature"], "K", 4, "HASI interpolado; Ta se conoce a 0.01 K"),
        ("z_mid_layer0", float(z_mid[0]), "m", 2, "malla inferida (H-01); 2 cifras bastan"),
        ("thickness_layer0", float(z_int[1] - z_int[0]), "m", 3, "idem"),
        ("d_tau_layer0", float(tau[0] - tau[1]), "-", 2, "tau_IR = 10 y f = 0.15 dados con 2 cifras"),
        ("delta_theta_0_4km", float(theta[top] - theta[0]), "K", 3, "valor usado en el texto (0.652); 3 cifras de HASI"),
    ]
    for key, value, units, sig, why in geometry:
        nums.append(entry(key, value, units, "n/a", "n/a", fig, f"{fig}::load_framework / state", sig, why))

    # Ec. 6, escala radiativa (T = 93 K del texto).
    cfg0 = ctx["config"]
    tau_rad = (cfg0.atmosphere.gas.specific_heat_cp * cfg0.atmosphere.vertical_grid.surface_pressure
               / (4.0 * cfg0.planet.gravity * STEFAN_BOLTZMANN * 93.0**3))
    nums += [
        entry("tau_rad_eq6", tau_rad, "s", "n/a", "n/a", script, "Ec. 6: cp p / (4 g sigma T^3), T = 93 K", 2,
              "T = 93 K del texto (2 cifras)"),
        entry("tau_rad_eq6_sols", tau_rad / cfg0.planet.sol_length, "Titan sols", "n/a", "n/a", script,
              "Ec. 6 / planet.sol_length", 2, "idem; v0 imprimio 469 truncado (H-29)"),
        entry("tau_rad_eq6_years", tau_rad / cfg0.planet.orbital_period, "Titan years", "n/a", "n/a", script,
              "Ec. 6 / planet.orbital_period", 2, "idem; v0 imprimio 0.69 truncado (H-29)"),
    ]

    # D4 y Tabla 4 (archivo; independientes de la frontera).
    t4 = __import__("pandas").read_csv(OUT / "full_absorption" / "table4_contrasts.csv")
    archive = [
        ("archive_contrast_max_sols2to5_baseline", contrasts["baseline"]["max_sols_2_to_5_K"], "K", 2,
         "lectura de TSK - T_PHY del archivo; D4"),
        ("archive_contrast_max_sols2to5_suppressed", contrasts["suppressed"]["max_sols_2_to_5_K"], "K", 2,
         "idem"),
        ("archive_contrast_band_low", contrasts["baseline"]["max_sols_2_to_5_K"], "K", 2, "banda D4 = [baseline, suprimida]"),
        ("archive_contrast_band_high", contrasts["suppressed"]["max_sols_2_to_5_K"], "K", 2, "idem"),
        ("archive_hfx_peak_baseline_sol3", float(t4["HFX_baseline_peak_W_m2"].iloc[1]), "W m-2", 3, "lectura del archivo"),
        ("archive_hfx_peak_suppressed_sol3", float(t4["HFX_suppressed_peak_W_m2"].iloc[1]), "W m-2", 3, "idem"),
        ("archive_hfx_reduction_sol3", (1.0 - float(t4["HFX_ratio"].iloc[1])) * 100.0, "%", 2, "cociente de dos lecturas"),
        ("archive_pbl_drift_suppressed_sols2to5", (float(t4["PBL_suppressed_max_m"].iloc[3]) / float(t4["PBL_suppressed_max_m"].iloc[0]) - 1.0) * 100.0,
         "%", 2, "cociente de dos lecturas"),
    ]
    for key, value, units, sig, why in archive:
        nums.append(entry(key, value, units, "n/a", "n/a", rel(OUT / "archive_contrasts_D4.csv") if "contrast" in key
                          else rel(OUT / "full_absorption" / "table4_contrasts.csv"),
                          f"{fig}::archive_contrast_maxima" if "contrast" in key else f"{fig}::table4_contrasts", sig, why))

    (OUT / "numbers_r1.json").write_text(
        json.dumps({"generated_by": script, "sol_length_s": sol, "numbers": nums}, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )

    # ---- Graficas de inspeccion -----------------------------------------
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig_, axes = plt.subplots(2, 3, figsize=(11, 7), sharey=True)
    for i, b in enumerate(BOUNDARIES):
        for j, v in enumerate(mod.SW_VARIANTS):
            ax = axes[i, j]
            for c in CONTRASTS_FIG2:
                sel = [r for r in profile_rows if r["lw_boundary"] == b and r["sw_variant"] == v and r["contrast_K"] == c]
                ax.plot([r["heating_K_per_sol"] for r in sel], [r["altitude_m"] for r in sel], label=f"dT={c:g} K")
            ax.axvline(0, color="k", lw=0.5)
            ax.set_title(f"{b} | {v}", fontsize=8)
            ax.set_xlabel("K per Titan sol", fontsize=8)
    axes[0, 0].legend(fontsize=7)
    axes[0, 0].set_ylabel("altitude [m]")
    axes[1, 0].set_ylabel("altitude [m]")
    fig_.suptitle("R1 quicklook: Fig. 2 profiles by LW boundary and SW variant (not a final figure)", fontsize=9)
    fig_.savefig(QUICKLOOK / "fig2_profiles_variants.png", dpi=120)
    plt.close(fig_)

    fig_, ax = plt.subplots(figsize=(6, 4))
    for b, style in zip(BOUNDARIES, ("-", "--")):
        for v in mod.SW_VARIANTS:
            sel = [r for r in mixing_rows if r["lw_boundary"] == b and r["sw_variant"] == v and r["contrast_K"] in (3.0, 3.5, 4.5)]
            ax.plot([r["contrast_K"] for r in sel], [r["mixing_time_sols"] for r in sel], style, marker="o", label=f"{b} | {v}")
    ax.axhspan(1, 5, color="green", alpha=0.15)
    ax.set_yscale("log")
    ax.set_xlabel("contrast [K]")
    ax.set_ylabel("mixing time [Titan sols]")
    ax.legend(fontsize=7)
    ax.set_title("R1 quicklook: mixing times (not a final figure)", fontsize=9)
    fig_.savefig(QUICKLOOK / "mixing_times_variants.png", dpi=120)
    plt.close(fig_)

    print(json.dumps(checks, indent=1))
    for r in thresholds_rows:
        print(r)
    return 0 if checks["isothermal_kirchhoff_pass"] and checks["layer0_lw_shift_pass"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
