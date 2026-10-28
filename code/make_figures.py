#!/usr/bin/env python
"""
make_figures.py -- figures of the revised Paper A.

Drawn ONLY from the single source (numbers/numbers_paper_A.json) and from the
CSV/JSON outputs of the analyses under build/ (including build/figdata/,
written by code/analysis/make_figure_data.py). No result value is written
here; only style parameters (figstyle.py).

Reglas:
- Estilo comun (figstyle.py): Arial, anchos Elsevier 90/140/190 mm,
  texto >= 7 pt al tamano impreso, Okabe-Ito, Kirchhoff-consistent continua y
  full absorption discontinua en todas las figuras.
- PDF vectorial (Arial TrueType incrustada) + PNG 1000 dpi; texto en ingles
  estadounidense con los terminos del manuscrito (figstyle.TERMS).
- Nombres finales del manuscrito: fig01 ... fig07, figS1, figS2 (FINAL_NAMES).
- Las claves "superseded" no se pueden usar (el script falla).
- El graphical abstract no se genera aqui: es un archivo del autor
  (figures/graphical_abstract.*).

Salidas: figures/ y figures/figures_manifest.csv (fuentes, claves y SHA-256).
"""

from __future__ import annotations

import csv
import json
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import figstyle as S  # noqa: E402

import matplotlib.patches as mpatches  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

FIG_DIR = REPO_ROOT / "figures"
FIGDATA = REPO_ROOT / "build" / "figdata"
NUMBERS = REPO_ROOT / "numbers" / "numbers_paper_A.json"
CONTRASTS = REPO_ROOT / "build" / "numbers" / "literature_contrasts.csv"
R1 = REPO_ROOT / "build" / "source_a"
R2_OUT = REPO_ROOT / "build" / "r2"
R3_OUT = REPO_ROOT / "build" / "r3"
BOUNDARIES = ("kirchhoff", "full_absorption")

# Nombre interno de cada funcion -> nombre final en el manuscrito.
FINAL_NAMES = {
    "fig_boundary_schematic": "fig01_boundary_schematic",
    "fig_heating_profiles": "fig02_heating_profiles",
    "fig_mixing_time": "fig03_mixing_time",
    "fig_flux_comparison": "fig04_flux_comparison",
    "fig_archive_contrasts": "fig05_archive_contrasts",
    "fig_stability_budget": "fig06_stability_budget",
    "fig_h33_flux_scaling": "fig07_myj_flux_scaling",
    "figS_sw_variants": "figS1_sw_variants",
    "figS_sensitivity": "figS2_sensitivity",
}


def rel(p: Path) -> str:
    return str(Path(p).resolve().relative_to(REPO_ROOT)).replace("\\", "/")


class Sources:
    """Reads numbers and files and records what each figure used."""

    def __init__(self):
        data = json.loads(NUMBERS.read_text(encoding="utf-8"))
        self.numbers = {n["key"]: n for n in data["numbers"]}
        self.reset()

    def reset(self):
        self.keys, self.files = [], {rel(NUMBERS)}

    def n(self, key: str) -> float:
        e = self.numbers[key]
        if e.get("status") == "superseded":
            raise SystemExit(f"clave superseded usada en una figura: {key}")
        self.keys.append(key)
        return e["value"]

    def p(self, key: str) -> str:
        self.n(key)
        return self.numbers[key]["print"]

    def csv(self, path: Path) -> list[dict]:
        self.files.add(rel(path))
        with path.open(newline="", encoding="utf-8") as fh:
            return list(csv.DictReader(fh))

    def json(self, path: Path):
        self.files.add(rel(path))
        return json.loads(path.read_text(encoding="utf-8"))


SRC = Sources()
MANIFEST: list[dict] = []


def finish(fig, name: str, formats=("pdf", "png")) -> None:
    paths = S.save(fig, FIG_DIR, FINAL_NAMES[name], formats)
    for p in paths:
        MANIFEST.append({"file": rel(p), "data_sources": " | ".join(sorted(SRC.files)),
                         "keys_used": " | ".join(sorted(set(SRC.keys))), "sha256": S.sha256(p)})
    SRC.reset()


def f(x) -> float:
    return float(x)


# ---------------------------------------------------------------------------
def arrow(ax, x, y0, y1, color, lw=1.6):
    ax.annotate("", xy=(x, y1), xytext=(x, y0),
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw, mutation_scale=9, shrinkA=0, shrinkB=0))


def schematic_panel(ax, kind: str, labels: dict | None = None) -> None:
    oi = S.OKABE_ITO
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.add_patch(mpatches.Rectangle((0.3, 5.6), 9.4, 1.6, fc="#E8F1F8", ec="k", lw=0.6))
    ax.add_patch(mpatches.Rectangle((0.3, 0.4), 9.4, 1.3, fc="#E3D3B5", ec="k", lw=0.6))
    lab = labels or {}
    ax.text(5.0, 6.4, lab.get("layer", "Layer, $T_a$"), ha="center", va="center")
    ax.text(5.0, 1.05, lab.get("surface", "Surface, $T_s$, emissivity $\\varepsilon$"), ha="center", va="center")
    if kind != "ga":
        # Emision de la capa hacia arriba y hacia abajo
        arrow(ax, 8.6, 7.2, 9.6, oi["black"])
        ax.text(8.4, 8.5, "$\\sigma T_a^4$", ha="right", va="center")
        arrow(ax, 8.6, 5.6, 1.7, oi["black"])
        ax.text(8.4, 3.0, "$\\sigma T_a^4$", ha="right", va="center")
    # F_down llega a la superficie
    x_dn = 6.6 if kind != "ga" else 9.0
    arrow(ax, x_dn, 5.6, 1.7, oi["blue"])
    ax.text(x_dn - 0.15, 4.6, "$F^{\\downarrow}$", ha="right", va="center", color=oi["blue"])
    if kind == "full":
        arrow(ax, 1.6, 1.7, 5.6, oi["vermillion"])
        ax.text(1.75, 3.7, "$F^{\\uparrow}=\\varepsilon\\sigma T_s^4$", ha="left", va="center", color=oi["vermillion"])
        ax.text(5.0, 0.0, "Surface absorbs $F^{\\downarrow}$ in full", ha="center", va="bottom", color=oi["blue"])
        ax.text(5.0, 9.0, "Isothermal residual:\n$(1-\\varepsilon)\\sigma T^4$ lost by the air", ha="center",
                va="center", color=oi["purple"])
    else:
        arrow(ax, 1.3, 1.7, 5.6, oi["vermillion"])
        arrow(ax, 3.5, 1.7, 5.6, oi["green"])
        if kind == "ga":
            ax.text(1.45, 3.65, lab["emitted"], ha="left", va="center", color=oi["vermillion"])
            ax.text(3.65, 2.6, lab["reflected"], ha="left", va="center", color=oi["green"])
            ax.text(5.0, 0.0, lab["absorbed"], ha="center", va="bottom", color=oi["blue"])
        else:
            ax.text(1.45, 4.4, "Emitted\n$\\varepsilon\\sigma T_s^4$", ha="left", va="center", color=oi["vermillion"])
            ax.text(3.65, 2.8, "Reflected\n$(1-\\varepsilon)F^{\\downarrow}$", ha="left", va="center", color=oi["green"])
            ax.text(5.0, 0.0, "Surface absorbs $\\varepsilon F^{\\downarrow}$", ha="center", va="bottom", color=oi["blue"])


def fig_boundary_schematic() -> None:
    fig, axes = S.figure("onehalf", 72, ncols=2)
    schematic_panel(axes[0], "full")
    schematic_panel(axes[1], "kirchhoff")
    S.panel_title(axes[0], "a", f"{S.TERMS['full_absorption']} (Chatain et al.; Hartwick et al.)")
    S.panel_title(axes[1], "b", S.TERMS["kirchhoff"])
    fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.02, wspace=0.05)
    finish(fig, "fig_boundary_schematic")


def profiles(rows, boundary, variant):
    out = {}
    for r in rows:
        if r["lw_boundary"] == boundary and r["sw_variant"] == variant:
            out.setdefault(f(r["contrast_K"]), []).append((f(r["altitude_m"]), f(r["heating_K_per_sol"])))
    return {c: np.array(sorted(v)) for c, v in sorted(out.items())}


def draw_profiles(ax, rows, boundary, variant, legend=True):
    for (c, arr), col in zip(profiles(rows, boundary, variant).items(), S.CONTRAST_COLORS):
        ax.plot(arr[:, 1], arr[:, 0] / 1000.0, ls=S.BOUNDARY_STYLE[boundary]["ls"], color=col,
                label=f"ΔT = {c:g} K")
    ax.axvline(0, color="k", lw=0.4)
    ax.set_ylim(0, 4.0)
    if legend:
        ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.92))


def fig_heating_profiles() -> None:
    rows = SRC.csv(R1 / "fig2_profiles_sw_variants.csv")
    fig, axes = S.figure("onehalf", 70, ncols=2, sharey=True)
    for ax, b, lab in zip(axes, ("full_absorption", "kirchhoff"), ("a", "b")):
        draw_profiles(ax, rows, b, "lw")
        S.panel_title(ax, lab, f"{S.BOUNDARY_STYLE[b]['label']}, {S.TERMS['longwave_only']}")
        ax.set_xlabel(f"Radiative heating rate ({S.TERMS['rate_unit']})")
    axes[0].set_ylabel("Altitude (km)")
    fig.tight_layout()
    finish(fig, "fig_heating_profiles")


def mixing_curves():
    rows = SRC.csv(FIGDATA / "mixing_time_curves.csv")
    out = {}
    for b in BOUNDARIES:
        rr = [r for r in rows if r["lw_boundary"] == b]
        x = np.array([f(r["contrast_K"]) for r in rr])
        y = np.array([f(r["mixing_time_sols"]) for r in rr])
        y[~np.isfinite(y)] = np.nan
        out[b] = (x, y)
    return out


def verified_contrasts():
    rows = SRC.csv(CONTRASTS)
    return [r for r in rows if r["surface_to_air_contrast"] == "yes"]


def fig_mixing_time() -> None:
    curves = mixing_curves()
    win = SRC.json(FIGDATA / "integration_window.json")
    dth, sig = SRC.n("delta_theta_0_4km"), SRC.n("hasi_sigma_dtheta_independent_brief")
    hatch_lw_default = plt.rcParams["hatch.linewidth"]
    plt.rcParams["hatch.linewidth"] = 0.35  # se restaura al final; ninguna otra figura usa achurado
    fig, ax = S.figure("single", 115)
    ax.axhspan(win["analysis_window_sols"], win["integration_length_sols"], color=S.OKABE_ITO["green"], alpha=0.25,
               lw=0, label="Integration window (4–5 Titan sols)")
    # Etiquetas del glosario (figstyle.TERMS / BOUNDARY_STYLE).
    curve_label = {b: S.BOUNDARY_STYLE[b]["label"] for b in BOUNDARIES}
    hasi_label = {b: f"HASI ±σ, {S.BOUNDARY_STYLE[b]['label_mid']}" for b in BOUNDARIES}
    for b in BOUNDARIES:
        x, y = curves[b]
        st = S.BOUNDARY_STYLE[b]
        # La entrada de leyenda es el propio poligono: mismo color y alfa que la envolvente.
        ax.fill_between(x, y * (dth - sig) / dth, y * (dth + sig) / dth, color=st["color"], alpha=0.15, lw=0,
                        label=hasi_label[b])
        ax.plot(x, y, ls=st["ls"], color=st["color"], label=curve_label[b])
        for n, mk in ((4, "o"), (5, "s")):
            ax.plot(SRC.n(f"dTreq_{n}sol_{b}_lw"), n, mk, color=st["color"], ms=3.5, mfc="white", mew=0.9)
    ax.plot([], [], "o", color="0.3", mfc="white", ms=3.5, label="ΔT$_{req}$(4)")
    ax.plot([], [], "s", color="0.3", mfc="white", ms=3.5, label="ΔT$_{req}$(5)")
    contrast_rows = SRC.csv(CONTRASTS)
    groups = {}
    for r in contrast_rows:
        if r["surface_to_air_contrast"] == "yes":
            groups.setdefault(float(re.search(r"([0-9.]+) m", r["Ta_height"]).group(1)), []).append(f(r["value_K"]))
    # Un color Okabe-Ito por altura de Ta, distinto de curvas (azul, bermellon), envolventes y ventana (verde).
    band_colors = (S.OKABE_ITO["yellow"], S.OKABE_ITO["purple"], S.OKABE_ITO["orange"])
    spans = {h: (min(vals) * 0.98, max(vals) * 1.02) for h, vals in groups.items()}
    for (h, (lo, hi)), col in zip(sorted(spans.items()), band_colors):
        # Una banda que se solapa con otra mas ancha se dibuja opaca y encima, para que conserve su color de leyenda.
        inside = any(o != h and lo < ohi and olo < hi and (hi / lo) < (ohi / olo) for o, (olo, ohi) in spans.items())
        ax.axvspan(lo, hi, color=col, alpha=1.0 if inside else 0.55, lw=0, zorder=1.5 if inside else 1.0,
                   label=f"Verified contrasts, T$_a$ at {h:g} m")
    bulk = [f(r["value_K"]) for r in contrast_rows if r["surface_to_air_contrast"] == "bulk"]
    ax.axvspan(min(bulk), max(bulk), facecolor="none", edgecolor="0.35", hatch="///", lw=0.5, zorder=1.2,
               label="Bulk estimates (Lorenz 2021; Han and Lora 2025)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.1, 20)
    ax.set_ylim(1, 1000)
    ax.set_xlabel("Surface-to-air contrast ΔT (K)")
    ax.set_ylabel("Time to mix the lowest 4 km (Titan sols)")
    handles, labels = ax.get_legend_handles_labels()
    order = ["Integration window (4–5 Titan sols)", curve_label["kirchhoff"], curve_label["full_absorption"],
             hasi_label["kirchhoff"], hasi_label["full_absorption"], "ΔT$_{req}$(4)", "ΔT$_{req}$(5)"]
    order += [l for l in labels if l.startswith("Verified")] + ["Bulk estimates (Lorenz 2021; Han and Lora 2025)"]
    by_label = dict(zip(labels, handles))
    ax.legend([by_label[l] for l in order], order, loc="upper center", bbox_to_anchor=(0.42, -0.14), ncol=1,
              handlelength=2.0)
    fig.tight_layout()
    finish(fig, "fig_mixing_time")
    plt.rcParams["hatch.linewidth"] = hatch_lw_default


def fig_flux_comparison() -> None:
    rows = SRC.csv(R3_OUT / "t3_flux_profiles.csv")
    SRC.n("flux_comparison_max_rms_rel_F_down_pct")
    SRC.n("flux_comparison_0p1flhc_max_rms_rel_F_down_pct")
    fig, axes = S.figure("double", 92, ncols=3, sharey=True)
    from matplotlib.lines import Line2D

    # Codificacion comun a los tres paneles: color = corrida (RUN_STYLE); F_down
    # linea fina y circulos abiertos; F_up linea gruesa translucida y circulos
    # llenos. Los flujos calculados usan la frontera full absorption (discontinua).
    enc = {"F_down": {"sym": "↓", "lw": 1.0, "alpha": 1.0, "filled": False},
           "F_up": {"sym": "↑", "lw": 2.2, "alpha": 0.45, "filled": True}}
    for ax, run, lab in ((axes[0], "TMN91", "a"), (axes[1], "0.1flhc", "b")):
        rr = [r for r in rows if r["run"] == run]
        z = np.array([f(r["z_interface_m"]) for r in rr]) / 1000.0
        col = S.RUN_STYLE[run]["color"]
        for fl, e in enc.items():
            ax.plot([f(r[f"{fl}_archive"]) for r in rr], z, "o", ms=2, color=col,
                    mfc=col if e["filled"] else "none", mew=0.6)
            ax.plot([f(r[f"{fl}_titan_pbl"]) for r in rr], z, ls=S.BOUNDARY_STYLE["full_absorption"]["ls"],
                    color=col, lw=e["lw"], alpha=e["alpha"])
        S.panel_title(ax, lab, S.RUN_STYLE[run]["label"].replace(" (", "\n("))
        ax.set_xlabel("Longwave flux (W m$^{-2}$)")
    grey = "0.3"
    proxies = []
    for fl, e in enc.items():
        proxies += [Line2D([], [], ls="none", marker="o", ms=2.5, color=grey, mfc=grey if e["filled"] else "none",
                           mew=0.6, label=f"{S.TERMS['archive']}, F{e['sym']}"),
                    Line2D([], [], ls=S.BOUNDARY_STYLE["full_absorption"]["ls"], color=grey, lw=e["lw"],
                           alpha=e["alpha"], label=f"{S.TERMS['this_work']}, F{e['sym']}")]
    axes[0].legend(handles=proxies, loc="lower left")
    axes[0].set_ylabel("Altitude (km)")
    for run in ("TMN91", "baseline5sol", "0.1flhc"):
        rr = [r for r in rows if r["run"] == run]
        z = np.array([f(r["z_interface_m"]) for r in rr]) / 1000.0
        col = S.RUN_STYLE[run]["color"]
        for fl, e in enc.items():
            axes[2].plot(100 * np.array([f(r[f"mean_rel_diff_{fl}"]) for r in rr]), z,
                         ls=S.BOUNDARY_STYLE["full_absorption"]["ls"], color=col, lw=e["lw"], alpha=e["alpha"])
    axes[2].axvline(0, color="k", lw=0.4)
    S.panel_title(axes[2], "c", "Relative difference")
    axes[2].set_xlabel(f"{S.TERMS['this_work']} − archive (%)")
    run_handles = [Line2D([], [], color=S.RUN_STYLE[r]["color"], lw=1.4, label=S.RUN_STYLE[r]["label"])
                   for r in ("TMN91", "baseline5sol", "0.1flhc")]
    fig.legend(handles=run_handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.08, 1, 1))
    finish(fig, "fig_flux_comparison")


def fig_archive_contrasts() -> None:
    rows = SRC.csv(FIGDATA / "archive_contrast_timeseries.csv")
    fig, ax = S.figure("onehalf", 82)
    zta = []
    for run in ("TMN91", "0.1flhc"):
        rr = [r for r in rows if r["run"] == run]
        ax.plot([f(r["tsol"]) for r in rr], [f(r["contrast_K"]) for r in rr], color=S.RUN_STYLE[run]["color"],
                lw=0.9, ls=":", label=S.RUN_STYLE[run]["label"])
        zta += [f(r["Ta_height_m"]) for r in rr]
    for b in BOUNDARIES:
        st = S.BOUNDARY_STYLE[b]
        ax.axhline(SRC.n(f"threshold_heating_layer0_{b}_lw"), ls=st["ls"], color="0.15", lw=0.9,
                   label=f"{S.TERMS['heating_threshold']}, {st['label_mid']}")
        ax.axhline(SRC.n(f"threshold_destabilisation_0_4km_{b}_lw"), ls=st["ls"], color=S.OKABE_ITO["orange"], lw=0.9,
                   label=f"{S.TERMS['destabilization_threshold']}, {st['label_mid']}")
    ax.set_xlim(2, 6)
    ax.set_ylim(bottom=min(-0.1, ax.get_ylim()[0]))
    ax.set_xlabel(S.TERMS["time_axis"])
    ax.set_ylabel(f"T$_s$ − T$_a$ at {np.mean(zta):.1f} m (K)")
    # Leyenda bajo el eje (con los nombres completos de las corridas no cabe
    # dentro): filas = corridas, umbrales Kirchhoff-consistent, full absorption.
    handles, labels = ax.get_legend_handles_labels()
    order = [0, 2, 4, 1, 3, 5]
    fig.legend([handles[i] for i in order], [labels[i] for i in order], loc="lower center", ncol=2,
               bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.17, 1, 1))
    finish(fig, "fig_archive_contrasts")


def fig_stability_budget() -> None:
    rows = SRC.csv(FIGDATA / "stability_budget_timeseries.csv")
    for z in ("1000", "2500", "4000"):
        SRC.n(f"mayor3_Q_S_{z}m")
    fig, axes = S.figure("double", 68, ncols=3, sharex=True)
    oi = S.OKABE_ITO
    for ax, z, lab in zip(axes, (1000.0, 2500.0, 4000.0), ("a", "b", "c")):
        rr = [r for r in rows if f(r["z_ref_m"]) == z]
        t = [f(r["tsol"]) for r in rr]
        # Purpura y no azul para la onda larga: el azul identifica a la corrida TMN91 en las Figs. 4, 5 y 7.
        for k, col, name, ls in (("cum_dS_total_K", "k", "Total", "-"), ("cum_dS_LW_K", oi["purple"], "Longwave", "--"),
                                 ("cum_dS_residual_K", oi["orange"], "Residual (non-longwave)", ":")):
            ax.plot(t, [f(r[k]) for r in rr], color=col, ls=ls, label=name)
        ax.axhline(0, color="k", lw=0.4)
        S.panel_title(ax, lab, f"z$_{{ref}}$ = {z / 1000:g} km")
        ax.set_xlabel(S.TERMS["time_axis"])
    axes[0].set_ylabel(f"Cumulative ΔS, {S.TERMS['suppressed_flux_mid']} − baseline (K)")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    finish(fig, "fig_stability_budget")


def fig_h33_flux_scaling() -> None:
    rows = SRC.csv(FIGDATA / "h33_timeseries.csv")
    for k in ("h33_median_r_TMN91", "h33_median_r_baseline5sol", "h33_median_r_0p1flhc"):
        SRC.n(k)
    fig, axes = S.figure("double", 78, ncols=3)
    for run in ("TMN91", "baseline5sol", "0.1flhc"):
        rr = [r for r in rows if r["run"] == run]
        t = np.array([f(r["tsol"]) for r in rr])
        st = S.RUN_STYLE[run]
        axes[0].plot(t, [f(r["HFX_W_m2"]) for r in rr], color=st["color"], lw=0.9, label=f"HFX, {st['label']}")
        axes[0].plot(t, [f(r["F_MYJ_W_m2"]) for r in rr], color=st["color"], lw=0.9, ls="-.", label=f"F$_{{MYJ}}$, {st['label']}")
        rv = [f(r["r"]) if r["r"] not in ("", None) else np.nan for r in rr]
        axes[1].plot(t, rv, ".", ms=1.5, color=st["color"], label=st["label"])
        axes[2].plot(t, [f(r["cum_resid_HFX_MJ_m2"]) for r in rr], color=st["color"], lw=1.0, label=f"with HFX, {st['label']}")
        axes[2].plot(t, [f(r["cum_resid_FMYJ_MJ_m2"]) for r in rr], color=st["color"], lw=1.0, ls="-.",
                     label=f"with F$_{{MYJ}}$, {st['label']}")
    axes[0].set_ylabel("Surface heat flux (W m$^{-2}$)")
    S.panel_title(axes[0], "a", "Archived HFX and flux entering MYJ")
    axes[1].set_ylabel("HFX/F$_{MYJ}$")
    axes[1].set_yscale("log")
    axes[1].set_ylim(0.05, 2)
    S.panel_title(axes[1], "b", "Ratio, daytime records")
    axes[2].set_ylabel("Cumulative column residual (MJ m$^{-2}$)")
    S.panel_title(axes[2], "c", "Column energy closure")
    for ax in axes:
        ax.set_xlabel(S.TERMS["time_axis"])
        ax.set_xlim(2, 6)
    # Leyenda comun: color = corrida; estilo = flujo (las dos baselines se superponen).
    from matplotlib.lines import Line2D
    runs = [Line2D([], [], color=S.RUN_STYLE[r]["color"], lw=1.2, label=S.RUN_STYLE[r]["label"])
            for r in ("TMN91", "baseline5sol", "0.1flhc")]
    fluxes = [Line2D([], [], color="0.3", lw=1.0, label="Archived HFX"),
              Line2D([], [], color="0.3", lw=1.0, ls="-.", label="F$_{MYJ}$ = c$_p$ρ$_1$AKHS(θ$_{sk}$ − θ$_1$)")]
    # Orden por columnas para que la primera fila sean las tres corridas y la segunda los dos flujos.
    handles = [runs[0], fluxes[0], runs[1], fluxes[1], runs[2]]
    fig.legend(handles=handles, loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    finish(fig, "fig_h33_flux_scaling")


def figS_sw_variants() -> None:
    rows = SRC.csv(R1 / "fig2_profiles_sw_variants.csv")
    fig, axes = S.figure("double", 120, nrows=2, ncols=3, sharex=True, sharey=True)
    names = {"lw": S.TERMS["longwave_only"], "noon": S.TERMS["lw_sw_noon"], "diurnal": S.TERMS["lw_sw_diurnal"]}
    letters = iter("abcdef")
    for i, b in enumerate(("full_absorption", "kirchhoff")):
        for j, v in enumerate(("lw", "noon", "diurnal")):
            ax = axes[i, j]
            draw_profiles(ax, rows, b, v, legend=(i == 0 and j == 0))
            S.panel_title(ax, next(letters), f"{S.BOUNDARY_STYLE[b]['label']},\n{names[v]}")
            if i == 1:
                ax.set_xlabel(f"Radiative heating rate ({S.TERMS['rate_unit']})")
            if j == 0:
                ax.set_ylabel("Altitude (km)")
    fig.tight_layout()
    finish(fig, "figS_sw_variants")


def figS_sensitivity() -> None:
    rows = [r for r in SRC.csv(R2_OUT / "t5_sensitivity_oat.csv") if r["role"] == "barrido declarado"]
    # Nombres de los parametros como en la Tabla 3 del manuscrito.
    names = {"radiation.surface_albedo": "Surface albedo α", "radiation.surface_emissivity": "Surface emissivity ε",
             "radiation.optical_depth_ir": "Infrared optical depth τ$_{IR}$",
             "radiation.ir_weighting_coefficient": "IR weighting coefficient f",
             "radiation.ir_top_flux": "Downward IR at model top I$_{IR}$ (W m$^{-2}$)",
             "S0_factor": "Solar constant factor"}
    factors = list(dict.fromkeys(r["factor"] for r in rows))
    fig, ax = S.figure("onehalf", 85)
    y, ticks, labels = 0, [], []
    for fac in factors:
        for r in [r for r in rows if r["factor"] == fac]:
            col = "k" if r["is_preset"] == "True" else S.OKABE_ITO["blue"]
            ax.plot(f(r["dT_req_4sol_K"]), y, "o", color=col, ms=3.5)
            ticks.append(y)
            labels.append(f"{names[fac]} = {f(r['value']):g}")
            y -= 1
        y -= 0.6
    oi = S.OKABE_ITO
    T = S.TERMS
    for key, col, ls, lab in (("dTreq_4sol_kirchhoff_lw", "k", "-", f"{T['table3_configuration']}, {T['kirchhoff_mid']}"),
                              ("corner_kirchhoff_lw_dTreq4", oi["green"], "-", T["most_favorable_corner"]),
                              ("worstcase_dTreq4", oi["vermillion"], "-", T["corner_minus_sigma"]),
                              ("archive_contrast_max_sols2to5_suppressed", oi["purple"], ":", T["max_verified_contrast"])):
        ax.axvline(SRC.n(key), color=col, ls=ls, lw=0.9, label=lab)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels)
    ax.set_xscale("log")
    ax.set_xlabel(f"ΔT$_{{req}}$(4), {T['kirchhoff_mid']}, {T['longwave_only']} (K)")
    ax.legend(loc="lower left", frameon=True, framealpha=0.95, edgecolor="none")
    fig.tight_layout()
    finish(fig, "figS_sensitivity")


def main(argv: list[str] | None = None) -> int:
    font = S.apply()
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for fn in (fig_boundary_schematic, fig_heating_profiles, fig_mixing_time, fig_flux_comparison,
               fig_archive_contrasts, fig_stability_budget, fig_h33_flux_scaling, figS_sw_variants,
               figS_sensitivity):
        fn()
    with (FIG_DIR / "figures_manifest.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["file", "data_sources", "keys_used", "sha256"], lineterminator="\n")
        w.writeheader()
        w.writerows(MANIFEST)
    (REPO_ROOT / "build").mkdir(exist_ok=True)
    (REPO_ROOT / "build" / "figures_info.json").write_text(json.dumps({
        "font_resolved": font, "font_preference": S.FONT_PREFERENCE, "n_files": len(MANIFEST)}, indent=2),
        encoding="utf-8")
    print("fuente:", font, "| archivos:", len(MANIFEST))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
