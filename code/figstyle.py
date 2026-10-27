"""
figstyle.py -- common style and glossary for the revised Paper A figures.

Only style parameters and figure wording live here: fonts, Elsevier column
widths, minimum text size, output resolution, the Okabe-Ito palette, the
line-style convention (Kirchhoff-consistent solid, full absorption dashed),
one fixed color per archive run, and the glossary of terms shared with the
manuscript text. No result value is written in this module.

Glosario (homologacion de figuras, 2026-10-26): todas las etiquetas
de las figuras toman sus terminos de TERMS, BOUNDARY_STYLE y RUN_STYLE, de
modo que la terminologia vive en un solo lugar y coincide con el manuscrito.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import font_manager  # noqa: E402

# Anchos de columna de Elsevier [mm]
WIDTH_MM = {"single": 90.0, "onehalf": 140.0, "double": 190.0}
MM = 1.0 / 25.4
MIN_FONT_PT = 7.0
BASE_FONT_PT = 7.5

# Fuente: Arial, una de las recomendadas por la guia de artwork de Elsevier
# (Arial/Helvetica, Courier, Symbol, Times). Era la que ya resolvia el estilo
# de R3, porque Liberation Sans no esta instalada; se deja como respaldo.
FONT_PREFERENCE = ["Arial", "Liberation Sans", "DejaVu Sans"]

# Resolucion de los PNG: la guia de Elsevier pide 1000 dpi para line art
# (graficas) y 500 dpi para arte combinado. Los PDF son vectoriales.
PNG_DPI = 1000

# Paleta Okabe-Ito
OKABE_ITO = {
    "black": "#000000", "orange": "#E69F00", "skyblue": "#56B4E9", "green": "#009E73",
    "yellow": "#F0E442", "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7",
}
CONTRAST_COLORS = [OKABE_ITO[c] for c in ("blue", "skyblue", "green", "orange", "vermillion")]

# Terminos del manuscrito, en sentence case. "*_mid" es la forma que va dentro
# de una etiqueta (despues de una coma), con minuscula inicial.
TERMS = {
    "full_absorption": "Full absorption",
    "full_absorption_mid": "full absorption",
    "kirchhoff": "Kirchhoff-consistent",
    "kirchhoff_mid": "Kirchhoff-consistent",
    "baseline": "Baseline",
    "suppressed_flux": "Suppressed flux",
    "suppressed_flux_mid": "suppressed flux",
    "this_work": "This work",
    "archive": "Archive (mtWRF)",
    "heating_threshold": "Heating threshold",
    "destabilization_threshold": "Destabilization threshold",
    "rate_unit": "K per Titan sol",
    "time_axis": "Time (Titan sols)",
    "longwave_only": "longwave only",
    "lw_sw_noon": "longwave + shortwave, local noon",
    "lw_sw_diurnal": "longwave + diurnal-mean shortwave",
    "table3_configuration": "Table 3 configuration",
    "most_favorable_corner": "Most favorable corner",
    "corner_minus_sigma": "Corner with Δθ − σ",
    "max_verified_contrast": "Max. verified contrast",
}

# Convencion de fronteras: Kirchhoff-consistent continua, full absorption discontinua.
BOUNDARY_STYLE = {
    "kirchhoff": {"ls": "-", "label": TERMS["kirchhoff"], "label_mid": TERMS["kirchhoff_mid"],
                  "color": OKABE_ITO["blue"]},
    "full_absorption": {"ls": "--", "label": TERMS["full_absorption"], "label_mid": TERMS["full_absorption_mid"],
                        "color": OKABE_ITO["vermillion"]},
}
# Un color fijo por corrida del archivo, igual en todas las figuras. Los nombres
# son los de las Tablas 6 y 8 del manuscrito (convencion aprobada por el autor).
RUN_STYLE = {
    "TMN91": {"color": OKABE_ITO["blue"], "label": "Baseline (1D_Huygens_TMN91_5sol)"},
    "baseline5sol": {"color": OKABE_ITO["green"], "label": "Baseline (1D_Huygens_baseline_5sol)"},
    "0.1flhc": {"color": OKABE_ITO["vermillion"], "label": "Suppressed flux (1D_Huygens_0.1flhc_5sol)"},
}


def panel_title(ax, letter: str, text: str) -> None:
    """Uniform panel label: "(a) Text", left-aligned above the axes."""
    ax.set_title(f"({letter}) {text}", loc="left")


def resolved_font() -> str:
    for name in FONT_PREFERENCE:
        try:
            path = font_manager.findfont(font_manager.FontProperties(family=name), fallback_to_default=False)
            if path:
                return name
        except ValueError:
            continue
    return "DejaVu Sans"


def apply() -> str:
    font = resolved_font()
    matplotlib.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": [font] + FONT_PREFERENCE,
        "font.size": BASE_FONT_PT, "axes.labelsize": BASE_FONT_PT, "axes.titlesize": BASE_FONT_PT,
        "xtick.labelsize": MIN_FONT_PT, "ytick.labelsize": MIN_FONT_PT, "legend.fontsize": MIN_FONT_PT,
        "mathtext.fontset": "custom", "mathtext.rm": font, "mathtext.it": f"{font}:italic",
        "mathtext.bf": f"{font}:bold", "axes.linewidth": 0.6, "lines.linewidth": 1.0,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.direction": "in", "ytick.direction": "in",
        "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600,
    })
    return font


def figure(width: str = "single", height_mm: float = 70.0, **kw):
    return plt.subplots(figsize=(WIDTH_MM[width] * MM, height_mm * MM), **kw)


def check_min_font(fig) -> float:
    """Smallest text size in the figure [pt] (printed size, figures are made at final width)."""
    sizes = [t.get_fontsize() for t in fig.findobj(matplotlib.text.Text) if t.get_text().strip()]
    smallest = min(sizes) if sizes else MIN_FONT_PT
    if smallest < MIN_FONT_PT - 1e-9:
        raise ValueError(f"texto de {smallest} pt (< {MIN_FONT_PT} pt)")
    return smallest


def save(fig, outdir: Path, name: str, formats=("pdf", "png")) -> list[Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    check_min_font(fig)
    paths = []
    for fmt in formats:
        p = outdir / f"{name}.{fmt}"
        kw = {"metadata": {"CreationDate": None}} if fmt == "pdf" else {}
        fig.savefig(p, format=fmt, dpi=PNG_DPI, **kw)
        paths.append(p)
    plt.close(fig)
    return paths


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
