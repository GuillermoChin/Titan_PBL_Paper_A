#!/usr/bin/env python
"""
source_a_core.py -- computational core of the original submission (v1.0.0, make_paper_A_figures.py).

Kept unchanged in its calculations because the analyses of code/analysis/
import its functions (thresholds, mixing times, heating profiles, archive
readings). Its own figures and tables are those of the original submission and
are no longer part of the manuscript; run_all.py calls it only with
--lw-boundary/--outdir to write the base tables under build/source_a/.

Genera todas las figuras y tablas del Paper A (version enviada):

    "Radiative Timescale Constraints on Boundary Layer Evolution in Models of
     Titan"

--------------------------------------------------------------------------
DOS FUENTES DE DATOS, DELIBERADAMENTE SEPARADAS
--------------------------------------------------------------------------
FUENTE A  Nuestro framework titan_pbl. Produce las Figuras 1 a 4 y las
          Tablas 1 y 2, que corresponden a las Secciones 2 y 3 del
          manuscrito.

FUENTE B  La salida archivada de Hartwick et al. (2026). Produce las Figuras
          5 y 6 y las Tablas 3 y 4, que corresponden a las Secciones 4 y 6.

La separacion es intencional y se refleja en la estructura del script. Las
Secciones 2 y 3 del manuscrito descansan en nuestra implementacion; las
Secciones 4 y 6 descansan en datos ajenos. Un lector debe poder ver cual es
cual sin leer el codigo.

Si una de las dos fuentes no esta disponible, el script produce lo que puede
y avisa de lo que falta.

--------------------------------------------------------------------------
USO
--------------------------------------------------------------------------
    conda activate titan_pbl_paper_A
    set TITAN_DATA_DIR=<directorio que contiene "Dry Simulations/">
    python code/analysis/source_a_core.py --lw-boundary kirchhoff --outdir build/source_a/kirchhoff

Salidas en el directorio de --outdir (o PAPER_A_OUT) como PDF, EPS y CSV.

Opciones por variable de entorno:
    TITAN_DATA_DIR   ruta al archivo de datos de Hartwick et al.
    PAPER_A_OUT      directorio de salida, por defecto ./paper_A
    SKIP_ARCHIVE     si vale 1, omite la Fuente B
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Guard de importacion (H-29)
# ---------------------------------------------------------------------------
#
# El entorno conda llego a tener titan_pbl instalado en modo editable desde
# OTRA copia del codigo, fuera del repositorio y con otro Tsol. Ejecutado como
# `python code/analysis/source_a_core.py`, sys.path[0] es data/ y esa copia
# ganaba la importacion sin aviso. Si el framework esta disponible, tiene que
# ser el de este repositorio; si no lo esta, la Fuente B sigue funcionando.

_REPO_ROOT = Path(__file__).resolve().parents[2]
_SUBSET = _REPO_ROOT / "code" / "titan_pbl_subset"
if str(_SUBSET) not in sys.path:
    sys.path.insert(0, str(_SUBSET))


def _check_framework_origin() -> None:
    """Abort if titan_pbl is not the subset shipped with this package."""
    try:
        import titan_pbl
    except ImportError:
        return
    origin = Path(titan_pbl.__file__).resolve()
    if _SUBSET not in origin.parents:
        raise SystemExit(
            f"titan_pbl se importo desde {origin}, fuera de {_SUBSET}. "
            "El paquete usa solo el subconjunto incluido."
        )


_check_framework_origin()

# ---------------------------------------------------------------------------
# Configuracion global
# ---------------------------------------------------------------------------

OUT_DIR = Path(os.environ.get("PAPER_A_OUT", "./paper_A"))
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Por defecto, data/ del paquete, donde se descomprime el archivo de Zenodo.
# TITAN_DATA_DIR gana, para maquinas que lo guarden en otro volumen.
DATA_DIR = os.environ.get("TITAN_DATA_DIR", str(_REPO_ROOT / "data"))
SKIP_ARCHIVE = os.environ.get("SKIP_ARCHIVE", "0") == "1"

# Constantes del archivo de Hartwick et al., tomadas de su script de figuras.
# Nota: aqui vivia R_ARCHIVE = 290.00, nunca usada y discrepante con la
# constante del gas del preset (296.8 J kg-1 K-1, N2 puro, M = 28.014 g/mol).
# Se elimina para que no quede un tercer valor del gas flotando en el proyecto.
CP = 1101.0
GRAVITY = 1.352
SOL_INDICES = {2: (87, 231), 3: (231, 375), 4: (375, 519), 5: (519, 663)}


def _sol_length() -> float:
    """Duracion del Tsol [s], leida del preset para no duplicar la constante.

    El preset es la unica fuente de verdad. Si el framework no esta instalado
    la Fuente B debe seguir funcionando, de modo que se recurre al valor
    documentado, que es el mismo del preset y esta cubierto por un test.
    """
    try:
        from titan_pbl.config import load_config

        return float(load_config("huygens_baseline").planet.sol_length)
    except Exception:  # noqa: BLE001 - la Fuente B no debe depender del framework
        # Respaldo: 86400 pasos por Tsol x DT = 15.95 s del archivo de mtWRF.
        return 1378080.0


SOL_LENGTH = _sol_length()

# ---------------------------------------------------------------------------
# Corridas del archivo
# ---------------------------------------------------------------------------
#
# El archivo de Zenodo trae varias corridas secas y ninguna se llama de forma
# inequivoca "la baseline del paper". La eleccion importa: el experimento de
# flujo suprimido (0.1flhc) se compara contra ella en las Figuras 5 y 6 y en
# las Tablas 3 y 4, y una comparacion solo es limpia si ambas corridas
# comparten condicion de frontera profunda.
#
# El discriminante es TMN, la temperatura del subsuelo profundo, que en estos
# archivos se impone fija en la capa mas honda (TSLB[-1] = TMN durante toda la
# integracion). Comparando el estado inicial campo a campo contra 0.1flhc:
#
# Las profundidades de PBL de la tabla siguiente son LECTURAS DEL ARCHIVO de
# Hartwick et al., no resultados de este framework. El paper solo imprime
# "2353 m"; los decimales salen del NetCDF. Ver H-26.
#
#   corrida                     TMN [K]    campos identicos en t=0   PBL Tsol3
#   1D_Huygens_TMN91_5sol       91.9833          9 de 11              2353.27 m
#   1D_Huygens_baseline_5sol    92.9640          9 de 11              2353.36 m
#
# Ambas difieren de 0.1flhc (TMN = 91.865 K) solo en TSLB y TMN, pero TMN91
# difiere en 0.118 K y baseline_5sol en 1.099 K, casi un orden de magnitud mas.
# Por eso el default es TMN91. Ver H-23 en data/README_DATA.md.
#
# Se deja baseline_5sol seleccionable porque no puede descartarse que el texto
# publicado se refiera a ella; el efecto sobre las conclusiones es nulo (la
# profundidad de PBL difiere en 0.09 m, cuatro ordenes de magnitud por debajo
# del efecto que discute el Paper A) pero la eleccion queda explicita y no
# enterrada en una constante.

ARCHIVE_BASELINES = {
    # nombre -> (subdirectorio, TMN [K], nota)
    "TMN91": (
        "Dry Simulations/1D_Huygens_TMN91_5sol",
        91.9833,
        "frontera profunda mas cercana a la del experimento suprimido",
    ),
    "baseline5sol": (
        "Dry Simulations/1D_Huygens_baseline_5sol",
        92.9640,
        "corrida rotulada 'baseline' en el archivo, frontera 1.1 K mas calida",
    ),
}

BASELINE_RUN = os.environ.get("TITAN_BASELINE_RUN", "TMN91")
if BASELINE_RUN not in ARCHIVE_BASELINES:
    raise SystemExit(
        f"TITAN_BASELINE_RUN='{BASELINE_RUN}' no es valido. "
        f"Opciones: {sorted(ARCHIVE_BASELINES)}"
    )

# Los nombres publicados de los archivos llevan "09:30:00". Los dos puntos son
# ilegales en nombres de archivo de Windows, de modo que en disco aparecen como
# "09_30_00". Se conserva el nombre publicado y resolve_archive() se encarga.
ARCHIVE_RUNS = {
    "baseline": (
        f"{ARCHIVE_BASELINES[BASELINE_RUN][0]}/wrfout_d01_2007-00001_09:30:00"
    ),
    "suppressed": (
        "Dry Simulations/1D_Huygens_0.1flhc_5sol/wrfout_d01_2007-00001_09:30:00"
    ),
}

# Umbrales calculados en la Seccion 3, recalculados abajo pero fijados aqui
# como respaldo si la Fuente A no esta disponible.
FALLBACK_THRESHOLDS = {"heating": 2.15, "destabilisation": 2.54, "analytic": 1.20}

# Estratificacion HASI entre superficie y 4 km, Seccion 3.5.
DELTA_THETA_4KM = 0.652

# Frontera LW inferior (H-27). Vacio = la que declara el preset, que es el
# comportamiento historico del script. main() lo llena con --lw-boundary.
LW_BOUNDARY_OVERRIDES: dict[str, str] = {}

# Variantes de onda corta para separar LW y SW (revision R1, T5):
#   "lw"      solo onda larga;
#   "noon"    onda larga + onda corta al mediodia del Tsol 1 (lo publicado);
#   "diurnal" onda larga + promedio diurno de la onda corta en el Tsol 1.
SW_VARIANTS = ("lw", "noon", "diurnal")
DIURNAL_SAMPLES = 48


def _framework_overrides() -> dict[str, str] | None:
    """Overrides for load_config; None keeps the historical call unchanged."""
    return dict(LW_BOUNDARY_OVERRIDES) or None


# ---------------------------------------------------------------------------
# Estilo de figura para publicacion
# ---------------------------------------------------------------------------


def configure_style() -> None:
    """Estilo sobrio, legible en escala de grises y apto para impresion."""
    import matplotlib as mpl

    mpl.use("Agg")
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
            "mathtext.fontset": "dejavuserif",
            "font.size": 9,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "legend.frameon": False,
            "axes.linewidth": 0.7,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "lines.linewidth": 1.3,
            "figure.dpi": 200,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
        }
    )


# Paleta segura para daltonismo, de Okabe e Ito.
PALETTE = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "cyan": "#56B4E9",
    "yellow": "#F0E442",
    "grey": "#5A5A5A",
}


def save(fig, name: str) -> None:
    """Guarda en PDF y EPS vectoriales, y en PNG a 600 dpi para Word.

    P&SS acepta vectorial para produccion; el PNG existe porque las figuras se
    insertan a mano en el .docx y Word maneja mal el EPS.
    """
    for extension in ("pdf", "eps"):
        fig.savefig(OUT_DIR / f"{name}.{extension}", format=extension)
    fig.savefig(OUT_DIR / f"{name}.png", format="png", dpi=600)
    print(f"    guardada  {name}.pdf  {name}.eps  {name}.png")


# ===========================================================================
# FUENTE A -- Nuestro framework
# ===========================================================================


def load_framework():
    """Importa el framework y devuelve un contexto reutilizable.

    Devuelve None si el paquete no esta instalado, de modo que el script
    pueda seguir con la Fuente B.
    """
    try:
        from titan_pbl.config import load_config
        from titan_pbl.core import build_grid, initialize_state
        from titan_pbl.schemes.radiation import GrayRadiation
    except ImportError as exc:
        print(f"  El framework titan_pbl no esta disponible: {exc}")
        return None

    config = load_config("huygens_baseline", overrides=_framework_overrides())
    state = initialize_state(config, build_grid(config))
    scheme = GrayRadiation(config)
    scheme.initialize(config, state)

    # Mediodia local.
    state.time = (12.0 - config.site.start_local_time) / 24.0 * config.planet.sol_length
    _, midpoints = state.geopotential_heights()

    return {
        "config": config,
        "state": state,
        "scheme": scheme,
        "altitude": midpoints,
        "air_temperature": float(state.atmospheric_temperature[0]),
        "load_config": load_config,
        "build_grid": build_grid,
        "initialize_state": initialize_state,
        "GrayRadiation": GrayRadiation,
    }


def heating_profile(ctx, contrast: float) -> np.ndarray:
    """Perfil de calentamiento radiativo [K/Tsol] para un contraste dado."""
    ctx["state"].skin_temperature = ctx["air_temperature"] + contrast
    return ctx["scheme"].heating_rate(ctx["state"]) * SOL_LENGTH


def destabilisation_rate(ctx, contrast: float, top_m: float = 4000.0) -> float:
    """Diferencia de calentamiento entre la superficie y una altura dada."""
    rates = heating_profile(ctx, contrast)
    top = int(np.argmin(np.abs(ctx["altitude"] - top_m)))
    return float(rates[0] - rates[top])


def compute_thresholds(ctx) -> dict[str, float]:
    """Localiza los umbrales de calentamiento y de desestabilizacion."""
    from scipy.optimize import brentq

    air = ctx["air_temperature"]
    emissivity = ctx["config"].radiation.surface_emissivity
    analytic = air * (emissivity ** (-0.25) - 1.0)

    thresholds = {"analytic": analytic}
    try:
        thresholds["heating"] = brentq(
            lambda d: heating_profile(ctx, d)[0], 0.05, 20.0, xtol=1e-6
        )
    except ValueError:
        thresholds["heating"] = np.nan
    try:
        thresholds["destabilisation"] = brentq(
            lambda d: destabilisation_rate(ctx, d), 0.05, 20.0, xtol=1e-6
        )
    except ValueError:
        thresholds["destabilisation"] = np.nan
    return thresholds


# --- Separacion de onda larga y onda corta (revision R1, T5) --------------
#
# Las funciones de arriba usan heating_rate, que suma onda larga y onda corta
# al mediodia. Las de abajo descomponen esa suma sin alterar las anteriores:
# la variante "noon" delega en heating_profile, de modo que reproduce lo
# publicado bit a bit.


def heating_components(ctx, contrast: float) -> dict[str, np.ndarray]:
    """Longwave and shortwave heating [K/Tsol] at the state's time of day.

    The contrast is imposed through the skin temperature only, exactly as in
    :func:`heating_profile`.
    """
    state, scheme, config = ctx["state"], ctx["scheme"], ctx["config"]
    state.skin_temperature = ctx["air_temperature"] + contrast
    gravity = config.planet.gravity
    cp = config.atmosphere.gas.specific_heat_cp
    dp = state.grid.atmosphere.delta_pressure
    lw_down, lw_up = scheme.longwave_fluxes(state)
    sw_down, sw_up = scheme.shortwave_fluxes(state)
    lw_net = lw_down - lw_up
    sw_net = sw_down - sw_up
    return {
        "lw": gravity * (lw_net[1:] - lw_net[:-1]) / (cp * dp) * SOL_LENGTH,
        "sw": gravity * (sw_net[1:] - sw_net[:-1]) / (cp * dp) * SOL_LENGTH,
    }


def diurnal_mean_shortwave_heating(ctx, samples: int = DIURNAL_SAMPLES) -> np.ndarray:
    """Shortwave heating [K/Tsol] averaged over Tsol 1 on the fixed state.

    ``samples`` equally spaced evaluations from the start of the sol. The
    shortwave depends on the solar geometry and the visible optical depth
    only, not on temperatures, so the result is cached in ``ctx``.
    """
    key = f"_sw_diurnal_{samples}"
    if key not in ctx:
        state, scheme, config = ctx["state"], ctx["scheme"], ctx["config"]
        saved = state.time
        gravity = config.planet.gravity
        cp = config.atmosphere.gas.specific_heat_cp
        dp = state.grid.atmosphere.delta_pressure
        total = np.zeros(state.grid.n_atmospheric_layers)
        for i in range(samples):
            state.time = i * config.planet.sol_length / samples
            sw_down, sw_up = scheme.shortwave_fluxes(state)
            sw_net = sw_down - sw_up
            total += gravity * (sw_net[1:] - sw_net[:-1]) / (cp * dp) * SOL_LENGTH
        state.time = saved
        ctx[key] = total / samples
    return ctx[key]


def heating_profile_variant(ctx, contrast: float, variant: str) -> np.ndarray:
    """Heating profile [K/Tsol] for one of the shortwave variants in SW_VARIANTS."""
    if variant == "noon":
        return heating_profile(ctx, contrast)
    components = heating_components(ctx, contrast)
    if variant == "lw":
        return components["lw"]
    if variant == "diurnal":
        return components["lw"] + diurnal_mean_shortwave_heating(ctx)
    raise ValueError(f"variant must be one of {SW_VARIANTS}, got {variant!r}")


def destabilisation_rate_variant(
    ctx, contrast: float, variant: str, top_m: float = 4000.0
) -> float:
    """Same definition as :func:`destabilisation_rate`, for a shortwave variant."""
    rates = heating_profile_variant(ctx, contrast, variant)
    top = int(np.argmin(np.abs(ctx["altitude"] - top_m)))
    return float(rates[0] - rates[top])


def compute_thresholds_variant(ctx, variant: str) -> dict[str, float]:
    """Same brackets and tolerance as :func:`compute_thresholds`, per variant."""
    from scipy.optimize import brentq

    air = ctx["air_temperature"]
    emissivity = ctx["config"].radiation.surface_emissivity
    thresholds = {"analytic": air * (emissivity ** (-0.25) - 1.0)}
    try:
        thresholds["heating"] = brentq(
            lambda d: heating_profile_variant(ctx, d, variant)[0], 0.05, 20.0, xtol=1e-6
        )
    except ValueError:
        thresholds["heating"] = np.nan
    try:
        thresholds["destabilisation"] = brentq(
            lambda d: destabilisation_rate_variant(ctx, d, variant), 0.05, 20.0, xtol=1e-6
        )
    except ValueError:
        thresholds["destabilisation"] = np.nan
    return thresholds


def mixing_times_variant(ctx, thresholds: dict[str, float], variant: str) -> list[dict]:
    """Same method as :func:`table_thresholds`: Delta_theta / rate at fixed contrasts.

    The rate is evaluated once and treated as constant over the whole sol.
    """
    rows = []
    for contrast in (thresholds["destabilisation"], 3.0, 3.5, 4.5):
        if not np.isfinite(contrast):
            rows.append({"contrast_K": contrast, "destabilisation_K_per_sol": np.nan,
                         "mixing_time_sols": np.nan})
            continue
        rate = destabilisation_rate_variant(ctx, float(contrast), variant)
        rows.append(
            {
                "contrast_K": contrast,
                "destabilisation_K_per_sol": rate,
                "mixing_time_sols": DELTA_THETA_4KM / rate if rate > 0 else np.inf,
            }
        )
    return rows


# --- Figura 1 -------------------------------------------------------------


def figure1_schematic() -> None:
    """Esquema del balance radiativo de una capa gris delgada.

    Es la unica figura del manuscrito que no representa datos. Su proposito
    es hacer visualmente evidente que la emision se duplica mientras cada
    termino de absorcion es simple, que es el nucleo del argumento de la
    Seccion 2.
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrow, Rectangle

    print("  Figura 1 -- esquema del balance radiativo")
    fig, ax = plt.subplots(figsize=(3.5, 4.1))

    # El trazado separa los cuatro terminos horizontalmente para que ninguna
    # etiqueta se superponga con otra.
    surface_top, layer_bottom, layer_top, sky = 0.28, 0.56, 0.70, 0.88

    ax.add_patch(Rectangle((0.06, 0.18), 0.88, 0.10, facecolor="#C8A882",
                           edgecolor="k", linewidth=0.7))
    ax.text(0.5, 0.23, r"Surface,  $T_s$,  emissivity $\varepsilon$",
            ha="center", va="center", fontsize=8)

    ax.add_patch(Rectangle((0.06, layer_bottom), 0.88, layer_top - layer_bottom,
                           facecolor="#DCE9F2", edgecolor="k", linewidth=0.7))
    ax.text(0.5, 0.63, r"Layer,  $T$,  thickness $d\tau$", ha="center",
            va="center", fontsize=8)

    arrow = dict(width=0.007, head_width=0.032, head_length=0.030,
                 length_includes_head=True)

    # Absorcion desde abajo, columna izquierda.
    ax.add_patch(FancyArrow(0.17, surface_top, 0.0, layer_bottom - surface_top,
                            color=PALETTE["red"], **arrow))
    ax.text(0.20, 0.42, r"$F^{\uparrow}=\varepsilon\sigma T_s^{4}$", fontsize=8,
            color=PALETTE["red"], ha="left", va="center")

    # Absorcion desde arriba, columna derecha.
    ax.add_patch(FancyArrow(0.86, sky, 0.0, layer_top - sky,
                            color=PALETTE["blue"], **arrow))
    ax.text(0.83, 0.80, r"$F^{\downarrow}\approx\sigma T_a^{4}$", fontsize=8,
            color=PALETTE["blue"], ha="right", va="center")

    # Emision hacia arriba, columna centro-izquierda.
    ax.add_patch(FancyArrow(0.40, layer_top, 0.0, sky - layer_top,
                            color=PALETTE["grey"], **arrow))
    ax.text(0.43, 0.80, r"$\sigma T^{4}$", fontsize=8, color=PALETTE["grey"],
            ha="left", va="center")

    # Emision hacia abajo, columna centro-derecha.
    ax.add_patch(FancyArrow(0.60, layer_bottom, 0.0, surface_top - layer_bottom,
                            color=PALETTE["grey"], **arrow))
    ax.text(0.63, 0.42, r"$\sigma T^{4}$", fontsize=8, color=PALETTE["grey"],
            ha="left", va="center")

    ax.text(0.5, 0.965,
            "Absorption from one direction each,\nemission into both",
            ha="center", va="top", fontsize=8, style="italic",
            color=PALETTE["grey"])

    ax.text(0.5, 0.07,
            r"$\dfrac{dE}{dt}=d\tau \left(F^{\uparrow}+F^{\downarrow}-2\sigma T^{4}\right)$",
            ha="center", va="center", fontsize=10)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    save(fig, "figure1_schematic")
    plt.close(fig)


# --- Figura 2 -------------------------------------------------------------


def figure2_heating_profiles(ctx) -> None:
    """Perfiles de calentamiento radiativo para varios contrastes termicos."""
    import matplotlib.pyplot as plt

    print("  Figura 2 -- perfiles de calentamiento")
    fig, ax = plt.subplots(figsize=(3.4, 4.0))

    contrasts = [0.0, 1.0, 2.0, 3.0, 4.0]
    colours = [PALETTE["blue"], PALETTE["cyan"], PALETTE["green"],
               PALETTE["orange"], PALETTE["red"]]
    # Estilo de linea ademas de color: las cinco series deben distinguirse en
    # escala de grises, que es como se imprime buena parte de la tirada.
    styles = ["-", "--", "-.", ":", (0, (5, 1, 1, 1))]
    altitude = ctx["altitude"]
    mask = altitude <= 4000.0

    for contrast, colour, style in zip(contrasts, colours, styles):
        rates = heating_profile(ctx, contrast)
        ax.plot(rates[mask], altitude[mask], color=colour, linestyle=style,
                label=rf"$\Delta T = {contrast:.0f}$ K")

    ax.axvline(0.0, color="k", linewidth=0.6, linestyle=":")
    ax.set_xlabel(r"Radiative heating rate [K per Titan sol]")
    ax.set_ylabel("Altitude [m]")
    ax.set_ylim(0, 4000)
    # Margen horizontal para que ninguna curva quede pegada al marco. No
    # recorta datos: solo ensancha el rango.
    ax.margins(x=0.06)
    ax.legend(loc="upper left")
    save(fig, "figure2_heating_profiles")
    plt.close(fig)


# --- Figura 3 -------------------------------------------------------------


def figure3_destabilisation(ctx, thresholds: dict[str, float]) -> None:
    """Tasa de desestabilizacion frente al contraste termico."""
    import matplotlib.pyplot as plt

    print("  Figura 3 -- tasa de desestabilizacion")
    fig, ax = plt.subplots(figsize=(3.6, 3.0))

    contrasts = np.linspace(0.0, 5.0, 60)
    rates = np.array([destabilisation_rate(ctx, float(c)) for c in contrasts])

    ax.plot(contrasts, rates, color=PALETTE["blue"])
    ax.axhline(0.0, color="k", linewidth=0.6, linestyle=":")

    zero = thresholds.get("destabilisation", np.nan)
    if np.isfinite(zero):
        ax.axvline(zero, color=PALETTE["red"], linewidth=0.9, linestyle="--")
        ax.annotate(rf"Destabilisation" "\n" rf"threshold, {zero:.2f} K",
                    xy=(zero, 0.0), xytext=(zero + 0.25, rates.min() * 0.55),
                    fontsize=7.5, color=PALETTE["red"],
                    arrowprops=dict(arrowstyle="-", color=PALETTE["red"],
                                    linewidth=0.6))

    # Banda de contrastes recuperados del archivo, Seccion 6.
    ax.axvspan(0.22, 0.69, color=PALETTE["grey"], alpha=0.16, linewidth=0)
    ax.annotate("Contrasts recovered\nfrom archived output",
                xy=(0.46, rates.max() * 0.55),
                xytext=(1.45, rates.max() * 0.86),
                fontsize=7.5, color=PALETTE["grey"], ha="left", va="top",
                arrowprops=dict(arrowstyle="-", color=PALETTE["grey"],
                                linewidth=0.6))

    ax.set_xlabel(r"Surface-to-air temperature contrast $\Delta T$ [K]")
    ax.set_ylabel("Destabilisation rate\n[K per Titan sol]")
    ax.set_xlim(0, 5)
    ax.margins(y=0.06)
    save(fig, "figure3_destabilisation")
    plt.close(fig)


# --- Figura 4 -------------------------------------------------------------


def figure4_mixing_time(ctx, thresholds: dict[str, float]) -> None:
    """Tiempo de mezcla frente al contraste termico, ejes logaritmicos."""
    import matplotlib.pyplot as plt

    print("  Figura 4 -- tiempo de mezcla")
    fig, ax = plt.subplots(figsize=(3.6, 3.0))

    zero = thresholds.get("destabilisation", FALLBACK_THRESHOLDS["destabilisation"])
    contrasts = np.linspace(zero + 0.05, 8.0, 90)

    times = []
    for contrast in contrasts:
        rate = destabilisation_rate(ctx, float(contrast))
        times.append(DELTA_THETA_4KM / rate if rate > 0 else np.nan)
    times = np.array(times)

    ax.plot(contrasts, times, color=PALETTE["blue"])
    ax.axhspan(1.0, 5.0, color=PALETTE["green"], alpha=0.18, linewidth=0)
    # La anotacion estaba en x = 7.7, donde la curva la atravesaba. Se mueve al
    # extremo izquierdo de la banda, que esta vacio: la curva solo existe a
    # partir del umbral de desestabilizacion, en torno a 2.6 K.
    ax.text(0.95, 1.35, "Typical 1D\nintegration length", ha="left", va="bottom",
            fontsize=7.5, color=PALETTE["green"])

    ax.axvspan(0.22, 0.69, color=PALETTE["grey"], alpha=0.16, linewidth=0)
    ax.set_yscale("log")
    ax.set_xlabel(r"Surface-to-air temperature contrast $\Delta T$ [K]")
    ax.set_ylabel("Time to mix the lowest 4 km\n[Titan sols]")
    ax.set_xlim(0, 8)
    ax.margins(x=0.02)
    save(fig, "figure4_mixing_time")
    plt.close(fig)


# --- Tabla 1 --------------------------------------------------------------


def table1_configuration(ctx) -> None:
    """Configuracion del modelo con su procedencia."""
    import pandas as pd

    print("  Tabla 1 -- configuracion del modelo")
    config = ctx["config"]
    radiation = config.radiation
    grid = config.atmosphere.vertical_grid

    rows = [
        ("Solar constant", "S0", radiation.solar_constant, "W m-2", "Hartwick+ 2026"),
        ("Downward IR at model top", "I_IR", radiation.ir_top_flux, "W m-2", "Hartwick+ 2026"),
        ("Surface albedo", "alpha", radiation.surface_albedo, "-", "Hartwick+ 2026"),
        ("Surface emissivity", "epsilon", radiation.surface_emissivity, "-", "Hartwick+ 2026"),
        ("Visible optical depth", "tau_vis", radiation.optical_depth_visible, "-", "Chatain+ 2022"),
        ("Infrared optical depth", "tau_IR", radiation.optical_depth_ir, "-", "Hartwick+ 2026"),
        ("Visible scaling exponent", "gamma", radiation.visible_scaling_exponent, "-", "Schneider+ 2012"),
        ("IR weighting coefficient", "f", radiation.ir_weighting_coefficient, "-", "Chatain+ 2022"),
        ("Single scattering albedo", "omega0", radiation.single_scattering_albedo, "-", "Chatain+ 2022"),
        ("Asymmetry factor", "g", radiation.asymmetry_factor, "-", "Chatain+ 2022"),
        ("Atmospheric layers", "N", grid.n_layers, "-", "Hartwick+ 2026"),
        ("Surface pressure", "p_s", grid.surface_pressure, "Pa", "Hartwick+ 2026"),
        ("Model top pressure", "p_t", grid.top_pressure, "Pa", "Hartwick+ 2026"),
        ("Specific heat", "c_p", config.atmosphere.gas.specific_heat_cp, "J kg-1 K-1", "Hartwick+ 2026"),
        ("Surface gravity", "g", config.planet.gravity, "m s-2", "standard"),
    ]
    table = pd.DataFrame(rows, columns=["parameter", "symbol", "value", "units", "source"])
    table.to_csv(OUT_DIR / "table1_configuration.csv", index=False)
    print(table.to_string(index=False))


# --- Tabla 2 --------------------------------------------------------------


def table2_parameter_sweep(ctx) -> None:
    """Barrido sobre los tres parametros del presupuesto de onda larga."""
    import pandas as pd

    print("\n  Tabla 2 -- barrido de parametros infrarrojos")
    load_config = ctx["load_config"]
    build_grid = ctx["build_grid"]
    initialize_state = ctx["initialize_state"]
    GrayRadiation = ctx["GrayRadiation"]

    rows = []
    for ir_top in (0.0, 2.0, 4.0, 8.0):
        for tau_ir in (10.0, 30.0):
            for weight in (0.15, 0.50):
                config = load_config(
                    "huygens_baseline",
                    overrides={
                        "radiation.ir_top_flux": ir_top,
                        "radiation.optical_depth_ir": tau_ir,
                        "radiation.ir_weighting_coefficient": weight,
                        **LW_BOUNDARY_OVERRIDES,
                    },
                )
                state = initialize_state(config, build_grid(config))
                scheme = GrayRadiation(config)
                scheme.initialize(config, state)
                state.time = (
                    (12.0 - config.site.start_local_time) / 24.0 * config.planet.sol_length
                )
                _, altitude = state.geopotential_heights()
                rates = scheme.heating_rate(state) * SOL_LENGTH
                below = altitude < 2500.0
                rows.append(
                    {
                        "I_IR_W_m2": ir_top,
                        "tau_IR": tau_ir,
                        "f": weight,
                        "max_heating_K_per_sol": float(np.max(rates[below])),
                        "min_heating_K_per_sol": float(np.min(rates[below])),
                        "surface_heating_K_per_sol": float(rates[0]),
                        "gradient_K_per_sol": float(rates[0] - rates[below][-1]),
                    }
                )

    table = pd.DataFrame(rows)
    table.to_csv(OUT_DIR / "table2_parameter_sweep.csv", index=False)
    print(table.to_string(index=False, float_format=lambda v: f"{v:10.5f}"))
    print(f"\n    Maximo calentamiento en todo el barrido: "
          f"{table['max_heating_K_per_sol'].max():+.5f} K/sol")
    print(f"    Todos los gradientes negativos: "
          f"{bool((table['gradient_K_per_sol'] < 0).all())}")


# --- Tabla de umbrales ----------------------------------------------------


def table_thresholds(ctx, thresholds: dict[str, float]) -> None:
    """Los tres umbrales de la Seccion 3.3."""
    import pandas as pd

    print("\n  Umbrales de la Seccion 3.3")
    rows = [
        ("Analytic heating threshold, Eq. (4)", thresholds["analytic"]),
        ("Numerical heating threshold, lowest layer", thresholds["heating"]),
        ("Numerical destabilisation threshold, 0-4 km", thresholds["destabilisation"]),
    ]
    table = pd.DataFrame(rows, columns=["quantity", "contrast_K"])
    table.to_csv(OUT_DIR / "table_thresholds.csv", index=False)
    print(table.to_string(index=False, float_format=lambda v: f"{v:8.3f}"))

    print("\n  Tiempos de mezcla de la Seccion 3.5")
    rows = []
    for contrast in (thresholds["destabilisation"], 3.0, 3.5, 4.5):
        rate = destabilisation_rate(ctx, float(contrast))
        rows.append(
            {
                "contrast_K": contrast,
                "destabilisation_K_per_sol": rate,
                "mixing_time_sols": DELTA_THETA_4KM / rate if rate > 0 else np.inf,
            }
        )
    table = pd.DataFrame(rows)
    table.to_csv(OUT_DIR / "table_mixing_times.csv", index=False)
    print(table.to_string(index=False, float_format=lambda v: f"{v:10.4f}"))


# ===========================================================================
# FUENTE B -- El archivo de Hartwick et al. (2026)
# ===========================================================================


def resolve_archive(relpath: str) -> Path:
    """Localiza un wrfout del archivo tolerando el renombrado de Windows.

    Los nombres publicados llevan dos puntos ("09:30:00"), ilegales en NTFS. Al
    descomprimir en Windows quedan como "09_30_00". Se intenta el nombre tal
    cual, luego la variante con guiones bajos, luego el directorio.
    """
    variants = [relpath.replace("\\", "/")]
    if ":" in relpath:
        variants.append(relpath.replace("\\", "/").replace(":", "_"))

    for variant in variants:
        candidate = Path(DATA_DIR) / variant
        if candidate.is_file():
            return candidate
        if candidate.is_dir():
            matches = sorted(candidate.glob("wrfout_d01*")) or sorted(
                candidate.rglob("wrfout_d01*")
            )
            if matches:
                return matches[0]
        parent = candidate.parent
        if parent.is_dir():
            matches = sorted(parent.glob(candidate.name + "*"))
            if matches:
                return matches[0]

    # Ultimo recurso: el directorio que contiene la corrida, ignorando el
    # nombre del archivo, que es lo unico que cambia entre plataformas.
    directory = Path(DATA_DIR) / Path(relpath.replace("\\", "/")).parent
    if directory.is_dir():
        matches = sorted(directory.glob("wrfout_d01*"))
        if matches:
            return matches[0]

    raise FileNotFoundError(
        f"No se encontro '{relpath}' bajo {DATA_DIR}. "
        "Revisa TITAN_DATA_DIR y que el archivo de Zenodo este descomprimido."
    )


def load_archive(label: str):
    import xarray as xr

    return xr.open_dataset(resolve_archive(ARCHIVE_RUNS[label]), decode_times=None)


def archive_noon_index(ds, sol: int) -> int:
    start, end = SOL_INDICES[sol]
    times = [t.decode("utf-8") if isinstance(t, bytes) else str(t) for t in ds.Times.values]
    matches = [i for i in range(start, end) if "_12:00" in times[i]]
    return matches[0] if matches else (start + end) // 2


def archive_heating(ds, index: int):
    """Perfiles de calentamiento LW y SW del archivo, con el signo correcto.

    El signo negativo se sigue de que la masa de una capa es dp/g y la
    convergencia de energia es la diferencia entre el flujo descendente que
    entra por su tope y el que sale por su base, con la presion decreciendo
    hacia arriba.
    """
    lw_net = np.asarray(
        ds.FLUX_LW_D[index, :, 0, 0] - ds.FLUX_LW_U[index, :, 0, 0], dtype=np.float64
    )
    sw_net = np.asarray(
        ds.FLUX_SW_D[index, :, 0, 0] - ds.FLUX_SW_U[index, :, 0, 0], dtype=np.float64
    )
    pressure = np.asarray(ds.P[index, :, 0, 0] + ds.PB[index, :, 0, 0], dtype=np.float64)
    altitude = np.asarray(
        (ds.PH[index, 1:, 0, 0] + ds.PHB[index, 1:, 0, 0]) / GRAVITY, dtype=np.float64
    )

    n = min(len(lw_net), len(pressure), len(altitude))
    dp = np.gradient(pressure[:n])
    lw_rate = -(GRAVITY / CP) * np.gradient(lw_net[:n]) / dp * SOL_LENGTH
    sw_rate = -(GRAVITY / CP) * np.gradient(sw_net[:n]) / dp * SOL_LENGTH
    return altitude[:n], lw_rate, sw_rate


# --- Figura 5 -------------------------------------------------------------


def figure5_validation(ctx) -> None:
    """Perfiles archivados frente a los calculados, para ambas corridas."""
    import matplotlib.pyplot as plt

    print("  Figura 5 -- validacion contra el archivo")
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.6), sharey=True)

    panel_labels = ("(a)", "(b)")
    for ax, label, title, panel in zip(
        axes, ("baseline", "suppressed"), ("Baseline", "Suppressed flux"),
        panel_labels,
    ):
        ds = load_archive(label)
        altitude, lw_rate, sw_rate = archive_heating(ds, archive_noon_index(ds, 3))
        mask = altitude <= 2500.0

        ax.plot(lw_rate[mask], altitude[mask], color=PALETTE["blue"],
                label="Archived, longwave")
        ax.plot((lw_rate + sw_rate)[mask], altitude[mask], color=PALETTE["cyan"],
                linestyle="-.", label="Archived, total")

        if ctx is not None:
            ours = heating_profile(ctx, 0.0)
            our_altitude = ctx["altitude"]
            our_mask = our_altitude <= 2500.0
            ax.plot(ours[our_mask], our_altitude[our_mask], color=PALETTE["red"],
                    linestyle="--", label="This work")

        ax.axvline(0.0, color="k", linewidth=0.6, linestyle=":")
        ax.set_xlabel("Radiative heating rate\n[K per Titan sol]")
        ax.set_title(title, fontsize=9)
        # Etiqueta de panel en coordenadas de ejes, fuera del area de datos.
        ax.text(0.02, 0.98, panel, transform=ax.transAxes, fontsize=9,
                fontweight="bold", ha="left", va="top")
        ax.margins(x=0.05)
        ds.close()

    axes[0].set_ylabel("Altitude [m]")
    axes[0].set_ylim(0, 2500)
    axes[0].legend(loc="lower left")
    save(fig, "figure5_validation")
    plt.close(fig)


# --- Figura 6 -------------------------------------------------------------


def figure6_timeseries() -> None:
    """Contraste termico y profundidad de PBL, ambas corridas, Tsols 2 a 5."""
    import matplotlib.pyplot as plt

    print("  Figura 6 -- series temporales del archivo")
    fig, axes = plt.subplots(2, 1, figsize=(7.0, 4.6), sharex=True)

    start = SOL_INDICES[2][0]
    end = SOL_INDICES[5][1]

    for label, colour, name, style in (
        ("baseline", PALETTE["blue"], "Baseline", "-"),
        ("suppressed", PALETTE["red"], "Suppressed flux", "--"),
    ):
        ds = load_archive(label)
        skin = np.asarray(ds.TSK[start:end, 0, 0], dtype=np.float64)
        air = np.asarray(ds.T_PHY[start:end, 0, 0, 0], dtype=np.float64)
        depth = np.asarray(ds.PBLH[start:end, 0, 0], dtype=np.float64)
        sols = 2.0 + np.arange(len(skin)) / 144.0

        axes[0].plot(sols, skin - air, color=colour, linestyle=style, label=name)
        axes[1].plot(sols, depth, color=colour, linestyle=style, label=name)
        ds.close()

    thresholds = FALLBACK_THRESHOLDS
    axes[0].axhline(thresholds["heating"], color=PALETTE["grey"], linestyle="--",
                    linewidth=0.8)
    # Estaba en x = 5.55 con el texto creciendo a la derecha, de modo que se
    # salia del lienzo y la linea discontinua lo tachaba. Se centra bajo la
    # linea, en una zona vacia: las curvas no pasan de 0.7 K.
    axes[0].text(4.0, thresholds["heating"] - 0.08, "Heating threshold",
                 fontsize=7.5, color=PALETTE["grey"], ha="center", va="top")
    axes[0].set_ylabel(r"$T_s - T_a$ [K]")
    axes[0].set_ylim(-0.2, 2.6)
    # A la izquierda chocaba con la etiqueta de panel (a). Arriba a la derecha
    # queda por encima de la linea de umbral y muy lejos de las curvas, que no
    # pasan de 0.7 K.
    axes[0].legend(loc="upper right")

    axes[1].set_ylabel("Boundary layer depth [m]")
    axes[1].set_xlabel("Time [Titan sols]")
    axes[1].set_xlim(2, 6)

    for ax, panel in zip(axes, ("(a)", "(b)")):
        ax.text(0.008, 0.96, panel, transform=ax.transAxes, fontsize=9,
                fontweight="bold", ha="left", va="top")

    save(fig, "figure6_timeseries")
    plt.close(fig)


# --- Tablas 3 y 4 ---------------------------------------------------------


def table3_validation(ctx) -> None:
    """Comparacion de tasas de calentamiento archivadas y calculadas."""
    import pandas as pd

    print("\n  Tabla 3 -- validacion")
    ds = load_archive("baseline")
    altitude, lw_rate, _ = archive_heating(ds, archive_noon_index(ds, 3))
    ds.close()

    ours = heating_profile(ctx, 0.0) if ctx is not None else None
    our_altitude = ctx["altitude"] if ctx is not None else None

    rows = []
    for target in (5.4, 170.6, 541.2, 1382.1):
        k = int(np.argmin(np.abs(altitude - target)))
        archived = float(lw_rate[k])
        entry = {
            "altitude_m": float(altitude[k]),
            "archived_K_per_sol": archived,
        }
        if ours is not None:
            j = int(np.argmin(np.abs(our_altitude - target)))
            mine = float(ours[j])
            entry["this_work_K_per_sol"] = mine
            entry["difference_percent"] = (
                abs(archived - mine) / abs(archived) * 100.0 if archived else np.nan
            )
        rows.append(entry)

    table = pd.DataFrame(rows)
    table.to_csv(OUT_DIR / "table3_validation.csv", index=False)
    print(table.to_string(index=False, float_format=lambda v: f"{v:10.5f}"))


def table4_contrasts() -> None:
    """Contrastes termicos recuperados del archivo."""
    import pandas as pd

    print("\n  Tabla 4 -- contrastes termicos")
    base = load_archive("baseline")
    supp = load_archive("suppressed")

    skin_base = np.asarray(base.TSK[:, 0, 0], dtype=np.float64)
    air_base = np.asarray(base.T_PHY[:, 0, 0, 0], dtype=np.float64)
    skin_supp = np.asarray(supp.TSK[:, 0, 0], dtype=np.float64)
    air_supp = np.asarray(supp.T_PHY[:, 0, 0, 0], dtype=np.float64)
    flux_base = np.asarray(base.HFX[:, 0, 0], dtype=np.float64)
    flux_supp = np.asarray(supp.HFX[:, 0, 0], dtype=np.float64)
    depth_base = np.asarray(base.PBLH[:, 0, 0], dtype=np.float64)
    depth_supp = np.asarray(supp.PBLH[:, 0, 0], dtype=np.float64)

    rows = []
    for sol, (i, j) in SOL_INDICES.items():
        rows.append(
            {
                "sol": sol,
                "contrast_baseline_max_K": float(np.max(skin_base[i:j] - air_base[i:j])),
                "contrast_suppressed_max_K": float(np.max(skin_supp[i:j] - air_supp[i:j])),
                "delta_Ts_between_runs_max_K": float(np.max(skin_supp[i:j] - skin_base[i:j])),
                "delta_Ts_between_runs_mean_K": float(np.mean(skin_supp[i:j] - skin_base[i:j])),
                "HFX_baseline_peak_W_m2": float(np.max(flux_base[i:j])),
                "HFX_suppressed_peak_W_m2": float(np.max(flux_supp[i:j])),
                "HFX_ratio": float(np.max(flux_supp[i:j]) / np.max(flux_base[i:j])),
                "PBL_baseline_max_m": float(np.max(depth_base[i:j])),
                "PBL_suppressed_max_m": float(np.max(depth_supp[i:j])),
            }
        )

    table = pd.DataFrame(rows)
    table.to_csv(OUT_DIR / "table4_contrasts.csv", index=False)
    print(table.to_string(index=False, float_format=lambda v: f"{v:10.4f}"))

    drift = (
        table["PBL_suppressed_max_m"].iloc[-1] / table["PBL_suppressed_max_m"].iloc[0] - 1
    ) * 100.0
    print(f"\n    Deriva de la PBL suprimida, Tsols 2 a 5: {drift:+.1f} %")

    base.close()
    supp.close()


def archive_contrast_maxima() -> dict[str, dict]:
    """Maximum surface-to-air contrast TSK - T_PHY[0] per Tsol and over Tsols 2-5.

    For both archive runs (revision R1, decision D4). Reads the archive only;
    independent of the radiation scheme and of the sol length.
    """
    out: dict[str, dict] = {}
    for label in ("baseline", "suppressed"):
        ds = load_archive(label)
        contrast = np.asarray(ds.TSK[:, 0, 0], dtype=np.float64) - np.asarray(
            ds.T_PHY[:, 0, 0, 0], dtype=np.float64
        )
        ds.close()
        per_sol = {sol: float(np.max(contrast[i:j])) for sol, (i, j) in SOL_INDICES.items()}
        best = max(per_sol, key=per_sol.get)
        out[label] = {
            "file": str(resolve_archive(ARCHIVE_RUNS[label])),
            "per_sol_max_K": per_sol,
            "max_sols_2_to_5_K": per_sol[best],
            "sol_of_max": best,
        }
    return out


# ===========================================================================
# Ejecucion
# ===========================================================================


def parse_args(argv: list[str] | None = None):
    """Command line options. The defaults reproduce the historical behaviour."""
    import argparse

    parser = argparse.ArgumentParser(description="Figuras y tablas del Paper A.")
    parser.add_argument(
        "--lw-boundary",
        choices=("full_absorption", "kirchhoff"),
        default=None,
        help="Frontera LW inferior (H-27). Por defecto, la que declara el preset.",
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Directorio de salida. Por defecto PAPER_A_OUT o ./paper_A, como siempre.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    global OUT_DIR
    args = parse_args(argv)
    if args.outdir is not None:
        OUT_DIR = Path(args.outdir)
        OUT_DIR.mkdir(parents=True, exist_ok=True)
    if args.lw_boundary is not None:
        LW_BOUNDARY_OVERRIDES["radiation.surface_lw_boundary"] = args.lw_boundary

    configure_style()
    print(f"Salidas en {OUT_DIR.resolve()}\n")
    if LW_BOUNDARY_OVERRIDES:
        print(f"Frontera LW inferior: {LW_BOUNDARY_OVERRIDES['radiation.surface_lw_boundary']}\n")

    # ---- Fuente A -------------------------------------------------------
    print("FUENTE A -- framework titan_pbl")
    print("-" * 60)
    ctx = load_framework()

    figure1_schematic()

    thresholds = FALLBACK_THRESHOLDS
    if ctx is not None:
        thresholds = compute_thresholds(ctx)
        figure2_heating_profiles(ctx)
        figure3_destabilisation(ctx, thresholds)
        figure4_mixing_time(ctx, thresholds)
        table1_configuration(ctx)
        table2_parameter_sweep(ctx)
        table_thresholds(ctx, thresholds)
    else:
        print("  Figuras 2 a 4 y Tablas 1 a 2 omitidas.")

    # ---- Fuente B -------------------------------------------------------
    print("\nFUENTE B -- archivo de Hartwick et al. (2026)")
    print("-" * 60)
    subdir, tmn, note = ARCHIVE_BASELINES[BASELINE_RUN]
    print(f"  Baseline: {BASELINE_RUN}  ({subdir})")
    print(f"            TMN = {tmn:.4f} K; {note}")
    print( "            Cambiala con TITAN_BASELINE_RUN="
          f"{'|'.join(sorted(ARCHIVE_BASELINES))}")

    if SKIP_ARCHIVE:
        print("  Omitida por SKIP_ARCHIVE=1")
    elif not Path(DATA_DIR).is_dir():
        print(f"  El directorio {DATA_DIR} no existe.")
        print("  Figuras 5 y 6 y Tablas 3 y 4 omitidas.")
    else:
        for name, function, needs_ctx in [
            ("Figura 5", figure5_validation, True),
            ("Figura 6", figure6_timeseries, False),
            ("Tabla 3", table3_validation, True),
            ("Tabla 4", table4_contrasts, False),
        ]:
            try:
                function(ctx) if needs_ctx else function()
            except Exception as exc:  # noqa: BLE001
                print(f"  {name} fallo: {type(exc).__name__}: {exc}")

    print("\n" + "=" * 60)
    print("COMPLETO")
    print("=" * 60)
    for path in sorted(OUT_DIR.iterdir()):
        print(f"  {path.name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
