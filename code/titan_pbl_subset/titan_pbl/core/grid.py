"""
Vertical grid construction for the atmosphere and the subsurface.

The framework uses two independent one-dimensional grids that share a common
interface at the surface:

Atmosphere
    A pressure based coordinate running from the surface upward. Index 0 is
    the lowest layer, index ``n_layers - 1`` the topmost. Interfaces are
    stored separately from midpoints, because fluxes live on interfaces while
    thermodynamic variables live on midpoints.

Subsurface
    A depth based coordinate running downward from the surface. Index 0 is the
    shallowest layer. Depths are positive downward.

Important
---------
Hartwick et al. (2026) prescribe 56 atmospheric layers spanning roughly 1.45
to 0.66 bar but do not state how those layers are distributed in the vertical.
The distribution is therefore exposed as the Tier 2 parameters
``atmosphere.vertical_grid.spacing_type`` and ``spacing_exponent``, and is a
documented source of potential discrepancy during validation. Layer thickness
near the ground directly controls how sharply the PBL top can be resolved.

Reference
---------
Hartwick, Rafkin & Soto (2026), PSJ 7:178, Sections 2.1 and 2.3, Table 1.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from titan_pbl.config.schema import TitanPBLConfig


# ---------------------------------------------------------------------------
# Grilla atmosferica
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AtmosphericGrid:
    """Pressure based vertical grid for the atmospheric column.

    Attributes
    ----------
    n_layers
        Number of atmospheric layers.
    eta_interface
        Normalised vertical coordinate at layer interfaces, shape
        ``(n_layers + 1,)``. Equals 1 at the surface and 0 at the model top.
        Ordered from the surface upward, so the array is decreasing.
    eta_midpoint
        Normalised coordinate at layer midpoints, shape ``(n_layers,)``.
    pressure_interface
        Pressure at layer interfaces [Pa], shape ``(n_layers + 1,)``.
    pressure_midpoint
        Pressure at layer midpoints [Pa], shape ``(n_layers,)``.
    surface_pressure
        Pressure at the surface [Pa].
    top_pressure
        Pressure at the model top [Pa].

    Notes
    -----
    Index 0 is always the layer in contact with the surface. Pressure
    decreases with index.
    """

    n_layers: int
    eta_interface: NDArray[np.float64]
    eta_midpoint: NDArray[np.float64]
    pressure_interface: NDArray[np.float64]
    pressure_midpoint: NDArray[np.float64]
    surface_pressure: float
    top_pressure: float
    spacing_type: str
    spacing_exponent: float

    # -----------------------------------------------------------------------
    # Propiedades derivadas puramente geometricas
    # -----------------------------------------------------------------------

    @property
    def delta_pressure(self) -> NDArray[np.float64]:
        """Pressure thickness of each layer [Pa], always positive.

        Defined as ``p_interface[k] - p_interface[k+1]``.
        """
        return self.pressure_interface[:-1] - self.pressure_interface[1:]

    def layer_mass(self, gravity: float) -> NDArray[np.float64]:
        """Mass per unit area contained in each layer [kg m-2].

        Follows from hydrostatic balance, ``dm = dp / g``.

        Parameters
        ----------
        gravity
            Surface gravity [m s-2]. A Tier 2 parameter; never assumed.
        """
        return self.delta_pressure / gravity

    def total_column_mass(self, gravity: float) -> float:
        """Total atmospheric mass per unit area in the domain [kg m-2]."""
        return float(np.sum(self.layer_mass(gravity)))

    # -----------------------------------------------------------------------
    # Altura geopotencial: requiere un perfil de temperatura
    # -----------------------------------------------------------------------

    def geopotential_heights(
        self,
        temperature_midpoint: NDArray[np.float64],
        gravity: float,
        gas_constant: float,
        surface_elevation: float = 0.0,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Integrate the hydrostatic equation to obtain layer heights.

        The hypsometric relation for an ideal gas gives, between two pressure
        levels within a layer of mean temperature ``T``,

        .. math:: \\Delta z = \\frac{R T}{g} \\ln\\!\\left(\\frac{p_1}{p_2}\\right)

        Parameters
        ----------
        temperature_midpoint
            Temperature at layer midpoints [K], shape ``(n_layers,)``.
        gravity
            Surface gravity [m s-2].
        gas_constant
            Specific gas constant of the atmosphere [J kg-1 K-1].
        surface_elevation
            Height assigned to the surface interface [m].

        Returns
        -------
        z_interface, z_midpoint
            Heights above the datum [m], shapes ``(n_layers + 1,)`` and
            ``(n_layers,)``.

        Notes
        -----
        Gravity is held constant across the domain. Over 20 km on Titan this
        introduces an error below 1.6 percent in ``g`` and is consistent with
        the shallow atmosphere approximation used by the reference model. The
        temperature used is the dry temperature; a virtual temperature
        correction for methane is deliberately not applied at v0.1.0 because
        the reference model does not apply one either.
        """
        temperature = np.asarray(temperature_midpoint, dtype=np.float64)
        if temperature.shape != (self.n_layers,):
            raise ValueError(
                f"temperature_midpoint must have shape ({self.n_layers},), "
                f"got {temperature.shape}."
            )
        if np.any(temperature <= 0.0):
            raise ValueError("temperature_midpoint contains non-positive values.")

        scale = gas_constant * temperature / gravity  # [m], altura de escala local

        # Espesor de cada capa entre sus interfaces.
        log_ratio_full = np.log(
            self.pressure_interface[:-1] / self.pressure_interface[1:]
        )
        dz_layer = scale * log_ratio_full

        z_interface = np.empty(self.n_layers + 1, dtype=np.float64)
        z_interface[0] = surface_elevation
        z_interface[1:] = surface_elevation + np.cumsum(dz_layer)

        # Altura del punto medio: se integra desde la interfaz inferior hasta
        # la presion del punto medio, no se promedia geometricamente.
        log_ratio_half = np.log(self.pressure_interface[:-1] / self.pressure_midpoint)
        z_midpoint = z_interface[:-1] + scale * log_ratio_half

        return z_interface, z_midpoint

    def layer_thickness(
        self,
        temperature_midpoint: NDArray[np.float64],
        gravity: float,
        gas_constant: float,
    ) -> NDArray[np.float64]:
        """Geometric thickness of each layer [m]."""
        z_interface, _ = self.geopotential_heights(
            temperature_midpoint, gravity, gas_constant
        )
        return np.diff(z_interface)

    def summary(self) -> dict[str, float | int | str]:
        """Return a compact description of the atmospheric grid."""
        dp = self.delta_pressure
        return {
            "n_layers": self.n_layers,
            "spacing_type": self.spacing_type,
            "spacing_exponent": self.spacing_exponent,
            "surface_pressure_Pa": round(self.surface_pressure, 1),
            "top_pressure_Pa": round(self.top_pressure, 1),
            "thinnest_layer_Pa": round(float(dp.min()), 3),
            "thickest_layer_Pa": round(float(dp.max()), 3),
            "thickness_ratio": round(float(dp.max() / dp.min()), 2),
        }


# ---------------------------------------------------------------------------
# Grilla del subsuelo
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SubsurfaceGrid:
    """Depth based vertical grid for the regolith column.

    Attributes
    ----------
    n_layers
        Number of subsurface layers.
    thickness
        Thickness of each layer [m], shape ``(n_layers,)``.
    midpoint_depth
        Depth of each layer midpoint below the surface [m].
    interface_depth
        Depth of each layer interface [m], shape ``(n_layers + 1,)``.
        ``interface_depth[0]`` is zero, the surface itself.

    Notes
    -----
    Index 0 is the shallowest layer. Depth is positive downward. The base of
    the deepest layer is where the fixed deep temperature boundary condition
    is applied.
    """

    n_layers: int
    thickness: NDArray[np.float64]
    midpoint_depth: NDArray[np.float64]
    interface_depth: NDArray[np.float64]

    @property
    def total_depth(self) -> float:
        """Depth of the lower boundary of the domain [m]."""
        return float(self.interface_depth[-1])

    @property
    def midpoint_spacing(self) -> NDArray[np.float64]:
        """Distance between consecutive layer midpoints [m], shape ``(n_layers - 1,)``.

        This is the length scale that appears in the conductive flux between
        adjacent layers, Equation 11 of the reference paper.
        """
        return np.diff(self.midpoint_depth)

    def layer_index_at_depth(self, depth: float) -> int:
        """Return the index of the layer containing a given depth.

        Parameters
        ----------
        depth
            Depth below the surface [m], positive downward.

        Raises
        ------
        ValueError
            If the depth lies outside the subsurface domain.
        """
        if depth < 0.0 or depth > self.total_depth:
            raise ValueError(
                f"Depth {depth} m lies outside the subsurface domain "
                f"[0, {self.total_depth}] m."
            )
        # searchsorted devuelve el primer indice tal que depth <= interface.
        index = int(np.searchsorted(self.interface_depth, depth, side="left")) - 1
        return int(np.clip(index, 0, self.n_layers - 1))

    def diurnal_skin_depth(
        self,
        thermal_conductivity: float,
        volumetric_heat_capacity: float,
        period: float,
    ) -> float:
        """Thermal skin depth for a forcing of a given period [m].

        .. math:: d = \\sqrt{\\frac{\\kappa P}{\\pi C_v}}

        Parameters
        ----------
        thermal_conductivity
            Effective conductivity [W m-1 K-1].
        volumetric_heat_capacity
            Effective volumetric heat capacity [J m-3 K-1].
        period
            Forcing period [s]; the sol length for the diurnal skin depth,
            the orbital period for the annual one.

        Notes
        -----
        The grid is only adequate if the shallowest layers resolve the diurnal
        skin depth and the deepest boundary lies below the annual one. This
        method exists so that condition can be checked quantitatively rather
        than assumed.
        """
        diffusivity = thermal_conductivity / volumetric_heat_capacity
        return float(np.sqrt(diffusivity * period / np.pi))

    def summary(self) -> dict[str, float | int]:
        """Return a compact description of the subsurface grid."""
        return {
            "n_layers": self.n_layers,
            "total_depth_m": round(self.total_depth, 4),
            "first_layer_thickness_m": round(float(self.thickness[0]), 5),
            "last_layer_thickness_m": round(float(self.thickness[-1]), 4),
            "first_midpoint_m": round(float(self.midpoint_depth[0]), 5),
            "last_midpoint_m": round(float(self.midpoint_depth[-1]), 4),
        }


# ---------------------------------------------------------------------------
# Contenedor conjunto
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ModelGrid:
    """The complete vertical discretisation of a single column."""

    atmosphere: AtmosphericGrid
    subsurface: SubsurfaceGrid
    metadata: dict[str, str] = field(default_factory=dict)

    @property
    def n_atmospheric_layers(self) -> int:
        return self.atmosphere.n_layers

    @property
    def n_subsurface_layers(self) -> int:
        return self.subsurface.n_layers

    def summary(self) -> dict[str, dict]:
        return {
            "atmosphere": self.atmosphere.summary(),
            "subsurface": self.subsurface.summary(),
        }


# ---------------------------------------------------------------------------
# Constructores
# ---------------------------------------------------------------------------


def _eta_interfaces(n_layers: int, spacing_type: str, exponent: float) -> NDArray[np.float64]:
    """Build the normalised vertical coordinate at layer interfaces.

    The returned array runs from 1.0 at the surface to 0.0 at the model top
    and has ``n_layers + 1`` entries.

    Strategies
    ----------
    sigma
        Power law stretching, ``eta = 1 - s**exponent`` with ``s = k / n``.
        An exponent greater than one concentrates layers near the surface;
        an exponent of one gives uniform spacing in sigma.
    log_pressure
        Uniform spacing in the logarithm of the normalised pressure, which
        yields layers of roughly constant thickness in scale heights. The
        exponent applies the same power law stretching to the log coordinate.
    geometric
        Layer thicknesses in eta follow a geometric progression whose common
        ratio is the exponent, mirroring the subsurface construction.
    """
    if n_layers < 2:
        raise ValueError(f"n_layers must be at least 2, got {n_layers}.")
    if exponent <= 0.0:
        raise ValueError(f"spacing_exponent must be positive, got {exponent}.")

    s = np.arange(n_layers + 1, dtype=np.float64) / n_layers

    if spacing_type == "sigma":
        eta = 1.0 - s**exponent

    elif spacing_type == "log_pressure":
        # Espaciamiento uniforme en log(eta) es singular en eta = 0, asi que
        # se estira el indice y luego se mapea de forma monotona.
        stretched = s**exponent
        eta = 1.0 - stretched

    elif spacing_type == "geometric":
        if np.isclose(exponent, 1.0):
            weights = np.ones(n_layers, dtype=np.float64)
        else:
            # Con razon r y n capas, la mas gruesa es r**(n-1) veces la mas
            # fina. Si ese contraste supera la precision util del flotante de
            # doble precision, las capas inferiores se pierden en el redondeo
            # y la coordenada deja de ser monotona. Se detecta aqui para dar
            # un diagnostico util en lugar de un fallo opaco mas abajo.
            contrast = exponent ** (n_layers - 1)
            if contrast > 1.0e12:
                raise ValueError(
                    f"spacing_type='geometric' with exponent={exponent} and "
                    f"n_layers={n_layers} implies a thickness contrast of "
                    f"{contrast:.3g} between the thinnest and thickest layer, "
                    "which cannot be represented reliably in double precision. "
                    f"Use an exponent below {1.0e12 ** (1.0 / (n_layers - 1)):.4f} "
                    "for this number of layers."
                )
            weights = exponent ** np.arange(n_layers, dtype=np.float64)
        cumulative = np.concatenate([[0.0], np.cumsum(weights)])
        eta = 1.0 - cumulative / cumulative[-1]

    else:  # pragma: no cover - el esquema Pydantic ya restringe los valores
        raise ValueError(f"Unknown spacing_type '{spacing_type}'.")

    # Blindaje numerico contra errores de redondeo en los extremos.
    eta[0] = 1.0
    eta[-1] = 0.0

    if not np.all(np.diff(eta) < 0.0):
        raise ValueError(
            f"spacing_type='{spacing_type}' with exponent={exponent} produced a "
            "non monotonic vertical coordinate. Choose a different exponent."
        )
    return eta


def build_atmospheric_grid(config: TitanPBLConfig) -> AtmosphericGrid:
    """Construct the atmospheric grid from a validated configuration.

    Every geometric choice is read from ``config``; nothing is hardcoded.
    """
    grid_cfg = config.atmosphere.vertical_grid
    p_surface = grid_cfg.surface_pressure
    p_top = grid_cfg.top_pressure

    eta_interface = _eta_interfaces(
        grid_cfg.n_layers, grid_cfg.spacing_type, grid_cfg.spacing_exponent
    )
    eta_midpoint = 0.5 * (eta_interface[:-1] + eta_interface[1:])

    # p = p_top + eta * (p_surface - p_top)
    span = p_surface - p_top
    pressure_interface = p_top + eta_interface * span
    pressure_midpoint = p_top + eta_midpoint * span

    return AtmosphericGrid(
        n_layers=grid_cfg.n_layers,
        eta_interface=eta_interface,
        eta_midpoint=eta_midpoint,
        pressure_interface=pressure_interface,
        pressure_midpoint=pressure_midpoint,
        surface_pressure=p_surface,
        top_pressure=p_top,
        spacing_type=grid_cfg.spacing_type,
        spacing_exponent=grid_cfg.spacing_exponent,
    )


def build_subsurface_grid(config: TitanPBLConfig) -> SubsurfaceGrid:
    """Construct the subsurface grid from a validated configuration.

    Layer thicknesses come from :attr:`TitanPBLConfig.subsurface_layer_thicknesses`,
    which honours both the geometric and the explicit construction modes.
    """
    thickness = np.asarray(config.subsurface_layer_thicknesses, dtype=np.float64)

    interface_depth = np.concatenate([[0.0], np.cumsum(thickness)])
    midpoint_depth = 0.5 * (interface_depth[:-1] + interface_depth[1:])

    return SubsurfaceGrid(
        n_layers=int(thickness.size),
        thickness=thickness,
        midpoint_depth=midpoint_depth,
        interface_depth=interface_depth,
    )


def build_grid(config: TitanPBLConfig) -> ModelGrid:
    """Construct the complete model grid from a validated configuration.

    Examples
    --------
    >>> from titan_pbl.config import load_config
    >>> from titan_pbl.core import build_grid
    >>> grid = build_grid(load_config("huygens_baseline"))
    >>> grid.n_atmospheric_layers, grid.n_subsurface_layers
    (56, 10)
    """
    return ModelGrid(
        atmosphere=build_atmospheric_grid(config),
        subsurface=build_subsurface_grid(config),
        metadata={
            "config_name": config.metadata.name,
            "config_version": config.metadata.config_version,
        },
    )


__all__ = [
    "AtmosphericGrid",
    "SubsurfaceGrid",
    "ModelGrid",
    "build_atmospheric_grid",
    "build_subsurface_grid",
    "build_grid",
]
