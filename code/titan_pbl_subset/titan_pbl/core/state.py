"""
Model state container and initialisation.

:class:`ModelState` holds every prognostic and diagnostic variable of the
column at a single instant. It is the object that every physics scheme
receives, reads from and writes to, and it is therefore the second half of the
module contract together with :class:`Tendencies`.

Contract
--------
A physics scheme must:

* read all of its parameters from the configuration object, never from
  literals in its own source;
* read the model variables it needs from a :class:`ModelState`;
* return its contribution as a :class:`Tendencies` object rather than mutating
  the state directly.

Only the core solver applies tendencies to the state. This separation is what
allows several schemes to be composed, swapped or disabled without any of them
knowing about the others.

Reference
---------
Hartwick, Rafkin & Soto (2026), PSJ 7:178.
"""

from __future__ import annotations

import copy as _copy
from dataclasses import dataclass, field, fields

import numpy as np
from numpy.typing import NDArray

from titan_pbl.config.schema import TitanPBLConfig
from titan_pbl.core.grid import ModelGrid
from titan_pbl.io.huygens_data import load_hasi_profile


# ---------------------------------------------------------------------------
# Tendencias
# ---------------------------------------------------------------------------


@dataclass
class Tendencies:
    """Time derivatives contributed by a single physics scheme.

    All fields default to ``None``, meaning "this scheme contributes nothing
    to that variable". The solver treats a missing field as zero, so a scheme
    only has to populate what it actually computes.

    Attributes
    ----------
    atmospheric_temperature
        Temperature tendency at atmospheric midpoints [K s-1].
    methane_mixing_ratio
        Methane mass mixing ratio tendency at atmospheric midpoints [s-1].
    turbulent_kinetic_energy
        TKE tendency at atmospheric midpoints [m2 s-3].
    subsurface_temperature
        Temperature tendency in subsurface layers [K s-1].
    skin_temperature
        Surface skin temperature tendency [K s-1].
    subsurface_methane
        Subsurface methane inventory tendency [kg m-2 s-1]. Unused while
        ``methane.infinite_reservoir`` is true; reserved for v0.2.0.
    """

    atmospheric_temperature: NDArray[np.float64] | None = None
    methane_mixing_ratio: NDArray[np.float64] | None = None
    turbulent_kinetic_energy: NDArray[np.float64] | None = None
    subsurface_temperature: NDArray[np.float64] | None = None
    skin_temperature: float | None = None
    subsurface_methane: NDArray[np.float64] | None = None

    def __add__(self, other: "Tendencies") -> "Tendencies":
        """Combine the contributions of two schemes."""
        if not isinstance(other, Tendencies):
            return NotImplemented
        combined = Tendencies()
        for f in fields(self):
            a = getattr(self, f.name)
            b = getattr(other, f.name)
            if a is None:
                setattr(combined, f.name, _copy.deepcopy(b))
            elif b is None:
                setattr(combined, f.name, _copy.deepcopy(a))
            else:
                setattr(combined, f.name, a + b)
        return combined

    def scaled(self, factor: float) -> "Tendencies":
        """Return a copy with every contribution multiplied by ``factor``."""
        scaled = Tendencies()
        for f in fields(self):
            value = getattr(self, f.name)
            setattr(scaled, f.name, None if value is None else value * factor)
        return scaled

    def is_empty(self) -> bool:
        """True when the object carries no contribution at all."""
        return all(getattr(self, f.name) is None for f in fields(self))


# ---------------------------------------------------------------------------
# Estado del modelo
# ---------------------------------------------------------------------------


@dataclass
class ModelState:
    """Complete instantaneous state of the modelled column.

    Ordering conventions
    --------------------
    Atmospheric arrays are indexed from the surface upward: index 0 is the
    layer in contact with the ground. Subsurface arrays are indexed downward:
    index 0 is the shallowest layer.

    Sign conventions
    ----------------
    Radiative fluxes are positive downward. Turbulent fluxes
    (``sensible_heat_flux``, ``latent_heat_flux``) are positive when energy
    leaves the surface and enters the atmosphere. ``ground_heat_flux`` is
    positive when energy is conducted downward, away from the surface layer.
    These match Equations 1 to 5 of the reference paper, where the surface
    residual is ``Qnet = Qrad - Qg - QH - QLH``.
    """

    # --- Referencias inmutables ------------------------------------------
    grid: ModelGrid
    config: TitanPBLConfig

    # --- Coordenada temporal ---------------------------------------------
    time: float = 0.0
    """Elapsed simulated time since the start of the run [s]."""

    step: int = 0
    """Number of completed time steps."""

    # --- Variables prognosticas atmosfericas ------------------------------
    atmospheric_temperature: NDArray[np.float64] = field(default=None)  # type: ignore[assignment]
    """Temperature at atmospheric midpoints [K]."""

    methane_mixing_ratio: NDArray[np.float64] = field(default=None)  # type: ignore[assignment]
    """Methane mass mixing ratio at atmospheric midpoints [kg kg-1]."""

    turbulent_kinetic_energy: NDArray[np.float64] = field(default=None)  # type: ignore[assignment]
    """Turbulent kinetic energy at atmospheric midpoints [m2 s-2]."""

    zonal_wind: NDArray[np.float64] = field(default=None)  # type: ignore[assignment]
    """Zonal wind component at atmospheric midpoints [m s-1]."""

    meridional_wind: NDArray[np.float64] = field(default=None)  # type: ignore[assignment]
    """Meridional wind component at atmospheric midpoints [m s-1]."""

    # --- Variables prognosticas de superficie y subsuelo -------------------
    skin_temperature: float = 0.0
    """Surface skin temperature Ts [K]."""

    subsurface_temperature: NDArray[np.float64] = field(default=None)  # type: ignore[assignment]
    """Temperature at subsurface layer midpoints [K]."""

    subsurface_methane: NDArray[np.float64] = field(default=None)  # type: ignore[assignment]
    """Liquid methane inventory per subsurface layer [kg m-2].

    Held fixed while ``methane.infinite_reservoir`` is true, which reproduces
    the assumption of the reference paper. Becomes prognostic in v0.2.0.
    """

    # --- Coeficientes de intercambio diagnosticados ------------------------
    friction_velocity: float = 0.0
    """Surface friction velocity u* [m s-1]."""

    exchange_coefficient_heat: float = 0.0
    """Bulk heat exchange coefficient psi_H [W m-2 K-1], Equation 4."""

    exchange_coefficient_moisture: float = 0.0
    """Bulk methane exchange coefficient psi_q [kg m-2 s-1], Equation 12."""

    surface_layer_heat_capacity: float = 0.0
    """Areal heat capacity of the uppermost regolith layer [J m-2 K-1].

    Published by the subsurface scheme and consumed by the surface scheme to
    convert the residual surface energy of Equation 1 into a skin temperature
    tendency. Communicating it through the state rather than importing the
    subsurface scheme is what keeps the two modules independent.
    """

    eddy_diffusivity_heat: NDArray[np.float64] = field(default=None)  # type: ignore[assignment]
    """Eddy diffusivity for heat at layer interfaces [m2 s-1]."""

    eddy_diffusivity_momentum: NDArray[np.float64] = field(default=None)  # type: ignore[assignment]
    """Eddy diffusivity for momentum at layer interfaces [m2 s-1]."""

    # --- Flujos superficiales diagnosticados -------------------------------
    shortwave_down: float = 0.0
    """Downward shortwave flux reaching the surface [W m-2]."""

    longwave_down: float = 0.0
    """Downward longwave flux reaching the surface [W m-2]."""

    longwave_up: float = 0.0
    """Total upward longwave flux at the surface [W m-2]: the emission, plus
    the reflected part when ``radiation.surface_lw_boundary='kirchhoff'``."""

    net_radiation: float = 0.0
    """Net surface radiative flux Qrad [W m-2], Equation 2."""

    ground_heat_flux: float = 0.0
    """Conductive flux into the subsurface Qg [W m-2], Equation 3."""

    sensible_heat_flux: float = 0.0
    """Sensible heat flux QH [W m-2], Equation 4."""

    latent_heat_flux: float = 0.0
    """Latent heat flux QLH [W m-2], Equation 5."""

    methane_flux: float = 0.0
    """Methane evaporative mass flux QCH4 [kg m-2 s-1]."""

    evaporation_depth: float = 0.0
    """Depth below the surface at which the phase change occurs [m].

    Zero for schemes that evaporate at the surface. When it is non-zero the
    latent cooling is a sink inside the regolith rather than a term of the
    surface energy budget, and it reaches the surface only by conduction,
    which is the temporal lag Section 2.4 of the reference paper describes.
    """

    net_surface_energy: float = 0.0
    """Residual surface energy Qnet [W m-2], Equation 1."""

    # --- Diagnosticos derivados --------------------------------------------
    pbl_depth: float = 0.0
    """Diagnosed planetary boundary layer depth [m]."""

    radiative_heating_rate: NDArray[np.float64] = field(default=None)  # type: ignore[assignment]
    """Atmospheric radiative heating rate [K s-1]."""

    # =======================================================================
    # Coordenada temporal derivada
    # =======================================================================

    @property
    def local_time(self) -> float:
        """Local solar time [h], in the interval [0, 24)."""
        sol_length = self.config.planet.sol_length
        elapsed_fraction = (self.time % sol_length) / sol_length
        start = self.config.site.start_local_time
        return float((start + 24.0 * elapsed_fraction) % 24.0)

    @property
    def sol_number(self) -> int:
        """Index of the current sol, starting at 1."""
        return int(self.time // self.config.planet.sol_length) + 1

    @property
    def is_past_spinup(self) -> bool:
        """True once the configured spin-up period has elapsed."""
        spinup_seconds = (
            self.config.simulation.spinup_sols * self.config.planet.sol_length
        )
        return self.time >= spinup_seconds

    # =======================================================================
    # Cantidades termodinamicas derivadas
    # =======================================================================

    @property
    def pressure(self) -> NDArray[np.float64]:
        """Pressure at atmospheric midpoints [Pa]. Read only alias of the grid."""
        return self.grid.atmosphere.pressure_midpoint

    @property
    def air_temperature(self) -> float:
        """Temperature of the lowest atmospheric layer Ta [K], Equation 4."""
        return float(self.atmospheric_temperature[0])

    @property
    def surface_air_temperature_difference(self) -> float:
        """Ts minus Ta [K]. Positive means the surface is warmer than the air."""
        return self.skin_temperature - self.air_temperature

    def potential_temperature(self) -> NDArray[np.float64]:
        """Potential temperature at atmospheric midpoints [K].

        Uses the reference pressure and gas properties from the configuration,
        so the result is directly comparable with the HASI profile processed
        through :meth:`~titan_pbl.io.HASIProfile.potential_temperature`.
        """
        gas = self.config.atmosphere.gas
        p0 = self.config.atmosphere.reference_state.reference_pressure
        return self.atmospheric_temperature * (p0 / self.pressure) ** gas.kappa

    def density(self) -> NDArray[np.float64]:
        """Air density at atmospheric midpoints [kg m-3] from the ideal gas law."""
        r = self.config.atmosphere.gas.specific_gas_constant
        return self.pressure / (r * self.atmospheric_temperature)

    def surface_air_density(self) -> float:
        """Air density in the lowest atmospheric layer [kg m-3]."""
        return float(self.density()[0])

    def geopotential_heights(self) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Interface and midpoint heights [m], integrated hydrostatically."""
        return self.grid.atmosphere.geopotential_heights(
            self.atmospheric_temperature,
            self.config.planet.gravity,
            self.config.atmosphere.gas.specific_gas_constant,
            self.config.site.surface_elevation,
        )

    def energy_balance_residual(self) -> float:
        """Closure error of the surface energy budget [W m-2].

        Equation 1 states ``Qnet = Qrad - Qg - QH - QLH``. This method returns
        the difference between the stored ``net_surface_energy`` and the value
        implied by the individual terms, which must stay at round-off level.
        A growing residual is the earliest symptom of a broken flux sign or a
        scheme writing to the wrong field.
        """
        latent = 0.0 if self.evaporation_depth > 0.0 else self.latent_heat_flux
        implied = (
            self.net_radiation
            - self.ground_heat_flux
            - self.sensible_heat_flux
            - latent
        )
        return float(self.net_surface_energy - implied)

    # =======================================================================
    # Integridad numerica
    # =======================================================================

    def check_finite(self, context: str = "") -> None:
        """Raise if any state variable has become non-finite.

        The error message names the offending field and the offending indices,
        which is what turns a silent NaN cascade into a one line diagnosis.

        Parameters
        ----------
        context
            Optional label identifying where the check was performed, for
            example the name of the scheme that just ran.
        """
        where = f" [{context}]" if context else ""

        arrays = {
            "atmospheric_temperature": self.atmospheric_temperature,
            "methane_mixing_ratio": self.methane_mixing_ratio,
            "turbulent_kinetic_energy": self.turbulent_kinetic_energy,
            "subsurface_temperature": self.subsurface_temperature,
            "eddy_diffusivity_heat": self.eddy_diffusivity_heat,
            "eddy_diffusivity_momentum": self.eddy_diffusivity_momentum,
        }
        for name, values in arrays.items():
            if values is None:
                continue
            bad = ~np.isfinite(values)
            if np.any(bad):
                indices = np.flatnonzero(bad).tolist()
                raise FloatingPointError(
                    f"Non-finite values in '{name}'{where} at step {self.step} "
                    f"(t = {self.time:.1f} s, LT = {self.local_time:.2f} h), "
                    f"indices {indices[:10]}"
                    + (" ..." if len(indices) > 10 else "")
                )

        scalars = {
            "skin_temperature": self.skin_temperature,
            "friction_velocity": self.friction_velocity,
            "sensible_heat_flux": self.sensible_heat_flux,
            "latent_heat_flux": self.latent_heat_flux,
            "ground_heat_flux": self.ground_heat_flux,
            "net_radiation": self.net_radiation,
            "pbl_depth": self.pbl_depth,
        }
        for name, value in scalars.items():
            if not np.isfinite(value):
                raise FloatingPointError(
                    f"Non-finite value in '{name}'{where} at step {self.step} "
                    f"(t = {self.time:.1f} s): {value}"
                )

    def check_physical_bounds(self, context: str = "") -> None:
        """Raise if a variable has left a physically meaningful range.

        This is deliberately stricter than :meth:`check_finite`: an integration
        that produces 40 K air on Titan has not failed numerically, but it has
        certainly failed physically, and catching it early saves hours.
        """
        where = f" [{context}]" if context else ""

        if np.any(self.atmospheric_temperature <= 0.0):
            raise ValueError(
                f"Non-positive atmospheric temperature{where} at step {self.step}."
            )
        if np.any(self.subsurface_temperature <= 0.0):
            raise ValueError(
                f"Non-positive subsurface temperature{where} at step {self.step}."
            )
        if self.skin_temperature <= 0.0:
            raise ValueError(
                f"Non-positive skin temperature{where} at step {self.step}: "
                f"{self.skin_temperature}"
            )
        if np.any(self.methane_mixing_ratio < 0.0):
            raise ValueError(
                f"Negative methane mixing ratio{where} at step {self.step}."
            )
        if np.any(self.turbulent_kinetic_energy < 0.0):
            raise ValueError(f"Negative turbulent kinetic energy{where} at step {self.step}.")

    # =======================================================================
    # Utilidades
    # =======================================================================

    def copy(self) -> "ModelState":
        """Return a deep copy of the state, sharing the grid and configuration.

        The grid and the configuration are immutable, so they are shared by
        reference; every array is duplicated.
        """
        clone = ModelState(grid=self.grid, config=self.config)
        for f in fields(self):
            if f.name in ("grid", "config"):
                continue
            value = getattr(self, f.name)
            setattr(clone, f.name, value.copy() if isinstance(value, np.ndarray) else value)
        return clone

    def scalar_diagnostics(self) -> dict[str, float]:
        """Return the scalar diagnostics of the current step as a flat dictionary."""
        return {
            "time": self.time,
            "step": float(self.step),
            "local_time": self.local_time,
            "sol_number": float(self.sol_number),
            "skin_temperature": self.skin_temperature,
            "air_temperature": self.air_temperature,
            "surface_air_difference": self.surface_air_temperature_difference,
            "shortwave_down": self.shortwave_down,
            "longwave_down": self.longwave_down,
            "longwave_up": self.longwave_up,
            "net_radiation": self.net_radiation,
            "ground_heat_flux": self.ground_heat_flux,
            "sensible_heat_flux": self.sensible_heat_flux,
            "latent_heat_flux": self.latent_heat_flux,
            "methane_flux": self.methane_flux,
            "evaporation_depth": self.evaporation_depth,
            "net_surface_energy": self.net_surface_energy,
            "pbl_depth": self.pbl_depth,
            "friction_velocity": self.friction_velocity,
        }

    def summary(self) -> dict[str, float | int]:
        """Return a compact human readable snapshot of the state."""
        return {
            "step": self.step,
            "time_s": round(self.time, 1),
            "sol": self.sol_number,
            "local_time_h": round(self.local_time, 3),
            "skin_temperature_K": round(self.skin_temperature, 4),
            "air_temperature_K": round(self.air_temperature, 4),
            "pbl_depth_m": round(self.pbl_depth, 2),
            "sensible_heat_flux_W_m2": round(self.sensible_heat_flux, 5),
            "latent_heat_flux_W_m2": round(self.latent_heat_flux, 5),
        }


# ---------------------------------------------------------------------------
# Inicializacion
# ---------------------------------------------------------------------------


def _initial_atmospheric_temperature(
    config: TitanPBLConfig, grid: ModelGrid
) -> NDArray[np.float64]:
    """Build the initial atmospheric temperature profile [K].

    The reference paper initialises the column with the HASI descent profile,
    which is interpolated onto the model pressure levels.
    """
    init = config.atmosphere.initialization
    pressure = grid.atmosphere.pressure_midpoint

    if init.profile_source == "isothermal":
        assert init.isothermal_temperature is not None
        return np.full(grid.n_atmospheric_layers, init.isothermal_temperature)

    if init.profile_source in ("hasi", "custom"):
        path = init.hasi_file if init.profile_source == "hasi" else init.custom_file
        assert path is not None
        profile = load_hasi_profile(path)

        p_min, p_max = float(profile.pressure.min()), float(profile.pressure.max())
        if pressure.min() < p_min or pressure.max() > p_max:
            raise ValueError(
                f"The model pressure range [{pressure.min():.1f}, {pressure.max():.1f}] Pa "
                f"is not fully covered by the observed profile "
                f"[{p_min:.1f}, {p_max:.1f}] Pa. Adjust "
                "atmosphere.vertical_grid.surface_pressure or top_pressure."
            )
        return profile.interpolate_to_pressure(pressure)["temperature"]

    raise ValueError(  # pragma: no cover - Pydantic ya restringe los valores
        f"Unsupported profile_source '{init.profile_source}'."
    )


def _initial_methane_profile(
    config: TitanPBLConfig, grid: ModelGrid
) -> NDArray[np.float64]:
    """Build the initial atmospheric methane mixing ratio profile [kg kg-1].

    Raises
    ------
    NotImplementedError
        When a published methane profile is requested but the corresponding
        digitised data file has not been supplied. The framework will not
        invent an observational profile.
    """
    n = grid.n_atmospheric_layers

    # En una corrida seca no hay metano que transportar.
    if config.methane.scheme == "none":
        return np.zeros(n, dtype=np.float64)

    source = config.methane.initial_atmospheric_profile.source

    if source == "zero":
        return np.zeros(n, dtype=np.float64)

    if source == "uniform":
        # Perfil bien mezclado. Es una simplificacion explicita y documentada,
        # no un perfil observacional inventado: se usa cuando no se dispone del
        # reanalisis GCMS digitalizado. Ver H-11.
        profile = config.methane.initial_atmospheric_profile
        if profile.relative_humidity is not None:
            from titan_pbl.schemes.methane.saturation import MethaneSchemeBase

            class _Probe(MethaneSchemeBase):
                def pore_space_mixing_ratio(self, state):  # pragma: no cover
                    return 0.0

                def evaporative_flux(self, state):  # pragma: no cover
                    return 0.0

            probe = _Probe(config)
            saturated = probe.saturation_mixing_ratio(
                config.subsurface.deep_temperature,
                config.atmosphere.vertical_grid.surface_pressure,
            )
            return np.full(n, profile.relative_humidity * saturated, dtype=np.float64)
        return np.full(n, profile.surface_mixing_ratio, dtype=np.float64)

    if source == "custom":
        path = config.methane.initial_atmospheric_profile.custom_file
        assert path is not None
        table = np.loadtxt(path, dtype=np.float64)
        if table.ndim != 2 or table.shape[1] < 2:
            raise ValueError(
                f"{path} must contain at least two columns: altitude [m] and "
                "methane mass mixing ratio [kg kg-1]."
            )
        altitude, mixing_ratio = table[:, 0], table[:, 1]
        order = np.argsort(altitude)
        temperature = _initial_atmospheric_temperature(config, grid)
        _, z_mid = grid.atmosphere.geopotential_heights(
            temperature,
            config.planet.gravity,
            config.atmosphere.gas.specific_gas_constant,
            config.site.surface_elevation,
        )
        return np.interp(z_mid, altitude[order], mixing_ratio[order])

    raise NotImplementedError(
        f"The '{source}' methane profile requires a digitised data file that is "
        "not bundled with the framework. Provide it through "
        "methane.initial_atmospheric_profile.source='custom' together with "
        "custom_file, or set the source to 'zero' for a methane free "
        "initialisation."
    )


def _initial_subsurface_temperature(
    config: TitanPBLConfig, grid: ModelGrid, skin_temperature: float
) -> NDArray[np.float64]:
    """Build the initial subsurface temperature profile [K].

    The reference paper prescribes a linear fit between the surface skin
    temperature and the fixed deep temperature at the base of the domain.
    """
    mode = config.subsurface.initial_profile_mode
    deep = config.subsurface.deep_temperature

    if mode == "isothermal":
        return np.full(grid.n_subsurface_layers, skin_temperature, dtype=np.float64)

    if mode == "linear":
        depth = grid.subsurface.midpoint_depth
        total = grid.subsurface.total_depth
        weight = depth / total
        return skin_temperature + weight * (deep - skin_temperature)

    raise NotImplementedError(
        f"initial_profile_mode='{mode}' is not implemented; use 'linear' or "
        "'isothermal'."
    )


def initialize_state(config: TitanPBLConfig, grid: ModelGrid) -> ModelState:
    """Build the initial :class:`ModelState` for a simulation.

    The initialisation reproduces Section 2.5 of the reference paper: the
    atmosphere is set from the HASI descent profile, the skin temperature is
    set equal to the temperature of the lowest atmospheric layer so that the
    initial sensible heat flux vanishes, and the subsurface is filled with a
    linear profile between the skin temperature and the fixed deep temperature.

    Parameters
    ----------
    config
        Validated configuration.
    grid
        Model grid built from the same configuration.

    Returns
    -------
    ModelState
        A state ready to be integrated, with every prognostic field populated
        and every diagnostic field zeroed.

    Examples
    --------
    >>> from titan_pbl.config import load_config
    >>> from titan_pbl.core import build_grid, initialize_state
    >>> cfg = load_config("huygens_baseline")
    >>> state = initialize_state(cfg, build_grid(cfg))
    >>> state.local_time
    9.5
    """
    n_atm = grid.n_atmospheric_layers
    n_sub = grid.n_subsurface_layers

    temperature = _initial_atmospheric_temperature(config, grid)

    # Ec. 4: la temperatura de piel iguala la del aire, de modo que el flujo
    # de calor sensible inicial es exactamente cero.
    skin_temperature = float(temperature[0])

    state = ModelState(
        grid=grid,
        config=config,
        time=0.0,
        step=0,
        atmospheric_temperature=temperature,
        methane_mixing_ratio=_initial_methane_profile(config, grid),
        turbulent_kinetic_energy=np.full(
            n_atm, config.turbulence.tke_minimum, dtype=np.float64
        ),
        zonal_wind=np.full(n_atm, config.atmosphere.initialization.zonal_wind),
        meridional_wind=np.full(n_atm, config.atmosphere.initialization.meridional_wind),
        skin_temperature=skin_temperature,
        subsurface_temperature=_initial_subsurface_temperature(
            config, grid, skin_temperature
        ),
        subsurface_methane=np.zeros(n_sub, dtype=np.float64),
        friction_velocity=config.surface_layer.initial_friction_velocity,
        eddy_diffusivity_heat=np.full(
            n_atm + 1, config.turbulence.diffusivity_min, dtype=np.float64
        ),
        eddy_diffusivity_momentum=np.full(
            n_atm + 1, config.turbulence.diffusivity_min, dtype=np.float64
        ),
        radiative_heating_rate=np.zeros(n_atm, dtype=np.float64),
    )

    state.check_finite("initialisation")
    state.check_physical_bounds("initialisation")
    return state


__all__ = [
    "ModelState",
    "Tendencies",
    "initialize_state",
]
