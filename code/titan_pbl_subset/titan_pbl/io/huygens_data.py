"""
Reader for the Huygens Atmospheric Structure Instrument (HASI) descent profile.

The HASI Level 4 atmospheric profile distributed through the NASA Planetary
Data System is a plain fixed-width table with five columns. The column
identification below was verified against the ideal gas law: reconstructing
density from the pressure and temperature columns reproduces the tabulated
density to within 1-2 percent, the residual being explained by the N2/CH4
mixture having a slightly higher mean molecular mass than pure N2.

Column layout
-------------
======  ========================  ==========
Index   Quantity                  Unit
======  ========================  ==========
0       Mission elapsed time      ms
1       Altitude above surface    m
2       Atmospheric pressure      Pa
3       Atmospheric temperature   K
4       Atmospheric density       kg m-3
======  ========================  ==========

The descent is recorded from the top down, so the raw file is ordered by
*decreasing* altitude. This reader returns profiles sorted by *increasing*
altitude, which is the convention used throughout the framework.

References
----------
Fulchignoni, M. et al. (2005), Nature 438, 785.
Hartwick, Rafkin & Soto (2026), PSJ 7:178, Figure 4.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray

from titan_pbl.paths import resolve_data_path

# Indices de columna en el archivo del PDS.
_COL_TIME = 0
_COL_ALTITUDE = 1
_COL_PRESSURE = 2
_COL_TEMPERATURE = 3
_COL_DENSITY = 4
_N_COLUMNS = 5

# Rangos de validacion fisica. Son deliberadamente generosos: solo sirven para
# detectar un archivo corrupto o un orden de columnas inesperado, no para
# filtrar ciencia.
_VALID_RANGES: dict[str, tuple[float, float]] = {
    "altitude": (-100.0, 2.0e6),  # [m]
    "pressure": (1.0e-3, 5.0e5),  # [Pa]
    "temperature": (50.0, 300.0),  # [K]
    "density": (1.0e-8, 50.0),  # [kg m-3]
}


@dataclass(frozen=True)
class HASIProfile:
    """An immutable HASI descent profile, sorted by increasing altitude.

    Attributes
    ----------
    time
        Mission elapsed time [s].
    altitude
        Altitude above the surface [m].
    pressure
        Atmospheric pressure [Pa].
    temperature
        Atmospheric temperature [K].
    density
        Atmospheric density [kg m-3].
    source_file
        Path the profile was read from.
    """

    time: NDArray[np.float64]
    altitude: NDArray[np.float64]
    pressure: NDArray[np.float64]
    temperature: NDArray[np.float64]
    density: NDArray[np.float64]
    source_file: Path

    # -----------------------------------------------------------------------
    # Propiedades basicas
    # -----------------------------------------------------------------------

    def __len__(self) -> int:
        return int(self.altitude.size)

    @property
    def surface_temperature(self) -> float:
        """Temperature at the lowest sampled altitude [K]."""
        return float(self.temperature[0])

    @property
    def surface_pressure(self) -> float:
        """Pressure at the lowest sampled altitude [Pa]."""
        return float(self.pressure[0])

    @property
    def altitude_range(self) -> tuple[float, float]:
        """Minimum and maximum sampled altitude [m]."""
        return float(self.altitude[0]), float(self.altitude[-1])

    # -----------------------------------------------------------------------
    # Cantidades derivadas
    # -----------------------------------------------------------------------

    def potential_temperature(
        self,
        reference_pressure: float,
        gas_constant: float,
        specific_heat_cp: float,
    ) -> NDArray[np.float64]:
        """Compute potential temperature from the measured profile.

        Uses the standard Poisson relation

        .. math:: \\theta = T \\left(\\frac{P_0}{P}\\right)^{R/c_p}

        Parameters
        ----------
        reference_pressure
            Reference pressure P0 [Pa]. Hartwick et al. (2026) use 1.5 bar.
        gas_constant
            Specific gas constant of the atmosphere R [J kg-1 K-1].
        specific_heat_cp
            Specific heat at constant pressure [J kg-1 K-1].

        Returns
        -------
        ndarray
            Potential temperature [K], same ordering as ``altitude``.

        Notes
        -----
        All three arguments are Tier 2 parameters and must be supplied by the
        caller from the configuration object. They are deliberately not given
        default values so that no reference state is ever silently assumed.
        """
        kappa = gas_constant / specific_heat_cp
        return self.temperature * (reference_pressure / self.pressure) ** kappa

    def derived_density(self, gas_constant: float) -> NDArray[np.float64]:
        """Reconstruct density from pressure and temperature via the ideal gas law.

        Comparing the result against :attr:`density` is the recommended sanity
        check that the file columns were interpreted correctly.
        """
        return self.pressure / (gas_constant * self.temperature)

    def scale_height(self, gravity: float, gas_constant: float) -> NDArray[np.float64]:
        """Local pressure scale height H = R*T/g [m]."""
        return gas_constant * self.temperature / gravity

    # -----------------------------------------------------------------------
    # Remuestreo
    # -----------------------------------------------------------------------

    def interpolate_to_altitude(
        self, target_altitude: NDArray[np.float64]
    ) -> dict[str, NDArray[np.float64]]:
        """Interpolate the profile onto an arbitrary altitude grid.

        Parameters
        ----------
        target_altitude
            Altitudes to interpolate onto [m], any ordering.

        Returns
        -------
        dict
            Keys ``pressure``, ``temperature`` and ``density``.

        Notes
        -----
        Pressure is interpolated in log space because it varies exponentially
        with height; temperature and density are interpolated linearly.
        Targets outside the sampled range are clamped to the nearest endpoint,
        so the caller should verify the requested range beforehand.
        """
        target = np.asarray(target_altitude, dtype=np.float64)
        log_pressure = np.interp(target, self.altitude, np.log(self.pressure))
        return {
            "pressure": np.exp(log_pressure),
            "temperature": np.interp(target, self.altitude, self.temperature),
            "density": np.interp(target, self.altitude, self.density),
        }

    def interpolate_to_pressure(
        self, target_pressure: NDArray[np.float64]
    ) -> dict[str, NDArray[np.float64]]:
        """Interpolate the profile onto a pressure grid.

        This is the routine used to initialise the model atmosphere, whose
        vertical coordinate is pressure based.
        """
        target = np.asarray(target_pressure, dtype=np.float64)
        # np.interp exige abscisas crecientes: la presion decrece con la
        # altitud, asi que invertimos el orden.
        p_ascending = self.pressure[::-1]
        order = np.argsort(p_ascending)
        p_sorted = p_ascending[order]
        return {
            "altitude": np.interp(target, p_sorted, self.altitude[::-1][order]),
            "temperature": np.interp(target, p_sorted, self.temperature[::-1][order]),
            "density": np.interp(target, p_sorted, self.density[::-1][order]),
        }

    def summary(self) -> dict[str, float | int | str]:
        """Return a compact dictionary describing the loaded profile."""
        z_min, z_max = self.altitude_range
        return {
            "n_points": len(self),
            "altitude_min_m": round(z_min, 2),
            "altitude_max_m": round(z_max, 2),
            "surface_temperature_K": round(self.surface_temperature, 4),
            "surface_pressure_Pa": round(self.surface_pressure, 1),
            "temperature_min_K": round(float(self.temperature.min()), 4),
            "temperature_max_K": round(float(self.temperature.max()), 4),
            "source_file": self.source_file.name,
        }


# ---------------------------------------------------------------------------
# Lectura
# ---------------------------------------------------------------------------


def load_hasi_profile(
    path: str | Path,
    max_altitude: float | None = None,
    validate: bool = True,
) -> HASIProfile:
    """Read a HASI Level 4 atmospheric descent profile.

    Parameters
    ----------
    path
        Path to the ``HASI_L4_ATMO_PROFILE_DESCEN`` table file. Relative paths
        are resolved with :func:`titan_pbl.paths.resolve_data_path`, so a
        preset loads the same file whatever the working directory.
    max_altitude
        Optional altitude ceiling [m]. Points above it are discarded, which is
        useful because the model domain only extends to about 20 km while the
        instrument sampled up to roughly 147 km.
    validate
        Run physical plausibility checks on the parsed columns.

    Returns
    -------
    HASIProfile
        Profile sorted by increasing altitude.

    Raises
    ------
    FileNotFoundError
        If the file does not exist.
    ValueError
        If the file cannot be parsed or fails validation.
    """
    # La resolucion se hace aqui, en el unico punto por el que pasa toda
    # lectura del perfil, para que ni el estado ni los tests dependan del
    # directorio de trabajo.
    source = resolve_data_path(path)

    # El archivo del PDS es de ancho fijo, sin encabezado y sin comentarios.
    # np.loadtxt maneja correctamente el espaciado variable.
    try:
        raw = np.loadtxt(source, dtype=np.float64)
    except Exception as exc:  # pragma: no cover - depende del archivo
        raise ValueError(f"Could not parse {source} as a numeric table: {exc}") from exc

    if raw.ndim != 2 or raw.shape[1] != _N_COLUMNS:
        raise ValueError(
            f"Expected a table with {_N_COLUMNS} columns, got shape {raw.shape}. "
            "The file may use a different HASI product format."
        )

    time_s = raw[:, _COL_TIME] / 1000.0  # ms -> s
    altitude = raw[:, _COL_ALTITUDE]
    pressure = raw[:, _COL_PRESSURE]
    temperature = raw[:, _COL_TEMPERATURE]
    density = raw[:, _COL_DENSITY]

    if validate:
        _validate_columns(altitude, pressure, temperature, density, source)

    # Ordenamos por altitud creciente (el descenso viene de arriba hacia abajo).
    order = np.argsort(altitude)
    time_s = time_s[order]
    altitude = altitude[order]
    pressure = pressure[order]
    temperature = temperature[order]
    density = density[order]

    if max_altitude is not None:
        keep = altitude <= max_altitude
        if not np.any(keep):
            raise ValueError(
                f"max_altitude={max_altitude} m removes every sample; "
                f"the profile spans {altitude.min():.1f}-{altitude.max():.1f} m."
            )
        time_s = time_s[keep]
        altitude = altitude[keep]
        pressure = pressure[keep]
        temperature = temperature[keep]
        density = density[keep]

    return HASIProfile(
        time=time_s,
        altitude=altitude,
        pressure=pressure,
        temperature=temperature,
        density=density,
        source_file=source,
    )


def _validate_columns(
    altitude: NDArray[np.float64],
    pressure: NDArray[np.float64],
    temperature: NDArray[np.float64],
    density: NDArray[np.float64],
    source: Path,
) -> None:
    """Check that each parsed column falls within a physically plausible range."""
    columns = {
        "altitude": altitude,
        "pressure": pressure,
        "temperature": temperature,
        "density": density,
    }
    for name, values in columns.items():
        if not np.all(np.isfinite(values)):
            raise ValueError(f"Column '{name}' in {source.name} contains non-finite values.")
        low, high = _VALID_RANGES[name]
        vmin, vmax = float(values.min()), float(values.max())
        if vmin < low or vmax > high:
            raise ValueError(
                f"Column '{name}' in {source.name} spans {vmin:.4g} to {vmax:.4g}, "
                f"outside the plausible range [{low:.4g}, {high:.4g}]. "
                "The column order of this product may differ from the expected layout."
            )

    # La presion debe crecer monotonicamente conforme baja la altitud.
    order = np.argsort(altitude)
    p_sorted = pressure[order]
    if not np.all(np.diff(p_sorted) <= 0.0):
        n_bad = int(np.sum(np.diff(p_sorted) > 0.0))
        raise ValueError(
            f"Pressure is not monotonically decreasing with altitude in "
            f"{source.name} ({n_bad} violations). Column identification may be wrong."
        )


def check_ideal_gas_consistency(
    profile: HASIProfile, gas_constant: float
) -> dict[str, float]:
    """Quantify how well the profile satisfies the ideal gas law.

    Reconstructs density from the measured pressure and temperature and
    compares it against the tabulated density. A systematic bias of a couple
    of percent is expected for Titan, because the atmosphere is a N2/CH4
    mixture while ``gas_constant`` is typically the pure N2 value.

    Parameters
    ----------
    profile
        Profile to check.
    gas_constant
        Specific gas constant used for the reconstruction [J kg-1 K-1].

    Returns
    -------
    dict
        Relative error statistics, expressed as fractions (not percentages).
    """
    reconstructed = profile.derived_density(gas_constant)
    relative_error = (reconstructed - profile.density) / profile.density
    return {
        "mean_relative_error": float(np.mean(relative_error)),
        "median_relative_error": float(np.median(relative_error)),
        "max_abs_relative_error": float(np.max(np.abs(relative_error))),
        "std_relative_error": float(np.std(relative_error)),
    }


__all__ = [
    "HASIProfile",
    "load_hasi_profile",
    "check_ideal_gas_consistency",
]
