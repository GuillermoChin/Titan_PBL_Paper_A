"""
Universal physical constants (Tier 1).

This module is the ONLY place in the entire framework where a numeric value
may be hardcoded. It contains exclusively universal constants whose values
are fixed by international convention (CODATA 2018) and therefore cannot be
treated as tunable model parameters.

Anything that describes Titan, the Huygens site, the regolith, the numerical
scheme or any other modelling choice belongs in a YAML configuration file
(Tier 2), never here.

References
----------
CODATA 2018 recommended values of the fundamental physical constants.
NIST Standard Reference Database 121.
"""

from __future__ import annotations

from typing import Final

# ---------------------------------------------------------------------------
# Constantes fundamentales (CODATA 2018)
# ---------------------------------------------------------------------------

STEFAN_BOLTZMANN: Final[float] = 5.670374419e-8
"""Stefan-Boltzmann constant [W m-2 K-4]. Exact by SI definition."""

BOLTZMANN: Final[float] = 1.380649e-23
"""Boltzmann constant [J K-1]. Exact by SI definition."""

AVOGADRO: Final[float] = 6.02214076e23
"""Avogadro constant [mol-1]. Exact by SI definition."""

UNIVERSAL_GAS_CONSTANT: Final[float] = 8.314462618
"""Molar gas constant R = N_A * k_B [J mol-1 K-1]."""

PLANCK: Final[float] = 6.62607015e-34
"""Planck constant [J s]. Exact by SI definition."""

SPEED_OF_LIGHT: Final[float] = 299792458.0
"""Speed of light in vacuum [m s-1]. Exact by SI definition."""

GRAVITATIONAL_CONSTANT: Final[float] = 6.67430e-11
"""Newtonian constant of gravitation [m3 kg-1 s-2]."""

# ---------------------------------------------------------------------------
# Constantes astronomicas de referencia
# ---------------------------------------------------------------------------

ASTRONOMICAL_UNIT: Final[float] = 1.495978707e11
"""Astronomical unit [m]. Exact by IAU 2012 definition."""

SOLAR_CONSTANT_1AU: Final[float] = 1361.0
"""Total solar irradiance at 1 AU [W m-2]. IAU/WMO reference value.

Note
----
This is the reference irradiance at Earth's mean orbital distance, provided
only for sanity checks and unit conversions. The insolation actually used by
the radiation scheme is the Tier 2 parameter ``radiation.solar_constant``,
which is prescribed at the top of the model domain for Titan.
"""

# ---------------------------------------------------------------------------
# Conversiones de unidades exactas
# ---------------------------------------------------------------------------

DEG_TO_RAD: Final[float] = 0.017453292519943295
"""Degrees to radians conversion factor [rad deg-1]."""

RAD_TO_DEG: Final[float] = 57.29577951308232
"""Radians to degrees conversion factor [deg rad-1]."""

PASCAL_PER_BAR: Final[float] = 1.0e5
"""Pascals per bar [Pa bar-1]. Exact by definition."""

SECONDS_PER_EARTH_DAY: Final[float] = 86400.0
"""Seconds in one Earth solar day [s]. Exact by definition."""

ZERO_CELSIUS: Final[float] = 273.15
"""Freezing point of water at standard pressure [K]. Exact by definition."""


__all__ = [
    "STEFAN_BOLTZMANN",
    "BOLTZMANN",
    "AVOGADRO",
    "UNIVERSAL_GAS_CONSTANT",
    "PLANCK",
    "SPEED_OF_LIGHT",
    "GRAVITATIONAL_CONSTANT",
    "ASTRONOMICAL_UNIT",
    "SOLAR_CONSTANT_1AU",
    "DEG_TO_RAD",
    "RAD_TO_DEG",
    "PASCAL_PER_BAR",
    "SECONDS_PER_EARTH_DAY",
    "ZERO_CELSIUS",
]
