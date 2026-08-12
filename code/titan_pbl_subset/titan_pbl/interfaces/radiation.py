"""
Interface for radiative transfer schemes.

A radiation scheme computes the shortwave and longwave fluxes through the
atmospheric column, the net radiative flux at the surface, and the resulting
atmospheric heating rate.

Equations owned by this category
--------------------------------
Equation 2 of Hartwick et al. (2026), the surface radiative balance, together
with the vertical transfer that produces the downwelling fluxes reaching the
surface and the heating rate in each layer.

Sign convention
---------------
All radiative fluxes are positive downward. ``net_surface_flux`` therefore
returns Qrad of Equation 2, which is positive when the surface gains energy.
"""

from __future__ import annotations

from abc import abstractmethod
from typing import ClassVar

import numpy as np
from numpy.typing import NDArray

from titan_pbl.core.state import ModelState
from titan_pbl.interfaces.base import PhysicsScheme


class RadiationScheme(PhysicsScheme):
    """Contract for any radiative transfer model."""

    category: ClassVar[str] = "radiation"

    # -----------------------------------------------------------------------
    # Geometria solar
    # -----------------------------------------------------------------------

    @abstractmethod
    def cosine_solar_zenith_angle(self, state: ModelState) -> float:
        """Cosine of the solar zenith angle [dimensionless].

        Returns zero at night rather than a negative value, so the caller may
        multiply by it directly. Computed from the latitude, the solar
        longitude and the local time, all of which are Tier 2 parameters.
        """

    @abstractmethod
    def top_of_atmosphere_insolation(self, state: ModelState) -> float:
        """Downward shortwave flux at the top of the model domain [W m-2].

        The product of the configured solar constant and the cosine of the
        solar zenith angle.
        """

    # -----------------------------------------------------------------------
    # Flujos verticales
    # -----------------------------------------------------------------------

    @abstractmethod
    def shortwave_fluxes(
        self, state: ModelState
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Downward and upward shortwave fluxes at layer interfaces [W m-2].

        Both arrays have shape ``(n_layers + 1,)``. Index 0 is the surface.
        """

    @abstractmethod
    def longwave_fluxes(
        self, state: ModelState
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Downward and upward longwave fluxes at layer interfaces [W m-2].

        Both arrays have shape ``(n_layers + 1,)``. Index 0 is the surface.
        """

    # -----------------------------------------------------------------------
    # Superficie y calentamiento atmosferico
    # -----------------------------------------------------------------------

    @abstractmethod
    def net_surface_flux(self, state: ModelState) -> float:
        """Net radiative flux absorbed by the surface Qrad [W m-2], Equation 2.

        Equals ``(1 - albedo) * SW_down + LW_down - emissivity * sigma * Ts**4``,
        with the albedo and emissivity read from the configuration and the
        Stefan-Boltzmann constant taken from
        :mod:`titan_pbl.config.physical_constants`.
        """

    @abstractmethod
    def heating_rate(self, state: ModelState) -> NDArray[np.float64]:
        """Radiative heating rate at atmospheric midpoints [K s-1].

        Obtained from the divergence of the net radiative flux divided by the
        heat capacity of each layer. Negative values denote radiative cooling.
        """

    def surface_emission(self, state: ModelState) -> float:
        """Upward longwave emission from the surface [W m-2].

        Broken out separately because it is a reported diagnostic and because
        the surface energy balance and the atmospheric transfer both need it.
        """
        raise NotImplementedError


__all__ = ["RadiationScheme"]
