"""
Interface for subsurface thermal schemes.

A subsurface scheme owns the regolith column: it computes the effective
thermal properties of each layer from the constituent mixture, conducts heat
between layers, and reports the conductive flux exchanged with the surface.

Equations owned by this category
--------------------------------
Equations 6 to 9 (effective thermal properties from the porosity and methane
saturation), Equation 10 (temperature tendency from flux divergence) and
Equation 11 (conductive flux) of Hartwick et al. (2026). Equation 3, the
ground heat flux entering the surface energy budget, is the value returned by
:meth:`SubsurfaceScheme.surface_conductive_flux`.

Sign convention
---------------
``surface_conductive_flux`` is positive downward, so a positive value removes
energy from the surface layer. This matches the ``- Qg`` term in Equation 1.
"""

from __future__ import annotations

from abc import abstractmethod
from typing import ClassVar

import numpy as np
from numpy.typing import NDArray

from titan_pbl.core.state import ModelState
from titan_pbl.interfaces.base import PhysicsScheme


class SubsurfaceScheme(PhysicsScheme):
    """Contract for any regolith thermal model."""

    category: ClassVar[str] = "subsurface"

    # -----------------------------------------------------------------------
    # Propiedades termicas efectivas
    # -----------------------------------------------------------------------

    @abstractmethod
    def effective_conductivity(self, state: ModelState) -> NDArray[np.float64]:
        """Effective thermal conductivity of each subsurface layer [W m-1 K-1].

        Implements the mixing rule of Equation 8, weighting the solid, liquid
        methane and pore air conductivities by the porosity and the pore space
        saturation.
        """

    @abstractmethod
    def effective_volumetric_heat_capacity(
        self, state: ModelState
    ) -> NDArray[np.float64]:
        """Effective volumetric heat capacity of each layer [J m-3 K-1].

        Implements Equations 6, 7 and 9: the specific heat and the density are
        mixed separately and then multiplied.
        """

    def effective_density(self, state: ModelState) -> NDArray[np.float64]:
        """Effective bulk density of each layer [kg m-3], Equation 7.

        Provided as a separate hook because the density enters the volumetric
        heat capacity and is also a reported diagnostic in its own right.
        """
        raise NotImplementedError

    def thermal_inertia(self, state: ModelState) -> NDArray[np.float64]:
        """Thermal inertia of each layer [J m-2 K-1 s-0.5].

        .. math:: I = \\sqrt{\\kappa \\, C_v}

        The default implementation composes the two required hooks, so a
        scheme gets this for free. Table 3 of the reference paper quotes
        354.44 for the dry baseline.
        """
        conductivity = self.effective_conductivity(state)
        heat_capacity = self.effective_volumetric_heat_capacity(state)
        return np.sqrt(conductivity * heat_capacity)

    # -----------------------------------------------------------------------
    # Flujos
    # -----------------------------------------------------------------------

    @abstractmethod
    def surface_conductive_flux(self, state: ModelState) -> float:
        """Conductive heat flux from the surface into the regolith [W m-2].

        This is Qg of Equation 3, positive downward. The solver writes the
        returned value into ``state.ground_heat_flux``.
        """

    @abstractmethod
    def interface_fluxes(self, state: ModelState) -> NDArray[np.float64]:
        """Conductive flux at every internal layer interface [W m-2].

        Shape ``(n_layers + 1,)``, positive downward. Index 0 is the flux
        across the surface and the last entry is the flux across the fixed
        deep boundary.
        """


__all__ = ["SubsurfaceScheme"]
