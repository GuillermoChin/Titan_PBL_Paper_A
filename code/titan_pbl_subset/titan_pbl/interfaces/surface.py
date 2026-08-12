"""
Interface for surface layer schemes.

A surface layer scheme owns the exchange between the ground and the lowest
atmospheric layer: the friction velocity, the bulk transfer coefficients, the
sensible heat flux, and the closure of the surface energy budget that evolves
the skin temperature.

Equations owned by this category
--------------------------------
Equation 1 (surface energy budget) and Equation 4 (sensible heat flux) of
Hartwick et al. (2026), together with the Monin-Obukhov similarity relations
that produce the exchange coefficients and the Beljaars gustiness correction.

Sign convention
---------------
``sensible_heat_flux`` is positive when energy leaves the surface and enters
the atmosphere, matching the ``- QH`` term in Equation 1. A negative value
therefore denotes a stable surface layer warming the ground.
"""

from __future__ import annotations

from abc import abstractmethod
from typing import ClassVar

from titan_pbl.core.state import ModelState
from titan_pbl.interfaces.base import PhysicsScheme


class SurfaceLayerScheme(PhysicsScheme):
    """Contract for any surface layer exchange model."""

    category: ClassVar[str] = "surface"

    # -----------------------------------------------------------------------
    # Estabilidad y velocidad de friccion
    # -----------------------------------------------------------------------

    @abstractmethod
    def bulk_richardson_number(self, state: ModelState) -> float:
        """Bulk Richardson number of the surface layer [dimensionless].

        Negative under unstable (convective) conditions, positive under stable
        (inversion) conditions. This is the quantity that selects which
        stability functions apply.
        """

    @abstractmethod
    def friction_velocity(self, state: ModelState) -> float:
        """Surface friction velocity u* [m s-1].

        With the zonal wind set to zero, as in the reference configuration,
        the mechanical contribution vanishes and the value is set by the
        Beljaars convective velocity scale under unstable conditions and by
        the configured floor otherwise.
        """

    # -----------------------------------------------------------------------
    # Coeficientes de intercambio
    # -----------------------------------------------------------------------

    @abstractmethod
    def heat_exchange_coefficient(self, state: ModelState) -> float:
        """Bulk heat exchange coefficient psi_H [W m-2 K-1], Equation 4.

        Larger under convective conditions and smaller under inversions, as
        described in Section 2.1 of the reference paper.
        """

    @abstractmethod
    def moisture_exchange_coefficient(self, state: ModelState) -> float:
        """Bulk methane exchange coefficient psi_q [kg m-2 s-1], Equation 12.

        Required by every methane scheme. A dry configuration may still call
        it, so an implementation must return a finite value rather than
        raising when no methane scheme is active.
        """

    # -----------------------------------------------------------------------
    # Flujos y cierre del balance
    # -----------------------------------------------------------------------

    @abstractmethod
    def sensible_heat_flux(self, state: ModelState) -> float:
        """Sensible heat flux QH [W m-2], Equation 4.

        Equals ``psi_H * (Ts - Ta)``, positive upward.
        """

    @abstractmethod
    def surface_energy_budget(self, state: ModelState) -> float:
        """Residual surface energy Qnet [W m-2], Equation 1.

        Equals ``Qrad - Qg - QH - QLH``, using the flux values already stored
        in the state by the radiation, subsurface and methane schemes. This
        method must therefore run after those three.
        """

    @abstractmethod
    def skin_temperature_tendency(self, state: ModelState) -> float:
        """Rate of change of the surface skin temperature [K s-1].

        The residual energy divides by the heat capacity of the uppermost
        regolith layer, so an implementation needs the effective volumetric
        heat capacity supplied by the subsurface scheme through the state.
        """


__all__ = ["SurfaceLayerScheme"]
