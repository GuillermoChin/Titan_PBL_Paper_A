"""
Interface for methane evaporation schemes.

A methane scheme computes the evaporative mass flux of methane from the
regolith into the atmosphere, the associated latent heat flux, and the depth
at which the phase change occurs.

Equations owned by this category
--------------------------------
Equation 5 (latent heat flux) plus one of the three parameterisations of
Hartwick et al. (2026):

SPS, saturated pore space
    Equations 12 and 13. The pore space is at saturation for the skin
    temperature and the flux is limited by the porosity and the turbulent
    exchange coefficient.

NSMT, near surface methane table
    Equations 14 to 17. Methane diffuses from a prescribed table depth through
    a dry porous medium, with an effective diffusivity reduced by tortuosity.

MRPS, methane rich pore space
    Equations 18 to 20. Regolith grains are coated with liquid methane and the
    pore space holds a mixture whose abundance is set by a moisture parameter.

Why the evaporation depth is part of the contract
-------------------------------------------------
The three schemes differ in *where* the latent cooling is applied. SPS and
MRPS apply it at the surface; NSMT applies it at the methane table depth, and
the reference paper notes that this introduces a temporal lag because the
cooling must be conducted upward before it reaches the surface energy budget.
:meth:`MethaneScheme.evaporation_depth` makes that difference explicit rather
than hiding it inside each implementation.

Sign convention
---------------
``evaporative_flux`` is positive when methane leaves the regolith and enters
the atmosphere. ``latent_heat_flux`` is correspondingly positive when the
surface loses energy, matching the ``- QLH`` term in Equation 1.
"""

from __future__ import annotations

from abc import abstractmethod
from typing import ClassVar

from titan_pbl.core.state import ModelState
from titan_pbl.interfaces.base import PhysicsScheme


class MethaneScheme(PhysicsScheme):
    """Contract for any subsurface methane evaporation model."""

    category: ClassVar[str] = "methane"

    # -----------------------------------------------------------------------
    # Termodinamica del metano
    # -----------------------------------------------------------------------

    @abstractmethod
    def saturation_mixing_ratio(self, temperature: float, pressure: float) -> float:
        """Saturation mass mixing ratio of methane [kg kg-1].

        Parameters
        ----------
        temperature
            Temperature at which saturation is evaluated [K].
        pressure
            Total pressure [Pa].

        Notes
        -----
        The reference paper writes this as ``qsat(Ts)`` in Equation 13 and
        ``qsat(Tz)`` in Equation 17 but does not state the formulation used.
        The choice is therefore exposed as the Tier 2 parameter
        ``methane.saturation.formulation`` and is a documented source of
        discrepancy during validation.
        """

    @abstractmethod
    def pore_space_mixing_ratio(self, state: ModelState) -> float:
        """Methane mass mixing ratio in the regolith pore space [kg kg-1].

        This is ``qp``, defined by Equation 13 for SPS and Equation 18 for
        MRPS. For NSMT it is the saturation value at the methane table depth.
        """

    @abstractmethod
    def atmospheric_mixing_ratio(self, state: ModelState) -> float:
        """Effective methane mixing ratio the pore space exchanges with [kg kg-1].

        This is ``qa`` for SPS and NSMT, and ``qeff`` of Equation 19 for MRPS,
        which weights the pore space and the lowest atmospheric layer by the
        porosity.
        """

    # -----------------------------------------------------------------------
    # Flujos
    # -----------------------------------------------------------------------

    @abstractmethod
    def evaporative_flux(self, state: ModelState) -> float:
        """Methane evaporative mass flux QCH4 [kg m-2 s-1].

        Implements Equation 12, 17 or 20 depending on the scheme. Positive
        upward, into the atmosphere.
        """

    def latent_heat_flux(self, state: ModelState) -> float:
        """Latent heat flux QLH [W m-2], Equation 5.

        The default implementation multiplies the evaporative flux by the
        configured latent heat of vaporisation, which is how the reference
        paper defines it for all three schemes. A scheme only needs to
        override this if it partitions the phase change differently.
        """
        # Ec. 5 de Hartwick et al. 2026: QLH = XLV * QCH4
        return self.config.methane.latent_heat_vaporization * self.evaporative_flux(
            state
        )

    # -----------------------------------------------------------------------
    # Geometria del cambio de fase
    # -----------------------------------------------------------------------

    @abstractmethod
    def evaporation_depth(self, state: ModelState) -> float:
        """Depth below the surface at which the phase change occurs [m].

        Zero for schemes that evaporate at the surface, such as SPS and MRPS.
        Equal to the configured methane table depth for NSMT, where the latent
        cooling is applied at depth and reaches the surface only by conduction.
        """

    def effective_porosity(self, state: ModelState) -> float:
        """Porosity used by this scheme [dimensionless].

        Defaults to :attr:`TitanPBLConfig.effective_porosity`, which returns
        the scheme specific override when present and the bulk regolith
        porosity otherwise.
        """
        return self.config.effective_porosity

    # -----------------------------------------------------------------------
    # Balance de masa, reservado para v0.2.0
    # -----------------------------------------------------------------------

    def reservoir_depletion_rate(self, state: ModelState) -> float:
        """Rate at which the subsurface methane inventory is depleted [kg m-2 s-1].

        Returns zero while ``methane.infinite_reservoir`` is true, which
        reproduces the assumption of the reference paper that the subsurface
        reservoir is static and therefore infinite. Schemes that track a
        finite inventory override this in v0.2.0.
        """
        if self.config.methane.infinite_reservoir:
            return 0.0
        raise NotImplementedError(
            f"{type(self).__name__} does not implement a finite methane budget. "
            "Set methane.infinite_reservoir to true, or use a scheme that "
            "tracks the subsurface inventory."
        )


__all__ = ["MethaneScheme"]
