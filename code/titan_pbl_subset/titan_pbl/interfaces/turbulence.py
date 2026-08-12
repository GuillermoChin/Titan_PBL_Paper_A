"""
Interface for boundary layer turbulence schemes.

A turbulence scheme owns the vertical mixing of the atmospheric column: the
turbulent kinetic energy budget, the eddy diffusivities derived from it, the
resulting tendencies of temperature and methane, and the diagnosis of the
boundary layer depth.

Equations owned by this category
--------------------------------
The Mellor-Yamada-Janjic closure described in Section 2.1 of Hartwick et al.
(2026), in which eddy diffusivities are diagnosed from the subgrid turbulent
kinetic energy evaluated through the Richardson number, and the boundary layer
depth is the highest level at which the TKE exceeds a critical value.

Why the PBL diagnosis lives here
--------------------------------
The reference paper defines the PBL depth in terms of the TKE field, which is
a prognostic variable of this scheme. Placing the diagnosis anywhere else
would force another module to know the internals of the closure. Alternative
diagnoses that use only the potential temperature profile belong in
:mod:`titan_pbl.diagnostics` instead, where they can be compared against this
one.
"""

from __future__ import annotations

from abc import abstractmethod
from typing import ClassVar

import numpy as np
from numpy.typing import NDArray

from titan_pbl.core.state import ModelState
from titan_pbl.interfaces.base import PhysicsScheme


class TurbulenceScheme(PhysicsScheme):
    """Contract for any boundary layer turbulence closure."""

    category: ClassVar[str] = "turbulence"

    # -----------------------------------------------------------------------
    # Estabilidad
    # -----------------------------------------------------------------------

    @abstractmethod
    def brunt_vaisala_frequency_squared(
        self, state: ModelState
    ) -> NDArray[np.float64]:
        """Squared buoyancy frequency at layer interfaces [s-2].

        .. math:: N^2 = \\frac{g}{\\theta} \\frac{\\partial \\theta}{\\partial z}

        Negative values indicate a statically unstable layer. Shape
        ``(n_layers + 1,)``.
        """

    @abstractmethod
    def richardson_number(self, state: ModelState) -> NDArray[np.float64]:
        """Gradient Richardson number at layer interfaces [dimensionless].

        In the one dimensional configuration the horizontal winds are
        invariant, so the shear term is prescribed rather than resolved and
        the variation is driven by buoyancy alone, as stated in Section 2.1.
        Shape ``(n_layers + 1,)``.
        """

    # -----------------------------------------------------------------------
    # Energia cinetica turbulenta y difusividades
    # -----------------------------------------------------------------------

    @abstractmethod
    def mixing_length(self, state: ModelState) -> NDArray[np.float64]:
        """Turbulent mixing length at layer interfaces [m]."""

    @abstractmethod
    def tke_tendency(self, state: ModelState, dt: float) -> NDArray[np.float64]:
        """Turbulent kinetic energy tendency at midpoints [m2 s-3].

        The budget combines buoyant production or consumption, shear
        production, vertical transport and dissipation.
        """

    @abstractmethod
    def eddy_diffusivities(
        self, state: ModelState
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Eddy diffusivities for heat and momentum at interfaces [m2 s-1].

        Both arrays have shape ``(n_layers + 1,)``. The returned values must
        already be clipped to the configured bounds, so that no downstream
        module has to guard against a diffusivity of zero or infinity.
        """

    # -----------------------------------------------------------------------
    # Mezcla vertical
    # -----------------------------------------------------------------------

    @abstractmethod
    def vertical_diffusion(
        self,
        state: ModelState,
        field: NDArray[np.float64],
        surface_flux: float,
        dt: float,
    ) -> NDArray[np.float64]:
        """Apply vertical eddy diffusion to an arbitrary scalar field.

        Parameters
        ----------
        state
            Current model state.
        field
            Scalar at layer midpoints, for example temperature or methane
            mixing ratio.
        surface_flux
            Flux of the scalar across the lower boundary, in the units of the
            field multiplied by kg m-2 s-1. This is the lower boundary
            condition supplied by the surface layer scheme.
        dt
            Time step [s].

        Returns
        -------
        ndarray
            The diffused field, same shape as ``field``.

        Notes
        -----
        Being generic in the field is what lets the same routine mix heat and
        methane without duplicating the solver, which is essential for keeping
        the two consistent when a methane scheme is active.
        """

    # -----------------------------------------------------------------------
    # Diagnostico de la PBL
    # -----------------------------------------------------------------------

    @abstractmethod
    def diagnose_pbl_depth(self, state: ModelState) -> float:
        """Diagnose the depth of the planetary boundary layer [m].

        Section 2.1 of the reference paper defines it as the height of the
        highest atmospheric level at which the turbulent kinetic energy
        exceeds the configured critical value ``turbulence.tke_critical``.

        An implementation must return zero, not a negative number or a NaN,
        when no level exceeds the threshold, which is the normal night time
        situation.
        """


__all__ = ["TurbulenceScheme"]
