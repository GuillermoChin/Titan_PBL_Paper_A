"""
Broadband gray radiative transfer.

Implements the scheme of Chatain et al. (2022), Appendix A (arXiv:2210.03278),
which is the scheme Hartwick et al. (2026) cite but do not reproduce. Shortwave
and longwave are treated as two independent broadband channels.

Equations implemented
---------------------
============  =========================================================
Ec. A1        Insolation at the top of the atmosphere with eccentricity
Ec. A2        Cosine of the solar zenith angle
Ec. A3        Visible optical depth profile
Ec. A5, A6    Two-stream equations for the diffuse intensity
Ec. A7, A8    Backscattering coefficients, Henyey-Greenstein truncated
Ec. A9, A10   Shortwave boundary conditions
Ec. A11, A12  Total shortwave fluxes
Ec. A13       Infrared optical depth profile
Ec. A14       Planck emission of the atmosphere and the surface
Ec. A15, A16  Upward and downward infrared integration
Ec. 2         Surface radiative balance of Hartwick et al. (2026)
============  =========================================================

Parameter discrepancies between the two papers
----------------------------------------------
Chatain et al. tuned the shortwave parameters against the Huygens DISR net
flux profile of Tomasko et al. (2008), so their values carry direct
observational support. Hartwick et al. state that their values follow Chatain,
but four of them differ:

==================  ==============  ==============
Parameter           Hartwick 2026   Chatain 2022
==================  ==============  ==============
Visible optical     10              8
depth
Surface emissivity  0.95            0.9
Downward IR at      2.0             1.55
model top
Surface albedo      0.30            0.3
==================  ==============  ==============

The defaults follow Hartwick, because replicating that paper is the goal of
v0.1.0. The Chatain values are available in the ``chatain2022_radiation``
preset so that the difference can be quantified rather than assumed harmless.

Analytic solution of the two-stream problem
-------------------------------------------
Writing Equations A5 and A6 as a linear system in the diffuse intensities,

.. math::
    \\frac{dI^+}{d\\tau} = \\gamma_1 I^+ - \\gamma_2 I^-
                          - \\frac{X^+}{\\bar\\mu} e^{-\\tau/\\mu_0}, \\qquad
    \\frac{dI^-}{d\\tau} = \\gamma_2 I^+ - \\gamma_1 I^-
                          + \\frac{X^-}{\\bar\\mu} e^{-\\tau/\\mu_0}

with :math:`\\gamma_1 = [1 - \\omega_0 (1 - b)] / \\bar\\mu` and
:math:`\\gamma_2 = \\omega_0 b / \\bar\\mu`, the homogeneous solutions vary as
:math:`e^{\\pm\\lambda\\tau}` with :math:`\\lambda = \\sqrt{\\gamma_1^2 -
\\gamma_2^2}`. The particular solution is proportional to
:math:`e^{-\\tau/\\mu_0}`. The two integration constants follow from the
boundary conditions A9 and A10. The result is exact for the truncated phase
function, so no vertical iteration is needed for the shortwave.

Sign convention
---------------
All fluxes are positive downward.

Surface longwave boundary
-------------------------
``radiation.surface_lw_boundary`` selects how the surface treats the downward
longwave flux. ``"full_absorption"`` reproduces Equation 2 of Hartwick et al.
(2026): the surface absorbs all of it and emits eps*sigma*Ts^4. ``"kirchhoff"``
makes the surface absorb eps*F_down and reflect (1 - eps)*F_down, so that an
isothermal column is in radiative equilibrium. Every preset declares
``"full_absorption"`` explicitly. See H-27.
"""

from __future__ import annotations

from typing import ClassVar

import numpy as np
from numpy.typing import NDArray

from titan_pbl.config.physical_constants import DEG_TO_RAD, STEFAN_BOLTZMANN
from titan_pbl.config.schema import TitanPBLConfig
from titan_pbl.core.state import ModelState, Tendencies
from titan_pbl.interfaces.radiation import RadiationScheme


class GrayRadiation(RadiationScheme):
    """Two-channel gray radiative transfer for Titan.

    Examples
    --------
    >>> from titan_pbl.config import load_config
    >>> from titan_pbl.core import build_grid, initialize_state
    >>> from titan_pbl.schemes.radiation import GrayRadiation
    >>> cfg = load_config("huygens_baseline")
    >>> state = initialize_state(cfg, build_grid(cfg))
    >>> scheme = GrayRadiation(cfg)
    >>> scheme.initialize(cfg, state)
    >>> 0.0 <= scheme.cosine_solar_zenith_angle(state) <= 1.0
    True
    """

    scheme_name: ClassVar[str] = "gray"

    def __init__(self, config: TitanPBLConfig) -> None:
        super().__init__(config)
        self._tau_visible: NDArray[np.float64] | None = None
        self._tau_ir: NDArray[np.float64] | None = None

    # =======================================================================
    # Ciclo de vida
    # =======================================================================

    def initialize(self, config: TitanPBLConfig, state: ModelState) -> None:
        """Precompute the optical depth profiles, which depend only on pressure."""
        self._tau_visible = self._visible_optical_depth(state)
        self._tau_ir = self._infrared_optical_depth(state)
        state.net_radiation = self.net_surface_flux(state)
        self._initialized = True

    # =======================================================================
    # Ecuaciones A1 y A2: geometria solar
    # =======================================================================

    def top_of_atmosphere_insolation(self, state: ModelState) -> float:
        """Downward shortwave flux at the top of the domain [W m-2].

        Combines the eccentricity corrected solar constant of Equation A1 with
        the cosine of the solar zenith angle of Equation A2.
        """
        return self.solar_constant_at_season(state) * self.cosine_solar_zenith_angle(
            state
        )

    def solar_constant_at_season(self, state: ModelState) -> float:
        """Solar constant corrected for orbital distance [W m-2], Equation A1.

        .. math::
            S_{TOA}(L_s) = S_0
            \\left(\\frac{1 + e \\cos(L_s - L_{s0})}{1 - e^2}\\right)^2
        """
        s0 = self.config.radiation.solar_constant
        if not self.config.radiation.apply_eccentricity_correction:
            return s0
        # Ec. A1 de Chatain et al. 2022
        e = self.config.planet.orbital_eccentricity
        ls = self.config.site.solar_longitude_ls * DEG_TO_RAD
        ls0 = self.config.planet.perihelion_solar_longitude * DEG_TO_RAD
        return s0 * ((1.0 + e * np.cos(ls - ls0)) / (1.0 - e**2)) ** 2

    def solar_declination(self) -> float:
        """Solar declination [rad], Equation A2.

        .. math:: \\delta = \\epsilon \\sin(L_s)
        """
        # Ec. A2 de Chatain et al. 2022
        obliquity = self.config.planet.obliquity * DEG_TO_RAD
        ls = self.config.site.solar_longitude_ls * DEG_TO_RAD
        return float(obliquity * np.sin(ls))

    def hour_angle(self, state: ModelState) -> float:
        """Hour angle [rad], Equation A2.

        .. math:: h = \\pi (1 - \\mathrm{local\\_hour} / 12)
        """
        # Ec. A2 de Chatain et al. 2022
        return float(np.pi * (1.0 - state.local_time / 12.0))

    def cosine_solar_zenith_angle(self, state: ModelState) -> float:
        """Cosine of the solar zenith angle, clipped at zero for night, Equation A2."""
        # Ec. A2 de Chatain et al. 2022
        latitude = self.config.site.latitude * DEG_TO_RAD
        declination = self.solar_declination()
        mu0 = np.sin(latitude) * np.sin(declination) + np.cos(latitude) * np.cos(
            declination
        ) * np.cos(self.hour_angle(state))
        return float(max(0.0, mu0))

    # =======================================================================
    # Ecuaciones A3 y A13: profundidades opticas
    # =======================================================================

    def _visible_optical_depth(self, state: ModelState) -> NDArray[np.float64]:
        """Visible optical depth at layer interfaces, Equation A3.

        .. math:: \\tau(z) = \\tau_0 \\left(\\frac{P(z)}{P_0}\\right)^\\gamma
        """
        # Ec. A3 de Chatain et al. 2022
        radiation = self.config.radiation
        grid = state.grid.atmosphere
        ratio = grid.pressure_interface / grid.surface_pressure
        return radiation.optical_depth_visible * ratio**radiation.visible_scaling_exponent

    def _infrared_optical_depth(self, state: ModelState) -> NDArray[np.float64]:
        """Infrared optical depth at layer interfaces, Equation A13.

        .. math::
            \\tau_{IR} = \\tau_{IR}^{surf}
            \\left[f \\frac{P}{P_0} + (1 - f)\\left(\\frac{P}{P_0}\\right)^2\\right]

        The linear term represents a well mixed absorber and the quadratic term
        collision induced absorption.
        """
        # Ec. A13 de Chatain et al. 2022
        radiation = self.config.radiation
        grid = state.grid.atmosphere
        ratio = grid.pressure_interface / grid.surface_pressure
        f = radiation.ir_weighting_coefficient
        return radiation.optical_depth_ir * (f * ratio + (1.0 - f) * ratio**2)

    def visible_optical_depth(self, state: ModelState) -> NDArray[np.float64]:
        """Visible optical depth at layer interfaces [dimensionless]."""
        if self._tau_visible is None:
            return self._visible_optical_depth(state)
        return self._tau_visible

    def infrared_optical_depth(self, state: ModelState) -> NDArray[np.float64]:
        """Infrared optical depth at layer interfaces [dimensionless]."""
        if self._tau_ir is None:
            return self._infrared_optical_depth(state)
        return self._tau_ir

    # =======================================================================
    # Ecuaciones A5 a A12: onda corta
    # =======================================================================

    def _backscattering_coefficients(self, mu0: float) -> tuple[float, float]:
        """Backscattering coefficients b(mu0) and b, Equation A8.

        .. math::
            b(\\mu) = \\tfrac{1}{2}(1 - 3 g \\mu \\bar\\mu), \\qquad
            b = \\tfrac{1}{2}(1 - 3 g \\bar\\mu^2)

        These follow from truncating the Legendre expansion of the
        Henyey-Greenstein phase function at the second term, which is what
        makes the two-stream problem analytically solvable.
        """
        # Ec. A8 de Chatain et al. 2022
        g = self.config.radiation.asymmetry_factor
        mu_bar = self.config.radiation.mean_cosine
        b_mu0 = 0.5 * (1.0 - 3.0 * g * mu0 * mu_bar)
        b = 0.5 * (1.0 - 3.0 * g * mu_bar**2)
        return float(b_mu0), float(b)

    def _two_stream_intensities(
        self, state: ModelState
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Solve Equations A5 and A6 analytically for the diffuse intensities.

        Returns the upward and downward diffuse intensities at every layer
        interface, ordered from the surface upward.
        """
        radiation = self.config.radiation
        mu_bar = radiation.mean_cosine
        omega0 = radiation.single_scattering_albedo
        albedo = radiation.surface_albedo

        mu0 = self.cosine_solar_zenith_angle(state)
        s_toa = self.solar_constant_at_season(state)
        tau = self.visible_optical_depth(state)
        tau_surface = float(tau[0])

        n = tau.size
        if mu0 <= 0.0 or s_toa <= 0.0:
            return np.zeros(n), np.zeros(n)

        b_mu0, b = self._backscattering_coefficients(mu0)

        # Coeficientes del sistema lineal.
        gamma1 = (1.0 - omega0 * (1.0 - b)) / mu_bar
        gamma2 = omega0 * b / mu_bar
        lam = np.sqrt(max(gamma1**2 - gamma2**2, 1.0e-30))

        # Terminos fuente de las Ecs. A5 y A6.
        x_plus = omega0 * s_toa / (2.0 * np.pi) * b_mu0
        x_minus = omega0 * s_toa / (2.0 * np.pi) * (1.0 - b_mu0)

        # Solucion particular proporcional a exp(-tau/mu0).
        denominator = 1.0 / mu0**2 - lam**2
        if abs(denominator) < 1.0e-12:
            # Resonancia entre el haz directo y un modo homogeneo: se
            # perturba mu0 infinitesimalmente para evitar la singularidad.
            mu0 *= 1.0 + 1.0e-6
            denominator = 1.0 / mu0**2 - lam**2

        a_plus = (
            -(x_plus / mu_bar) * (gamma1 - 1.0 / mu0) - gamma2 * (x_minus / mu_bar)
        ) / denominator
        a_minus = (
            -(gamma1 + 1.0 / mu0) * (x_minus / mu_bar) - gamma2 * (x_plus / mu_bar)
        ) / denominator

        # Solucion homogenea: autovectores (gamma2, gamma1 -/+ lambda).
        exp_top = 1.0  # exp(lambda * 0), tope del dominio
        beam_surface = np.exp(-tau_surface / mu0)
        exp_lam_s = np.exp(lam * tau_surface)
        exp_mlam_s = np.exp(-lam * tau_surface)

        # Ec. A10: I^-(tau = 0) = 0
        row1 = np.array([gamma1 - lam, gamma1 + lam])
        rhs1 = -a_minus * exp_top

        # Ec. A9: I^+(tau_s) = alpha [ I^-(tau_s) + mu0 S exp(-tau_s/mu0)/(2 pi mu_bar) ]
        direct_term = albedo * mu0 * s_toa * beam_surface / (2.0 * np.pi * mu_bar)
        row2 = np.array(
            [
                gamma2 * exp_lam_s - albedo * (gamma1 - lam) * exp_lam_s,
                gamma2 * exp_mlam_s - albedo * (gamma1 + lam) * exp_mlam_s,
            ]
        )
        rhs2 = albedo * a_minus * beam_surface - a_plus * beam_surface + direct_term

        constants = np.linalg.solve(np.vstack([row1, row2]), np.array([rhs1, rhs2]))
        c1, c2 = float(constants[0]), float(constants[1])

        exp_lam = np.exp(np.clip(lam * tau, -700.0, 700.0))
        exp_mlam = np.exp(np.clip(-lam * tau, -700.0, 700.0))
        beam = np.exp(-tau / mu0)

        intensity_up = gamma2 * (c1 * exp_lam + c2 * exp_mlam) + a_plus * beam
        intensity_down = (
            c1 * (gamma1 - lam) * exp_lam
            + c2 * (gamma1 + lam) * exp_mlam
            + a_minus * beam
        )
        return intensity_down, intensity_up

    def shortwave_fluxes(
        self, state: ModelState
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Downward and upward shortwave fluxes at interfaces [W m-2].

        Equations A11 and A12: the downward flux adds the attenuated direct
        beam to the diffuse component, while the upward flux is purely diffuse.
        """
        radiation = self.config.radiation
        mu0 = self.cosine_solar_zenith_angle(state)
        s_toa = self.solar_constant_at_season(state)
        tau = self.visible_optical_depth(state)
        n = tau.size

        if mu0 <= 0.0:
            return np.zeros(n), np.zeros(n)

        if radiation.solar_scheme == "none":
            return np.zeros(n), np.zeros(n)

        if radiation.solar_scheme == "beer_lambert":
            # Extincion pura sin dispersion, util como referencia analitica.
            downward = mu0 * s_toa * np.exp(-tau / mu0)
            upward = np.zeros(n)
            upward[:] = radiation.surface_albedo * downward[0]
            return downward, upward

        mu_bar = radiation.mean_cosine
        intensity_down, intensity_up = self._two_stream_intensities(state)

        # Ecs. A11 y A12 de Chatain et al. 2022
        upward = 2.0 * np.pi * mu_bar * intensity_up
        downward = 2.0 * np.pi * mu_bar * intensity_down + mu0 * s_toa * np.exp(
            -tau / mu0
        )
        return downward, upward

    # =======================================================================
    # Ecuaciones A14 a A16: onda larga
    # =======================================================================

    def surface_emission(self, state: ModelState) -> float:
        """Upward longwave emission from the surface [W m-2], Equation A14."""
        # Ec. A14 de Chatain et al. 2022
        emissivity = self.config.radiation.surface_emissivity
        return float(emissivity * STEFAN_BOLTZMANN * state.skin_temperature**4)

    def longwave_fluxes(
        self, state: ModelState
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        """Downward and upward longwave fluxes at interfaces [W m-2].

        Integrates Equations A15 and A16 layer by layer: the downward flux
        starts from the prescribed flux at the model top and the upward flux
        from the surface. At the surface the upward flux is the emission
        ``eps*sigma*Ts^4``, plus the reflected ``(1 - eps)*F_down(0)`` when
        ``radiation.surface_lw_boundary`` is ``"kirchhoff"`` (see H-27).
        """
        radiation = self.config.radiation
        tau = self.infrared_optical_depth(state)
        n = tau.size

        if radiation.ir_scheme == "none":
            return np.zeros(n), np.zeros(n)

        # Ec. A14: emision de Planck de cada capa.
        emission = STEFAN_BOLTZMANN * state.atmospheric_temperature**4

        # El espesor optico decrece hacia arriba, de modo que la diferencia
        # entre interfaces contiguas es positiva por construccion.
        d_tau = tau[:-1] - tau[1:]
        transmission = np.exp(-d_tau)

        upward = np.empty(n, dtype=np.float64)
        downward = np.empty(n, dtype=np.float64)

        # Ec. A16: integracion descendente desde el tope del modelo.
        # Chatain fija F^-(0) = 0; Hartwick prescribe un flujo I_IR no nulo
        # que representa la atmosfera por encima del dominio. Sin dispersion
        # el flujo descendente no depende del ascendente, de modo que se
        # integra primero: la frontera de Kirchhoff necesita F^-(0).
        downward[-1] = radiation.ir_top_flux
        for k in range(state.grid.n_atmospheric_layers - 1, -1, -1):
            downward[k] = downward[k + 1] * transmission[k] + emission[k] * (
                1.0 - transmission[k]
            )

        # Ec. A15: integracion ascendente desde la superficie. Con
        # 'full_absorption' la superficie solo emite; con 'kirchhoff' refleja
        # ademas la fraccion (1 - eps) del flujo descendente (H-27).
        upward[0] = self.surface_emission(state)
        if radiation.surface_lw_boundary == "kirchhoff":
            upward[0] = upward[0] + self._reflected_longwave(float(downward[0]))
        for k in range(state.grid.n_atmospheric_layers):
            upward[k + 1] = upward[k] * transmission[k] + emission[k] * (
                1.0 - transmission[k]
            )

        return downward, upward

    def _reflected_longwave(self, downward_surface: float) -> float:
        """Reflected part (1 - eps) * F_down(0) of the downward longwave [W m-2]."""
        emissivity = self.config.radiation.surface_emissivity
        return float((1.0 - emissivity) * downward_surface)

    def surface_reflected_longwave(self, state: ModelState) -> float:
        """Downward longwave reflected by the surface [W m-2].

        Equals ``(1 - eps) * F_down(0)`` when ``radiation.surface_lw_boundary``
        is ``"kirchhoff"`` and zero when it is ``"full_absorption"``, where the
        surface absorbs all of the downward longwave. Provided as a method
        rather than a :meth:`diagnostics` key so that existing diagnostic
        tables keep their columns.
        """
        if self.config.radiation.surface_lw_boundary != "kirchhoff":
            return 0.0
        downward, _ = self.longwave_fluxes(state)
        return self._reflected_longwave(float(downward[0]))

    # =======================================================================
    # Ecuacion 2: balance radiativo en superficie
    # =======================================================================

    def net_surface_flux(self, state: ModelState) -> float:
        """Net radiative flux absorbed by the surface Qrad [W m-2], Equation 2.

        With ``radiation.surface_lw_boundary = "full_absorption"`` (Equation 2
        of Hartwick et al. 2026, the default of every preset):

        .. math::
            Q_{rad} = (1 - \\alpha) SW^\\downarrow + LW^\\downarrow
                      - \\epsilon \\sigma T_s^4

        With ``"kirchhoff"`` the surface absorbs only the fraction epsilon of
        the downward longwave and reflects the rest:

        .. math::
            Q_{rad} = (1 - \\alpha) SW^\\downarrow + \\epsilon LW^\\downarrow
                      - \\epsilon \\sigma T_s^4
        """
        radiation = self.config.radiation
        shortwave_down, _ = self.shortwave_fluxes(state)
        longwave_down, _ = self.longwave_fluxes(state)

        absorbed_shortwave = (1.0 - radiation.surface_albedo) * float(shortwave_down[0])
        if radiation.surface_lw_boundary == "kirchhoff":
            # Kirchhoff: absorbe eps LW_down y refleja (1 - eps) LW_down (H-27).
            absorbed_longwave = radiation.surface_emissivity * float(longwave_down[0])
            return absorbed_shortwave + absorbed_longwave - self.surface_emission(state)

        # Ec. 2 de Hartwick et al. 2026
        return (
            absorbed_shortwave + float(longwave_down[0]) - self.surface_emission(state)
        )

    # =======================================================================
    # Calentamiento atmosferico
    # =======================================================================

    def net_flux(self, state: ModelState) -> NDArray[np.float64]:
        """Net downward radiative flux at interfaces [W m-2], shortwave plus longwave."""
        shortwave_down, shortwave_up = self.shortwave_fluxes(state)
        longwave_down, longwave_up = self.longwave_fluxes(state)
        return (shortwave_down - shortwave_up) + (longwave_down - longwave_up)

    def heating_rate(self, state: ModelState) -> NDArray[np.float64]:
        """Radiative heating rate at layer midpoints [K s-1].

        A layer warms when more net downward flux enters through its top than
        leaves through its bottom,

        .. math::
            \\frac{\\partial T}{\\partial t}
            = \\frac{g}{c_p} \\frac{F_{k+1} - F_k}{\\Delta p_k}
        """
        gravity = self.config.planet.gravity
        cp = self.config.atmosphere.gas.specific_heat_cp
        net = self.net_flux(state)
        delta_pressure = state.grid.atmosphere.delta_pressure
        return gravity * (net[1:] - net[:-1]) / (cp * delta_pressure)

    # =======================================================================
    # Contrato
    # =======================================================================

    def compute(self, state: ModelState, dt: float) -> Tendencies:
        """Compute the radiative heating tendency and publish surface diagnostics."""
        self.require_initialized()

        shortwave_down, _ = self.shortwave_fluxes(state)
        longwave_down, longwave_up = self.longwave_fluxes(state)
        heating = self.heating_rate(state)

        state.shortwave_down = float(shortwave_down[0])
        state.longwave_down = float(longwave_down[0])
        state.longwave_up = float(longwave_up[0])
        state.net_radiation = self.net_surface_flux(state)
        state.radiative_heating_rate = heating

        return Tendencies(atmospheric_temperature=heating)

    def diagnostics(self, state: ModelState) -> dict[str, float]:
        """Return the scalar diagnostics of the radiation scheme."""
        shortwave_down, shortwave_up = self.shortwave_fluxes(state)
        longwave_down, longwave_up = self.longwave_fluxes(state)
        heating = self.heating_rate(state)
        return {
            "cosine_solar_zenith": self.cosine_solar_zenith_angle(state),
            "insolation_top": self.top_of_atmosphere_insolation(state),
            "solar_constant_at_season": self.solar_constant_at_season(state),
            "shortwave_down_surface": float(shortwave_down[0]),
            "shortwave_up_surface": float(shortwave_up[0]),
            "longwave_down_surface": float(longwave_down[0]),
            "longwave_up_surface": float(longwave_up[0]),
            "surface_emission": self.surface_emission(state),
            "net_radiation": self.net_surface_flux(state),
            "net_longwave_surface": float(longwave_down[0] - longwave_up[0]),
            "outgoing_longwave_top": float(longwave_up[-1]),
            "max_heating_rate": float(np.max(heating)),
            "min_heating_rate": float(np.min(heating)),
        }


__all__ = ["GrayRadiation"]
