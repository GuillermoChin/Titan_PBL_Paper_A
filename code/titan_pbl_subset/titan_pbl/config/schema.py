"""
Configuration schema (Tier 2 validation layer).

Every tunable parameter of the framework is declared here as a Pydantic model
field. Nothing in the physics code is allowed to read a raw dictionary: all
configuration must pass through :class:`TitanPBLConfig` first, so that invalid
combinations are rejected before a single time step is integrated.

Design rules
------------
1. No physics module may hardcode a value that appears in this schema.
2. Every field carries its physical unit in the description.
3. Validators enforce physical plausibility, not just type correctness.
4. Adding a new scheme means adding a new nested model, never editing the core.

Reference
---------
Hartwick, Rafkin & Soto (2026), PSJ 7:178, DOI 10.3847/PSJ/ae809d
"""

from __future__ import annotations

from typing import Annotated, Any, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    field_validator,
    model_validator,
)

# ---------------------------------------------------------------------------
# Tipos anotados reutilizables
# ---------------------------------------------------------------------------

PositiveFloat = Annotated[float, Field(gt=0.0)]
NonNegativeFloat = Annotated[float, Field(ge=0.0)]
Fraction = Annotated[float, Field(ge=0.0, le=1.0)]
PositiveInt = Annotated[int, Field(gt=0)]


class _Base(BaseModel):
    """Modelo base: prohibe campos desconocidos para atrapar errores de tipeo."""

    model_config = ConfigDict(
        extra="forbid",
        validate_assignment=True,
        frozen=False,
    )


# ---------------------------------------------------------------------------
# Metadatos
# ---------------------------------------------------------------------------


class MetadataConfig(_Base):
    """Descriptive information about the configuration itself."""

    name: str = Field(description="Short identifier for this configuration")
    description: str = Field(default="", description="Human readable description")
    reference: str = Field(default="", description="Bibliographic reference")
    config_version: str = Field(default="0.1.0", description="Schema version")


# ---------------------------------------------------------------------------
# Planeta
# ---------------------------------------------------------------------------


class PlanetConfig(_Base):
    """Bulk planetary properties."""

    name: str = Field(default="Titan", description="Planet or moon name")
    gravity: PositiveFloat = Field(description="Surface gravity [m s-2]")
    radius: PositiveFloat = Field(description="Mean radius [m]")
    sol_length: PositiveFloat = Field(description="Solar day length [s]")
    minute_length: PositiveFloat = Field(
        description="Length of one planetary minute [s], used as base time unit"
    )
    orbital_period: PositiveFloat = Field(description="Orbital period [s]")
    obliquity: float = Field(
        default=26.73, description="Obliquity to the Sun [deg]"
    )
    orbital_eccentricity: Fraction = Field(
        default=0.0565, description="Orbital eccentricity [dimensionless]"
    )
    perihelion_solar_longitude: Annotated[float, Field(ge=0.0, lt=360.0)] = Field(
        default=280.0, description="Solar longitude at perihelion Ls0 [deg]"
    )

    @model_validator(mode="after")
    def _check_time_consistency(self) -> "PlanetConfig":
        # El "minuto planetario" debe dividir razonablemente al sol.
        ratio = self.sol_length / self.minute_length
        if not (100.0 < ratio < 100000.0):
            raise ValueError(
                f"sol_length / minute_length = {ratio:.1f} is implausible; "
                "expected roughly 1440 planetary minutes per sol."
            )
        if self.orbital_period <= self.sol_length:
            raise ValueError("orbital_period must exceed sol_length.")
        return self


# ---------------------------------------------------------------------------
# Sitio
# ---------------------------------------------------------------------------


class SiteConfig(_Base):
    """Geographic location and season of the simulated column."""

    latitude: Annotated[float, Field(ge=-90.0, le=90.0)] = Field(
        description="Latitude [deg], positive north"
    )
    longitude: Annotated[float, Field(ge=-360.0, le=360.0)] = Field(
        description="Longitude [deg], positive east"
    )
    solar_longitude_ls: Annotated[float, Field(ge=0.0, lt=360.0)] = Field(
        description="Solar longitude Ls [deg], season indicator"
    )
    start_local_time: Annotated[float, Field(ge=0.0, lt=24.0)] = Field(
        description="Local solar time at simulation start [h]"
    )
    surface_elevation: float = Field(
        default=0.0, description="Surface elevation above datum [m]"
    )


# ---------------------------------------------------------------------------
# Atmosfera
# ---------------------------------------------------------------------------


class GasConfig(_Base):
    """Thermodynamic properties of the background atmospheric gas."""

    specific_gas_constant: PositiveFloat = Field(
        description="Specific gas constant R [J kg-1 K-1]"
    )
    specific_heat_cp: PositiveFloat = Field(
        description="Specific heat at constant pressure [J kg-1 K-1]"
    )
    reference_density: PositiveFloat = Field(
        description="Reference near surface density [kg m-3]"
    )
    thermal_conductivity: PositiveFloat = Field(
        description="Molecular thermal conductivity [W m-1 K-1]"
    )

    @property
    def kappa(self) -> float:
        """Poisson constant R/cp [dimensionless], used for potential temperature."""
        return self.specific_gas_constant / self.specific_heat_cp


class VerticalGridConfig(_Base):
    """Definition of the atmospheric vertical grid."""

    n_layers: PositiveInt = Field(description="Number of atmospheric layers")
    surface_pressure: PositiveFloat = Field(description="Surface pressure [Pa]")
    top_pressure: PositiveFloat = Field(description="Model top pressure [Pa]")
    top_altitude: PositiveFloat = Field(description="Model top altitude [m]")
    spacing_type: Literal["sigma", "log_pressure", "geometric"] = Field(
        default="sigma", description="Vertical coordinate spacing strategy"
    )
    spacing_exponent: PositiveFloat = Field(
        default=1.0,
        description="Stretching exponent; >1 refines layers near the surface",
    )

    @model_validator(mode="after")
    def _check_pressures(self) -> "VerticalGridConfig":
        if self.top_pressure >= self.surface_pressure:
            raise ValueError(
                "top_pressure must be lower than surface_pressure "
                f"(got {self.top_pressure} >= {self.surface_pressure})."
            )
        return self


class ReferenceStateConfig(_Base):
    """Reference state used to define potential temperature."""

    reference_pressure: PositiveFloat = Field(
        description="Reference pressure P0 for potential temperature [Pa]"
    )
    reference_temperature: PositiveFloat = Field(
        description="Reference temperature [K]"
    )


class AtmosphereInitConfig(_Base):
    """Initial atmospheric state."""

    profile_source: Literal["hasi", "isothermal", "custom"] = Field(
        default="hasi", description="Source of the initial thermal profile"
    )
    hasi_file: str | None = Field(
        default=None, description="Path to the HASI descent profile file"
    )
    custom_file: str | None = Field(
        default=None, description="Path to a user supplied profile"
    )
    isothermal_temperature: PositiveFloat | None = Field(
        default=None, description="Uniform temperature [K] if profile_source=isothermal"
    )
    zonal_wind: float = Field(default=0.0, description="Initial zonal wind [m s-1]")
    meridional_wind: float = Field(
        default=0.0, description="Initial meridional wind [m s-1]"
    )

    @model_validator(mode="after")
    def _check_source_requirements(self) -> "AtmosphereInitConfig":
        if self.profile_source == "hasi" and not self.hasi_file:
            raise ValueError("profile_source='hasi' requires hasi_file to be set.")
        if self.profile_source == "custom" and not self.custom_file:
            raise ValueError("profile_source='custom' requires custom_file to be set.")
        if self.profile_source == "isothermal" and self.isothermal_temperature is None:
            raise ValueError(
                "profile_source='isothermal' requires isothermal_temperature."
            )
        return self


class AtmosphereConfig(_Base):
    """Full atmospheric configuration block."""

    gas: GasConfig
    vertical_grid: VerticalGridConfig
    reference_state: ReferenceStateConfig
    initialization: AtmosphereInitConfig


# ---------------------------------------------------------------------------
# Subsuelo
# ---------------------------------------------------------------------------


class SubsurfaceConfig(_Base):
    """Subsurface layer structure and lower boundary condition."""

    scheme: str = Field(
        default="conduction", description="Registered subsurface scheme name"
    )
    n_layers: PositiveInt = Field(description="Number of subsurface layers")
    layer_thickness_mode: Literal["geometric", "explicit"] = Field(
        default="geometric", description="How layer thicknesses are generated"
    )
    first_layer_thickness: PositiveFloat = Field(
        description="Thickness of the uppermost subsurface layer [m]"
    )
    thickness_growth_factor: PositiveFloat = Field(
        default=2.0, description="Geometric ratio between consecutive layers"
    )
    explicit_thicknesses: list[PositiveFloat] | None = Field(
        default=None, description="Explicit layer thicknesses [m]"
    )
    deep_temperature: PositiveFloat = Field(
        description="Fixed temperature at the base of the domain, TMN [K]"
    )
    initial_profile_mode: Literal["linear", "isothermal", "custom"] = Field(
        default="linear", description="Initial subsurface temperature profile"
    )
    surface_coupling: Literal["skin_is_top_layer", "separate_skin"] = Field(
        default="skin_is_top_layer",
        description=(
            "How the skin temperature relates to the regolith column. With "
            "'skin_is_top_layer' the skin is the uppermost subsurface layer, "
            "so the surface energy budget and the conduction to the next layer "
            "share one heat capacity, which is what Figures 3(c) and 3(e) of "
            "the reference paper describe. With 'separate_skin' the skin is an "
            "independent film above the column"
        ),
    )
    interface_averaging: Literal["harmonic", "arithmetic", "geometric"] = Field(
        default="harmonic",
        description=(
            "How conductivity is averaged across a layer interface. Harmonic "
            "is the exact result for conduction through layers in series and "
            "is the default; the reference paper does not state its choice"
        ),
    )
    geothermal_flux: float = Field(
        default=0.0, description="Prescribed geothermal heat flux [W m-2]"
    )

    @model_validator(mode="after")
    def _check_thickness_mode(self) -> "SubsurfaceConfig":
        if self.layer_thickness_mode == "explicit":
            if self.explicit_thicknesses is None:
                raise ValueError(
                    "layer_thickness_mode='explicit' requires explicit_thicknesses."
                )
            if len(self.explicit_thicknesses) != self.n_layers:
                raise ValueError(
                    f"explicit_thicknesses has {len(self.explicit_thicknesses)} "
                    f"entries but n_layers is {self.n_layers}."
                )
        if self.layer_thickness_mode == "geometric":
            if self.thickness_growth_factor < 1.0:
                raise ValueError(
                    "thickness_growth_factor must be >= 1.0 for geometric spacing."
                )
        return self


# ---------------------------------------------------------------------------
# Regolito
# ---------------------------------------------------------------------------


class MaterialConfig(_Base):
    """Thermophysical properties of a single regolith constituent."""

    specific_heat: PositiveFloat = Field(
        description="Specific heat capacity [J kg-1 K-1]"
    )
    density: PositiveFloat = Field(description="Density [kg m-3]")
    thermal_conductivity: PositiveFloat = Field(
        description="Thermal conductivity [W m-1 K-1]"
    )


class RegolithComponentsConfig(_Base):
    """The three constituents mixed by Equations 6-9 of the reference paper."""

    solid_particle: MaterialConfig
    liquid_methane: MaterialConfig
    pore_air: MaterialConfig


class RegolithConfig(_Base):
    """Bulk regolith description.

    Effective thermal properties are never prescribed directly: they are
    derived from the constituent properties and the porosity/saturation pair
    through Equations 6-9 of Hartwick et al. (2026).
    """

    porosity: Fraction = Field(description="Bulk porosity phi [dimensionless]")
    methane_saturation: Fraction = Field(
        description="Pore space methane saturation beta [dimensionless]"
    )
    density_override: PositiveFloat | None = Field(
        default=None,
        description="Optional bulk density override [kg m-3]; null derives it",
    )
    conductivity_override: PositiveFloat | None = Field(
        default=None,
        description=(
            "Optional bulk thermal conductivity override [W m-1 K-1]; null "
            "derives it from Equation 8. Required to reproduce sensitivity "
            "tests 2-5 of Table 3, whose conductivities cannot be produced by "
            "any physical porosity"
        ),
    )
    volumetric_heat_capacity_override: PositiveFloat | None = Field(
        default=None,
        description=(
            "Optional volumetric heat capacity override [J m-3 K-1]; null "
            "derives it from Equation 9"
        ),
    )
    components: RegolithComponentsConfig


# ---------------------------------------------------------------------------
# Radiacion
# ---------------------------------------------------------------------------


class RadiationConfig(_Base):
    """Radiative transfer scheme parameters."""

    scheme: str = Field(default="gray", description="Registered radiation scheme name")
    solar_scheme: Literal["two_stream", "beer_lambert", "none"] = Field(
        default="two_stream",
        description=(
            "Shortwave solver. 'two_stream' is the analytic multiple "
            "scattering solution of Chatain et al. (2022) Appendix A.2; "
            "'beer_lambert' is pure extinction with no scattering"
        ),
    )
    ir_scheme: Literal["schneider", "none"] = Field(
        default="schneider",
        description="Longwave solver; 'schneider' follows Equations A13 to A16",
    )
    apply_eccentricity_correction: bool = Field(
        default=True,
        description="Apply Equation A1, the orbital eccentricity correction to S0",
    )
    mean_cosine: PositiveFloat = Field(
        default=0.5773502691896258,
        description=(
            "Mean cosine mu-bar of the two-stream approximation "
            "[dimensionless]; 1/sqrt(3) for a collimated beam source"
        ),
    )
    solar_constant: NonNegativeFloat = Field(
        description="Insolation at the top of the model domain S0 [W m-2]"
    )
    ir_top_flux: NonNegativeFloat = Field(
        description="Downward IR flux prescribed at model top I_IR [W m-2]"
    )
    surface_albedo: Fraction = Field(description="Surface albedo alpha [dimensionless]")
    surface_emissivity: Fraction = Field(
        description="Surface emissivity epsilon [dimensionless]"
    )
    optical_depth_visible: NonNegativeFloat = Field(
        description="Visible optical depth tau_vis [dimensionless]"
    )
    optical_depth_ir: NonNegativeFloat = Field(
        description="Infrared optical depth tau_IR [dimensionless]"
    )
    visible_scaling_exponent: float = Field(
        description="Visible optical depth scaling exponent gamma [dimensionless]"
    )
    ir_weighting_coefficient: float = Field(
        description="Infrared weighting coefficient f [dimensionless]"
    )
    single_scattering_albedo: Fraction = Field(
        description="Atmospheric single scattering albedo omega0 [dimensionless]"
    )
    asymmetry_factor: Annotated[float, Field(ge=-1.0, le=1.0)] = Field(
        description="Scattering asymmetry factor g [dimensionless]"
    )
    update_interval_steps: PositiveInt = Field(
        default=1, description="Radiation update frequency in time steps"
    )


# ---------------------------------------------------------------------------
# Capa superficial
# ---------------------------------------------------------------------------


class SurfaceLayerConfig(_Base):
    """Surface layer turbulent exchange parameters (MOST)."""

    scheme: str = Field(default="most", description="Registered surface scheme name")
    roughness_length: PositiveFloat = Field(
        description="Aerodynamic roughness length z0 [m]"
    )
    initial_friction_velocity: PositiveFloat = Field(
        description="Friction velocity u* imposed at the first step [m s-1]"
    )
    beljaars_gustiness: NonNegativeFloat = Field(
        description="Beljaars gustiness coefficient beta [dimensionless]"
    )
    von_karman_constant: PositiveFloat = Field(
        default=0.4, description="von Karman constant [dimensionless]"
    )
    min_friction_velocity: PositiveFloat = Field(
        description="Numerical floor for u* [m s-1]"
    )
    roughness_length_heat: PositiveFloat | None = Field(
        default=None,
        description=(
            "Thermal roughness length z0h [m]; null makes it equal to the "
            "aerodynamic roughness length"
        ),
    )
    stability_functions: Literal["businger_dyer", "beljaars_holtslag", "neutral"] = (
        Field(
            default="businger_dyer",
            description=(
                "Monin-Obukhov stability functions. The reference paper cites "
                "Janjic (1994) but does not state the functional forms"
            ),
        )
    )
    unstable_coefficient: PositiveFloat = Field(
        default=16.0, description="Coefficient of the unstable similarity functions"
    )
    stable_coefficient: PositiveFloat = Field(
        default=5.0, description="Coefficient of the stable similarity functions"
    )
    apply_beljaars_correction: bool = Field(
        default=True,
        description="Add the Beljaars (1995) convective velocity scale under free convection",
    )
    minimum_convective_depth: PositiveFloat = Field(
        default=100.0,
        description=(
            "Depth used for the convective velocity scale before the boundary "
            "layer has been diagnosed [m]"
        ),
    )
    max_iterations: PositiveInt = Field(
        default=30, description="Iteration cap for the implicit surface layer solve"
    )
    convergence_tolerance: PositiveFloat = Field(
        default=1.0e-8, description="Relative convergence tolerance for u*"
    )
    max_stability_parameter: PositiveFloat = Field(
        default=10.0, description="Cap on |z/L| to keep the similarity functions bounded"
    )
    implicit_skin_update: bool = Field(
        default=True,
        description=(
            "Linearise the surface energy budget about the current skin "
            "temperature before advancing it. Required because the exchange "
            "coefficient can make the surface to air coupling relax faster "
            "than one time step, which an explicit update turns into a "
            "relaxation oscillation"
        ),
    )
    min_denominator_fraction: Fraction = Field(
        default=0.25,
        description=(
            "Floor on the similarity denominator as a fraction of its neutral "
            "value. Required because the reference configuration places the "
            "lowest level only about 13 roughness lengths above the surface, "
            "where the unstable stability correction can otherwise cancel the "
            "logarithmic term and make the exchange coefficient diverge"
        ),
    )


# ---------------------------------------------------------------------------
# Turbulencia
# ---------------------------------------------------------------------------


class TurbulenceConfig(_Base):
    """Boundary layer turbulence scheme parameters."""

    scheme: str = Field(default="myj", description="Registered turbulence scheme name")
    tke_critical: PositiveFloat = Field(
        description="TKE threshold used to diagnose PBL top [m2 s-2]"
    )
    tke_minimum: PositiveFloat = Field(description="Numerical floor for TKE [m2 s-2]")
    mixing_length_max: PositiveFloat = Field(
        description="Asymptotic maximum mixing length [m]"
    )
    richardson_critical: PositiveFloat = Field(
        description="Critical Richardson number [dimensionless]"
    )
    diffusivity_min: PositiveFloat = Field(
        description="Minimum eddy diffusivity [m2 s-1]"
    )
    diffusivity_max: PositiveFloat = Field(
        description="Maximum eddy diffusivity [m2 s-1]"
    )
    length_scale: Literal["blackadar", "constant"] = Field(
        default="blackadar", description="Master mixing length formulation"
    )
    apply_stable_length_limit: bool = Field(
        default=True,
        description="Limit the mixing length by the buoyancy scale in stable layers",
    )
    stable_length_coefficient: PositiveFloat = Field(
        default=0.53,
        description="Coefficient of the stable mixing length limit l <= c q / N",
    )
    dissipation_constant_b1: PositiveFloat = Field(
        default=16.6,
        description="Mellor-Yamada dissipation constant B1 [dimensionless]",
    )
    neutral_momentum_stability: PositiveFloat = Field(
        default=0.3933,
        description="Neutral value of the momentum stability function S_M",
    )
    neutral_heat_stability: PositiveFloat = Field(
        default=0.4939,
        description="Neutral value of the heat stability function S_H",
    )
    galperin_heat_denominator: PositiveFloat = Field(
        default=34.6764, description="Galperin et al. (1988) coefficient for S_H"
    )
    galperin_momentum_numerator: PositiveFloat = Field(
        default=3.0858, description="Galperin et al. (1988) numerator coefficient for S_M"
    )
    galperin_momentum_denominator: PositiveFloat = Field(
        default=6.1272,
        description="Galperin et al. (1988) second denominator coefficient for S_M",
    )
    stability_limit_stable: float = Field(
        default=0.0233,
        description="Upper bound on the Galperin stability parameter G_H",
    )
    stability_limit_unstable: float = Field(
        default=-0.28,
        description="Lower bound on the Galperin stability parameter G_H",
    )
    tke_maximum: PositiveFloat = Field(
        default=10.0, description="Numerical ceiling for the turbulent kinetic energy [m2 s-2]"
    )
    implicit_dissipation: bool = Field(
        default=True,
        description=(
            "Treat the dissipation term implicitly. The buoyancy production "
            "and the diffusivity form a positive feedback whose explicit "
            "dissipation sink overshoots at the reference time step"
        ),
    )
    include_tke_transport: bool = Field(
        default=True, description="Include vertical turbulent transport of TKE"
    )
    tke_transport_prandtl: PositiveFloat = Field(
        default=1.0, description="Ratio of the TKE diffusivity to the momentum one"
    )
    prescribed_shear: NonNegativeFloat = Field(
        default=0.0,
        description=(
            "Background vertical wind shear [s-1] used for the shear production "
            "term. Zero reproduces the reference configuration, in which the "
            "horizontal winds are invariant and buoyancy is the only source"
        ),
    )

    @model_validator(mode="after")
    def _check_diffusivity_range(self) -> "TurbulenceConfig":
        if self.diffusivity_min >= self.diffusivity_max:
            raise ValueError("diffusivity_min must be smaller than diffusivity_max.")
        if self.tke_minimum >= self.tke_critical:
            raise ValueError("tke_minimum must be smaller than tke_critical.")
        return self


# ---------------------------------------------------------------------------
# Metano
# ---------------------------------------------------------------------------


class MethaneSaturationConfig(_Base):
    """Saturation vapour pressure formulation for methane."""

    formulation: Literal["antoine", "clausius_clapeyron", "tabulated"] = Field(
        default="antoine", description="Saturation curve formulation"
    )
    antoine_A: float = Field(default=3.9895, description="Antoine coefficient A")
    antoine_B: float = Field(default=443.028, description="Antoine coefficient B")
    antoine_C: float = Field(default=-0.49, description="Antoine coefficient C")
    tabulated_file: str | None = Field(
        default=None, description="Path to a tabulated saturation curve"
    )

    @model_validator(mode="after")
    def _check_tabulated(self) -> "MethaneSaturationConfig":
        if self.formulation == "tabulated" and not self.tabulated_file:
            raise ValueError("formulation='tabulated' requires tabulated_file.")
        return self


class MethaneInitProfileConfig(_Base):
    """Initial atmospheric methane distribution."""

    source: Literal["gautier2024", "niemann2010", "zero", "uniform", "custom"] = Field(
        default="gautier2024", description="Reference methane profile"
    )
    relative_humidity: Fraction | None = Field(
        default=None,
        description=(
            "Optional methane relative humidity used with the 'uniform' "
            "source; when set it overrides surface_mixing_ratio and the "
            "profile is built from the saturation curve"
        ),
    )
    custom_file: str | None = Field(
        default=None, description="Path to a user supplied methane profile"
    )
    surface_mixing_ratio: NonNegativeFloat = Field(
        default=0.0, description="Near surface methane mass mixing ratio [kg kg-1]"
    )

    @model_validator(mode="after")
    def _check_custom(self) -> "MethaneInitProfileConfig":
        if self.source == "custom" and not self.custom_file:
            raise ValueError("source='custom' requires custom_file.")
        return self


class SPSConfig(_Base):
    """Saturated Pore Space scheme parameters (Eqs. 12-13)."""

    porosity: Fraction | None = Field(
        default=None, description="Scheme porosity; null inherits regolith.porosity"
    )


class NSMTConfig(_Base):
    """Near Surface Methane Table scheme parameters (Eqs. 14-17)."""

    methane_table_depth: PositiveFloat = Field(
        description="Depth to the methane table z [m]"
    )
    porosity: Fraction | None = Field(
        default=None, description="Scheme porosity; null inherits regolith.porosity"
    )


class MRPSConfig(_Base):
    """Methane Rich Pore Space scheme parameters (Eqs. 18-20)."""

    pore_methane_abundance: Fraction = Field(
        description="Regolith moisture parameter beta [dimensionless]"
    )
    porosity: Fraction | None = Field(
        default=None, description="Scheme porosity; null inherits regolith.porosity"
    )


class MethaneConfig(_Base):
    """Methane evaporation configuration block."""

    scheme: Literal["none", "sps", "nsmt", "mrps"] = Field(
        default="none", description="Active methane evaporation scheme"
    )
    latent_heat_vaporization: PositiveFloat = Field(
        description="Latent heat of vaporization XLV [J kg-1]"
    )
    molecular_diffusivity_air: PositiveFloat = Field(
        description="Methane diffusivity through air D_a [m2 s-1]"
    )
    tortuosity_exponent: PositiveFloat = Field(
        default=4.0 / 3.0, description="Millington-Quirk exponent [dimensionless]"
    )
    molar_mass: PositiveFloat = Field(
        default=0.016043, description="Molar mass of methane [kg mol-1]"
    )
    allow_atmospheric_condensation: bool = Field(
        default=False,
        description="Enable atmospheric condensation; False replicates the paper",
    )
    infinite_reservoir: bool = Field(
        default=True,
        description="Assume an infinite subsurface reservoir; True replicates the paper",
    )
    initial_atmospheric_profile: MethaneInitProfileConfig
    saturation: MethaneSaturationConfig
    sps: SPSConfig = Field(default_factory=SPSConfig)
    nsmt: NSMTConfig | None = None
    mrps: MRPSConfig | None = None

    @model_validator(mode="after")
    def _check_active_scheme_present(self) -> "MethaneConfig":
        if self.scheme == "nsmt" and self.nsmt is None:
            raise ValueError("scheme='nsmt' requires the methane.nsmt block.")
        if self.scheme == "mrps" and self.mrps is None:
            raise ValueError("scheme='mrps' requires the methane.mrps block.")
        return self


# ---------------------------------------------------------------------------
# Simulacion
# ---------------------------------------------------------------------------


class SimulationConfig(_Base):
    """Temporal integration settings."""

    n_sols: PositiveFloat = Field(description="Total simulated duration [sols]")
    spinup_sols: NonNegativeFloat = Field(
        description="Leading sols excluded from the analysis [sols]"
    )
    timestep: PositiveFloat = Field(description="Integration time step [s]")
    time_integrator: Literal["rk4", "euler", "crank_nicolson"] = Field(
        default="rk4", description="Time integration scheme"
    )
    diffusion_solver: Literal["crank_nicolson", "explicit", "implicit"] = Field(
        default="crank_nicolson", description="Vertical diffusion solver"
    )
    max_iterations_implicit: PositiveInt = Field(
        default=50, description="Iteration cap for implicit solvers"
    )
    convergence_tolerance: PositiveFloat = Field(
        default=1.0e-10, description="Convergence tolerance for implicit solvers"
    )

    @model_validator(mode="after")
    def _check_spinup(self) -> "SimulationConfig":
        if self.spinup_sols >= self.n_sols:
            raise ValueError(
                f"spinup_sols ({self.spinup_sols}) must be smaller than "
                f"n_sols ({self.n_sols}); otherwise no output remains for analysis."
            )
        return self


# ---------------------------------------------------------------------------
# Salida
# ---------------------------------------------------------------------------


class OutputConfig(_Base):
    """Output writing settings."""

    directory: str = Field(default="data/outputs", description="Output directory")
    format: Literal["netcdf", "hdf5", "none"] = Field(
        default="netcdf", description="Output file format"
    )
    write_interval_steps: PositiveInt = Field(
        default=1, description="Write frequency in time steps"
    )
    variables: list[str] = Field(
        default_factory=list, description="Variables to store"
    )


# ---------------------------------------------------------------------------
# Configuracion raiz
# ---------------------------------------------------------------------------


class TitanPBLConfig(_Base):
    """Root configuration object for a titan_pbl simulation.

    This is the single source of truth for every tunable parameter. Physics
    modules receive an instance of this class and must read all their
    constants from it.

    Examples
    --------
    >>> from titan_pbl.config import load_config
    >>> cfg = load_config("huygens_baseline")
    >>> cfg.subsurface.n_layers
    10
    >>> cfg = load_config("huygens_baseline",
    ...                   overrides={"regolith.porosity": 0.02})
    >>> cfg.regolith.porosity
    0.02
    """

    metadata: MetadataConfig
    planet: PlanetConfig
    site: SiteConfig
    atmosphere: AtmosphereConfig
    subsurface: SubsurfaceConfig
    regolith: RegolithConfig
    radiation: RadiationConfig
    surface_layer: SurfaceLayerConfig
    turbulence: TurbulenceConfig
    methane: MethaneConfig
    simulation: SimulationConfig
    output: OutputConfig = Field(default_factory=OutputConfig)

    # -----------------------------------------------------------------------
    # Validaciones cruzadas entre bloques
    # -----------------------------------------------------------------------

    @model_validator(mode="after")
    def _check_cross_block_consistency(self) -> "TitanPBLConfig":
        # El paso de tiempo no puede exceder la duracion total simulada.
        total_seconds = self.simulation.n_sols * self.planet.sol_length
        if self.simulation.timestep >= total_seconds:
            raise ValueError(
                "simulation.timestep exceeds the total simulated duration."
            )

        # Un regolito con saturacion de metano no nula exige un esquema activo.
        if self.regolith.methane_saturation > 0.0 and self.methane.scheme == "none":
            raise ValueError(
                "regolith.methane_saturation > 0 but methane.scheme is 'none'; "
                "either enable a methane scheme or set the saturation to zero."
            )

        # NSMT: la tabla de metano debe caer dentro del dominio del subsuelo.
        if self.methane.scheme == "nsmt" and self.methane.nsmt is not None:
            total_depth = self.subsurface_total_depth
            if self.methane.nsmt.methane_table_depth > total_depth:
                raise ValueError(
                    f"methane.nsmt.methane_table_depth "
                    f"({self.methane.nsmt.methane_table_depth} m) lies below the "
                    f"subsurface domain ({total_depth:.2f} m)."
                )
        return self

    # -----------------------------------------------------------------------
    # Propiedades derivadas
    # -----------------------------------------------------------------------

    @property
    def subsurface_layer_thicknesses(self) -> list[float]:
        """Thickness of each subsurface layer [m], derived from the config."""
        sub = self.subsurface
        if sub.layer_thickness_mode == "explicit":
            assert sub.explicit_thicknesses is not None
            return list(sub.explicit_thicknesses)
        # Progresion geometrica: dz_k = dz_0 * r^k
        return [
            sub.first_layer_thickness * sub.thickness_growth_factor**k
            for k in range(sub.n_layers)
        ]

    @property
    def subsurface_total_depth(self) -> float:
        """Total depth of the subsurface domain [m]."""
        return float(sum(self.subsurface_layer_thicknesses))

    @property
    def effective_porosity(self) -> float:
        """Porosity actually used by the active methane scheme [dimensionless].

        Each scheme may override the bulk regolith porosity; when its own
        value is null it inherits ``regolith.porosity``.
        """
        scheme = self.methane.scheme
        override: float | None = None
        if scheme == "sps":
            override = self.methane.sps.porosity
        elif scheme == "nsmt" and self.methane.nsmt is not None:
            override = self.methane.nsmt.porosity
        elif scheme == "mrps" and self.methane.mrps is not None:
            override = self.methane.mrps.porosity
        return self.regolith.porosity if override is None else override

    @property
    def total_steps(self) -> int:
        """Total number of time steps in the simulation."""
        total_seconds = self.simulation.n_sols * self.planet.sol_length
        return int(round(total_seconds / self.simulation.timestep))

    @property
    def spinup_steps(self) -> int:
        """Number of leading time steps discarded from the analysis."""
        spinup_seconds = self.simulation.spinup_sols * self.planet.sol_length
        return int(round(spinup_seconds / self.simulation.timestep))

    def summary(self) -> dict[str, Any]:
        """Return a compact dictionary describing the configuration."""
        return {
            "name": self.metadata.name,
            "atmospheric_layers": self.atmosphere.vertical_grid.n_layers,
            "subsurface_layers": self.subsurface.n_layers,
            "subsurface_depth_m": round(self.subsurface_total_depth, 3),
            "methane_scheme": self.methane.scheme,
            "effective_porosity": self.effective_porosity,
            "n_sols": self.simulation.n_sols,
            "timestep_s": self.simulation.timestep,
            "total_steps": self.total_steps,
        }


__all__ = [
    "TitanPBLConfig",
    "MetadataConfig",
    "PlanetConfig",
    "SiteConfig",
    "AtmosphereConfig",
    "GasConfig",
    "VerticalGridConfig",
    "ReferenceStateConfig",
    "AtmosphereInitConfig",
    "SubsurfaceConfig",
    "RegolithConfig",
    "RegolithComponentsConfig",
    "MaterialConfig",
    "RadiationConfig",
    "SurfaceLayerConfig",
    "TurbulenceConfig",
    "MethaneConfig",
    "MethaneSaturationConfig",
    "MethaneInitProfileConfig",
    "SPSConfig",
    "NSMTConfig",
    "MRPSConfig",
    "SimulationConfig",
    "OutputConfig",
]
