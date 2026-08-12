# Pies de figura y tablas — Paper A

Texto listo para pegar. Los pies van en inglés, como el manuscrito.

---

## Figure captions

**Figure 1.** Radiative energy balance of a thin gray layer. The layer absorbs the upward flux from the surface and the downward flux from above, one from each direction, but emits into both hemispheres. The doubling of the emission term against the single absorption terms is what sets the heating threshold of Equation (4).

*Archivo:* `figures/fig01_schematic.png` — *sección:* §2.1

---

**Figure 2.** Radiative heating rate against altitude for surface-to-air temperature contrasts of 0 to 4 K. The profile is negative everywhere at small contrast; positive values appear near the ground only at large contrast, and the vertical gradient changes sign at a larger contrast still.

*Archivo:* `figures/fig02_heating_profiles.png` — *sección:* §3.2

---

**Figure 3.** Destabilisation rate, defined as the difference in radiative heating rate between the surface layer and 4 km, against surface-to-air contrast. The zero crossing at 2.53 K marks the destabilisation threshold. The shaded band spans the contrasts recovered from the archived output of Hartwick et al. (2026).

*Archivo:* `figures/fig03_destabilisation.png` — *sección:* §3.3

---

**Figure 4.** Time required to erode the observed static stability over the lowest four kilometres, against surface-to-air contrast, on logarithmic axes. The horizontal band spans typical one-dimensional integration lengths; the vertical band spans contrasts reported in published Titan studies. The two do not intersect.

*Archivo:* `figures/fig04_mixing_time.png` — *sección:* §3.5

---

**Figure 5.** Radiative heating rate against altitude, in K per Titan sol, for (a) the baseline run and (b) the suppressed-flux run of Hartwick et al. (2026). Solid line: archived longwave heating rate. Dash-dotted line: archived total, longwave plus shortwave. Dashed line: the profile computed here. Agreement in the lowest layer is 0.90 %.

*Archivo:* `figures/fig05_validation.png` — *sección:* §3.6

---

**Figure 6.** Time series over Titan sols 2 to 5 from the archived output of Hartwick et al. (2026), for the baseline run (solid) and the suppressed-flux run (dashed). (a) Surface-to-air temperature contrast in K, with the numerical heating threshold of 2.15 K marked. The contrast stays below it throughout, peaking at 0.67 K. (b) Boundary layer depth in m. The suppressed-flux layer grows by 55 % across the four sols and does not level off.

*Archivo:* `figures/fig06_timeseries.png` — *sección:* §4.4

---

## Table captions and contents

**Table 1.** Model configuration. Parameter, symbol, value, units and source for every quantity entering the radiative transfer calculation.

*Archivo:* `tables/table01_configuration.csv` — *sección:* §3.1

| parameter                | symbol   |        value | units      | source          |
|:-------------------------|:---------|-------------:|:-----------|:----------------|
| Solar constant           | S0       |    15.04     | W m-2      | Hartwick+ 2026  |
| Downward IR at model top | I_IR     |     2        | W m-2      | Hartwick+ 2026  |
| Surface albedo           | alpha    |     0.3      | -          | Hartwick+ 2026  |
| Surface emissivity       | epsilon  |     0.95     | -          | Hartwick+ 2026  |
| Visible optical depth    | tau_vis  |     8        | -          | Chatain+ 2022   |
| Infrared optical depth   | tau_IR   |    10        | -          | Hartwick+ 2026  |
| Visible scaling exponent | gamma    |     0.21     | -          | Schneider+ 2012 |
| IR weighting coefficient | f        |     0.15     | -          | Chatain+ 2022   |
| Single scattering albedo | omega0   |     0.9      | -          | Chatain+ 2022   |
| Asymmetry factor         | g        |     0.85     | -          | Chatain+ 2022   |
| Atmospheric layers       | N        |    56        | -          | Hartwick+ 2026  |
| Surface pressure         | p_s      |     1.45e+05 | Pa         | Hartwick+ 2026  |
| Model top pressure       | p_t      | 66000        | Pa         | Hartwick+ 2026  |
| Specific heat            | c_p      |  1101        | J kg-1 K-1 | Hartwick+ 2026  |
| Surface gravity          | g        |     1.352    | m s-2      | standard        |

---

**Table 2.** Sweep over the three uncertain longwave parameters: the downward infrared flux at model top, the infrared optical depth and its weighting coefficient. No combination produces net heating in the lower atmosphere or reverses the sign of the vertical heating gradient.

*Archivo:* `tables/table02_parameter_sweep.csv` — *sección:* §3.4

|   I_IR_W_m2 |   tau_IR |    f |   max_heating_K_per_sol |   min_heating_K_per_sol |   surface_heating_K_per_sol |   gradient_K_per_sol |
|------------:|---------:|-----:|------------------------:|------------------------:|----------------------------:|---------------------:|
|           0 |       10 | 0.15 |               0.0045334 |               -0.083657 |                   -0.083657 |            -0.08819  |
|           0 |       10 | 0.5  |               0.0011225 |               -0.074562 |                   -0.074562 |            -0.075685 |
|           0 |       30 | 0.15 |               0.02111   |               -0.17017  |                   -0.17017  |            -0.19128  |
|           0 |       30 | 0.5  |               0.020783  |               -0.1453   |                   -0.1453   |            -0.16608  |
|           2 |       10 | 0.15 |               0.0059396 |               -0.08343  |                   -0.08343  |            -0.08937  |
|           2 |       10 | 0.5  |               0.0031088 |               -0.074125 |                   -0.074125 |            -0.077233 |
|           2 |       30 | 0.15 |               0.02111   |               -0.17017  |                   -0.17017  |            -0.19128  |
|           2 |       30 | 0.5  |               0.020783  |               -0.1453   |                   -0.1453   |            -0.16608  |
|           4 |       10 | 0.15 |               0.0073458 |               -0.083203 |                   -0.083203 |            -0.090549 |
|           4 |       10 | 0.5  |               0.005095  |               -0.073687 |                   -0.073687 |            -0.078782 |
|           4 |       30 | 0.15 |               0.02111   |               -0.17017  |                   -0.17017  |            -0.19128  |
|           4 |       30 | 0.5  |               0.020783  |               -0.1453   |                   -0.1453   |            -0.16608  |
|           8 |       10 | 0.15 |               0.010158  |               -0.08275  |                   -0.08275  |            -0.092908 |
|           8 |       10 | 0.5  |               0.0090675 |               -0.072812 |                   -0.072812 |            -0.081879 |
|           8 |       30 | 0.15 |               0.02111   |               -0.17017  |                   -0.17017  |            -0.19128  |
|           8 |       30 | 0.5  |               0.020784  |               -0.1453   |                   -0.1453   |            -0.16609  |

---

**Table 3.** Archived against computed longwave heating rates at four altitudes, with the relative difference at each. Agreement is best in the lowest layers, which is where the argument of Sections 2 and 3 operates.

*Archivo:* `tables/table03_validation.csv` — *sección:* §3.6

|   altitude_m |   archived_K_per_sol |   this_work_K_per_sol |   difference_percent |
|-------------:|---------------------:|----------------------:|---------------------:|
|       5.3879 |            -0.084191 |             -0.08343  |              0.90399 |
|     172.13   |            -0.072987 |             -0.073547 |              0.76723 |
|     556.4    |            -0.053481 |             -0.048319 |              9.6517  |
|    1415.8    |            -0.025836 |             -0.013686 |             47.028   |

---

**Table 4.** Thermal contrasts, sensible heat fluxes and boundary layer depths recovered from the archive, by Titan sol, for both runs.

*Archivo:* `tables/table04_contrasts.csv` — *sección:* §4.4

|   sol |   contrast_baseline_max_K |   contrast_suppressed_max_K |   delta_Ts_between_runs_max_K |   delta_Ts_between_runs_mean_K |   HFX_baseline_peak_W_m2 |   HFX_suppressed_peak_W_m2 |   HFX_ratio |   PBL_baseline_max_m |   PBL_suppressed_max_m |
|------:|--------------------------:|----------------------------:|------------------------------:|-------------------------------:|-------------------------:|---------------------------:|------------:|---------------------:|-----------------------:|
|     2 |                   0.23097 |                     0.68601 |                        0.8304 |                        0.4786  |                   1.8272 |                     1.0546 |     0.57715 |               2352.7 |                 3837.8 |
|     3 |                   0.21982 |                     0.66973 |                        0.9861 |                        0.6357  |                   1.7931 |                     1.0576 |     0.58982 |               2353.3 |                 4190.2 |
|     4 |                   0.21928 |                     0.65745 |                        1.1126 |                        0.76257 |                   1.7858 |                     1.0465 |     0.586   |               2663.8 |                 4901   |
|     5 |                   0.21861 |                     0.63766 |                        1.2216 |                        0.87322 |                   1.8031 |                     1.0504 |     0.58258 |               2664.3 |                 5960   |

---

