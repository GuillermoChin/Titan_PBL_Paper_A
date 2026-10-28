# Changelog

## [2.0.1] - 2026-10-28

### Fixed
- Corrected the Zenodo DOI. Version 2.0.0 listed doi:10.5281/zenodo.23743549,
  which was never registered; the archived version is doi:10.5281/zenodo.23747852.
  No code, data, figures, tables or numbers changed.

## [2.0.0] - 2026-10-28

Package for the revised version of the manuscript submitted to *Planetary and
Space Science* (PSS-D-26-00233). Version 1.0.0 remains available at the tag
`v1.0.0` and at [doi:10.5281/zenodo.21909362](https://doi.org/10.5281/zenodo.21909362).

Values of the revised text are not repeated here; they are in
`numbers/numbers_paper_A.json`, and the keys named below point to them.

### Scientific changes

- **Surface longwave boundary condition.** The gray scheme accepts
  `radiation.surface_lw_boundary`. `full_absorption` is the condition of
  v1.0.0 (upward surface flux ε σ Ts⁴). `kirchhoff` adds the reflected fraction
  (1 − ε) of the downward flux, as required for an opaque surface of emissivity
  ε < 1. Every framework diagnostic is computed for both conditions.
- **Exact heating threshold.** The surface-to-air contrast at which the lowest
  layer starts to warm is computed in closed form (Eq. 5 of the manuscript) for
  both conditions, and checked against the numerical zero of the scheme over
  384 combinations of profile, emissivity, infrared optical depth, weighting
  coefficient, radiation treatment and boundary condition
  (`code/analysis/r1_threshold_verification.py`,
  `tables/aux_threshold_verification_384.csv`, `code/tests/`).
- **Shortwave treatments.** Thresholds, mixing times and required contrasts are
  given for longwave only, longwave plus shortwave at local noon, and longwave
  plus the diurnal-mean shortwave term (Table 4, Figure S1).
- **Sensitivity and uncertainty.** One-at-a-time sensitivity of the required
  contrast to the surface albedo, surface emissivity, infrared optical depth,
  weighting coefficient, top-of-model infrared flux and solar constant; the
  most favorable corner of that parameter space; and the HASI temperature
  uncertainty propagated to the stratification Δθ (Table 5, Figure S2;
  `code/analysis/r2_t5_sw_dtreq.py`, `code/analysis/r3_t4_robustness.py`).
- **Flux comparison in three runs.** The longwave fluxes of the scheme,
  evaluated on the archived model state of Hartwick et al. (2026), are compared
  with the archived fluxes at all 56 interfaces and all records of sols 2–5, for
  the two baseline runs and the suppressed-flux run (Table 6, Figure 4; keys
  `flux_comparison_*`).
- **Suppressed-flux experiment.** Stability budget of the suppressed-flux run
  against the baseline at 1, 2.5 and 4 km, with the longwave contribution
  computed offline with the scheme (Figure 6); energy closure of the full
  atmospheric column in the three runs, with the archived sensible heat flux
  and with the flux entering the Mellor–Yamada–Janjic scheme (Table 8,
  Figure 7); and a reading of the lower boundary of that scheme in the public
  WRF source (`code/analysis/r2_t3_mayor3_budget.py`,
  `code/analysis/r3_t1_h33_wrf.py`).
- **Han and Lora (2025).** Surface-to-air contrasts and the lower longwave
  boundary read from version 2 of their archived TAM output, with an MD5 check
  of every file (Table 7; `code/analysis/r2_t1_tam_boundary.py`).
- **Radiative timescale.** Table 1 lists the published estimates with their
  definitions, and two evaluations with the same inputs: the bulk form
  c_p p / (g σ T³) and the linearized form c_p p / (4 g σ T³)
  (`code/analysis/r3_t5_tau_rad.py`).
- **Single source of numbers.** Every number quoted in the text is an entry of
  `numbers/numbers_paper_A.json` with its full value, the printed form proposed
  by the rounding rule, its uncertainty when there is one, the section where it
  is used, and the script and function that produce it.
  `code/audit_paper_A_numbers.py` regenerates the file and fails on any change.

### Corrections of values in v1.0.0

- **Comparison with the archived heating rate (0.90 %).** v1.0.0 compared the
  archived longwave heating rate of the lowest layer with the total (longwave
  plus shortwave) rate of the scheme, on layers of different thickness. The
  comparison is withdrawn and replaced by the flux comparison above.
- **Values from an out-of-date local copy.** Six values of the original text
  (the heating-rate profile quoted in Section 3.2 and one mixing time, quoted
  twice) had been taken from an out-of-date local copy of the outputs, not from
  the code released in v1.0.0. All values now come from the single source.
- **Maximum archive contrasts.** The contrasts given as the maxima over sols
  2–5, 0.22 K (baseline) and 0.67 K (suppressed flux), were the sol-3 maxima.
  The maxima over sols 2–5 are the keys
  `archive_contrast_max_sols2to5_baseline` and
  `archive_contrast_max_sols2to5_suppressed`.
- **Section 2.4, radiative timescale.** The values 6.4 × 10⁸ s, 469 sols and
  0.69 Titan years were truncated instead of rounded, and the length of the
  Titan year used was not stated. The formula used is the linearized form,
  a factor of 4 smaller than the bulk form. Both forms are now evaluated
  (keys `tau_rad_eq6`, `tau_rad_eq6_sols`; Table 1), and the conversion to Titan
  years uses 29.5 Earth years (Strobel et al., 2009).

### Package

- Figures follow the final numbering of the manuscript (fig01–fig07, figS1,
  figS2) as vector PDF and PNG; the graphical abstract is included as TIFF,
  PNG and PDF; `figures/figures_manifest.csv` lists the data and keys behind
  each file with its SHA-256.
- Tables 1–8 of the manuscript and eleven auxiliary tables, in CSV, with the
  source of every value.
- `code/run_all.py` runs the whole pipeline; `code/tests/` tests the subset
  without external data; `code/check_columns.py` is kept, updated for the
  new paths.
- `code/make_paper_A_figures.py` became `code/analysis/source_a_core.py` (the
  framework diagnostics); figures and tables are drawn by
  `code/make_figures.py` and `code/make_tables.py`.
- `code/audit_paper_A_numbers.py` was rewritten: it regenerates the single
  source and checks the SHA-256 of the released files listed in
  `PROVENANCE.md`.
- New environment variables `TITAN_TAM_DIR`, `WRF_SRC` and `SKIP_ARCHIVE`;
  `TITAN_DATA_DIR` keeps its meaning.
- `environment.yml` re-pinned (environment `titan_pbl_paper_A_v2`) to the
  conda-forge packages and PyPI wheels of the environment that produced the
  outputs;
  `.gitattributes` fixes the line endings; `PROVENANCE.md` records the origin
  of the package and the SHA-256 of every file.
- Documentation rewritten in English.

### Removed

- The six figures of v1.0.0 and their EPS files: `fig01_schematic`,
  `fig02_heating_profiles` (EPS; PDF and PNG replaced), `fig03_destabilisation`,
  `fig04_mixing_time`, `fig05_validation`, `fig06_timeseries`.
- The tables of v1.0.0: `table01_configuration.csv`,
  `table02_parameter_sweep.csv`, `table03_validation.csv`,
  `table04_contrasts.csv`, `table_aux_mixing_times.csv`,
  `table_aux_thresholds.csv`.
- `manuscript/captions.md` and `manuscript/equations.md`.
- `code/make_paper_A_figures.py` (see Package).

## [1.0.0] - 2026-08-12

Initial release with the original submission.
