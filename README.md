# Paper A: reproduction package (v2.0.1)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23747852.svg)](https://doi.org/10.5281/zenodo.23747852)

**Radiative Timescale Constraints on Boundary Layer Evolution in Models of Titan**

Guillermo Adrian Chin-Canche
Instituto Tecnológico Superior de Calkiní (ITESCAM), Tecnológico Nacional de México, Calkiní, Campeche, Mexico
gcchin@itescam.edu.mx · ORCID [0000-0003-2104-6625](https://orcid.org/0000-0003-2104-6625)

Revised version submitted to *Planetary and Space Science* (PSS-D-26-00233).

This package regenerates every number quoted in the revised manuscript, every
figure and every table, from the code in `code/` and the published model output
cited below. It contains a deliberate subset of the `titan_pbl` framework
(configuration subsystem, vertical grid, model state and the gray radiation
scheme); the full framework is under development and will be released
separately.

## What changed since v1.0.0

Version 2.0.0 accompanies the revised manuscript. The gray scheme now offers two
surface longwave boundary conditions, the original full-absorption condition and
a Kirchhoff-consistent one in which the surface reflects a fraction (1 − ε) of
the downward flux, and the heating threshold is computed exactly for both. The
package adds a one-at-a-time sensitivity study, the comparison of computed and
archived longwave fluxes in all three archive runs, the energy budget of the
suppressed-flux experiment, the readings of the output of Han and Lora (2025),
and a single versioned file with every number of the text
(`numbers/numbers_paper_A.json`). Figures and tables follow the final numbering
of the manuscript, and several values of v1.0.0 were corrected. The details, and
the list of files removed, are in [CHANGELOG.md](CHANGELOG.md). Version 1.0.0
remains available at its tag and at
[doi:10.5281/zenodo.21909362](https://doi.org/10.5281/zenodo.21909362).

## Contents

```
code/
  titan_pbl_subset/        the framework subset (pure Python; HASI profile included)
  analysis/                the analyses that produce the numbers of the text
    inputs/                documented inputs: Zenodo MD5 list of the TAM files,
                           literature quantities, transcription of the WRF source
  build_numbers_paper_A.py builds numbers/ from the analysis outputs
  make_figures.py          figures/ (with figstyle.py and analysis/make_figure_data.py)
  make_tables.py           tables/
  run_all.py               runs everything above, in order
  audit_paper_A_numbers.py verifies numbers/, figures/ and tables/
  check_columns.py         verifies that the two columns of each archive run are identical
  tests/                   tests of the subset (no external data needed)
numbers/                   numbers_paper_A.json (single source) and numbers_review.md
figures/                   fig01-fig07, figS1, figS2 (vector PDF and PNG), graphical
                           abstract (TIFF, PNG, PDF), figures_manifest.csv
tables/                    table01-table08 and the auxiliary tables aux_*.csv
data/README_DATA.md        the external data and how to obtain them
manuscript/figure_table_map.md   the script and data behind each figure and table
PROVENANCE.md              origin of the package and SHA-256 of every file
environment.yml            exact versions of the environment that produced the outputs
```

`build/` is created by `run_all.py` for intermediate files and is not versioned.

## Data

The analyses read published output; neither set of simulations is rerun.

| Data | Used for | DOI |
|---|---|---|
| Archived mtWRF output of Hartwick et al. (2026), `Dry Simulations/` | flux comparison, surface-to-air contrasts, suppressed-flux experiment, column energy closure (Sections 3.6, 4 and 5; Figures 4-7; Tables 6-8) | [10.5281/zenodo.17246644](https://doi.org/10.5281/zenodo.17246644) |
| Archived TAM output of Han and Lora (2025), version 2: the eight `atmos_diurnal_*.nc` files | surface-to-air contrasts (Table 7) | [10.5281/zenodo.16791001](https://doi.org/10.5281/zenodo.16791001) |
| Huygens HASI Level 4 atmospheric descent profile (NASA Planetary Data System) | initial state of every framework diagnostic; included in `code/titan_pbl_subset/titan_pbl/data/huygens/` | [10.17189/bcm3-4061](https://doi.org/10.17189/bcm3-4061) |
| Public WRF source code (optional) | lower boundary of the MYJ scheme (Section 5.4) | https://github.com/wrf-model/WRF |

[data/README_DATA.md](data/README_DATA.md) describes the directory layout, the
MD5 check of the TAM files and the archive runs used.

## Reproduction on Windows with conda

1. Create the environment, activate it, and install the numpy wheel from PyPI
   over the conda build (see the comment at the top of `environment.yml`):

   ```
   conda env create -f environment.yml
   conda activate titan_pbl_paper_A_v2
   python -m pip install --no-deps --force-reinstall numpy==2.4.6
   ```

2. Download the two Zenodo records (see `data/README_DATA.md`) and point the
   package at them:

   ```
   set TITAN_DATA_DIR=D:\path\to\archive
   set TITAN_TAM_DIR=D:\path\to\TAM
   ```

3. From the package root, run the tests, then the full pipeline:

   ```
   python -m pytest code/tests -q
   python code/run_all.py
   ```

   `run_all.py` deletes `build/`, runs every analysis, and rewrites `numbers/`,
   `figures/` and `tables/`. It takes about one minute on a laptop; each step
   writes its log to `build/logs/`.

4. Verify:

   ```
   python code/audit_paper_A_numbers.py
   ```

The scripts can also be run one at a time, in the order of the list `STEPS` in
`code/run_all.py`. The scripts that use `titan_pbl` import it from
`code/titan_pbl_subset/`, check that they did, and stop otherwise.

### Environment variables

| Variable | Meaning | Default |
|---|---|---|
| `TITAN_DATA_DIR` | directory that contains `Dry Simulations/` of the archive of Hartwick et al. (2026) | `<package>/data` |
| `TITAN_TAM_DIR` | directory with the eight `atmos_diurnal_*.nc` files of Han and Lora (2025) | `$TITAN_DATA_DIR/TAM Simulations_zenodo_16791001` |
| `WRF_SRC` | optional: directory with `<tag>/module_*.F` for the WRF tags V3.0, V3.6, V3.9.1.1 and v4.8.0; without it, the transcription recorded in `code/analysis/inputs/wrf_myj_transcription.json` is used | not set |
| `SKIP_ARCHIVE=1` | run only the framework diagnostics that need neither the archive nor the TAM files, and stop with a notice | not set |

### File names with colons (NTFS)

The archived output files are published as `wrfout_d01_2007-00001_09:30:00`.
Colons are not allowed in NTFS file names, so on Windows the files unpack as
`wrfout_d01_2007-00001_09_30_00`. The scripts accept either name; nothing needs
to be renamed.

### Length of the sol

Rates per Titan sol use a sol of 1 378 080 s, the value implied by the archived
integration of Hartwick et al. (2026) and the length of the Titan day in the TAM
output. It is set as `planet.sol_length` in the configuration preset
(`code/titan_pbl_subset/titan_pbl/config/presets/huygens_baseline.yaml`) and
recorded as `sol_length_s` in `numbers/numbers_paper_A.json`.

### Bit-for-bit reproduction

With the environment of step 1 on Windows, `run_all.py` reproduces every file
in `numbers/`, `figures/` and `tables/` byte for byte, which is what the audit
checks. Two things matter for that:

- the numpy wheel from PyPI: its linear algebra (bundled OpenBLAS) and that of
  the conda-forge build (MKL) give results that differ by about 1e-13 in
  relative terms, enough to change the bytes of the single source, of some
  tables and of a vector figure, although no printed value changes;
- the Arial font (present on Windows) and the FreeType bundled with the
  matplotlib wheel, for the figures.

On another system, or with other builds, the audit may report such last-digit
differences; compare the printed values (`print` in the single source) in that
case.

## Verification

`code/audit_paper_A_numbers.py`

1. regenerates the single source of numbers from the analysis outputs in
   `build/` into `build/audit/`, without touching `numbers/` (run `run_all.py`
   first, or pass `--full` to have the audit run it);
2. compares the result with `numbers/numbers_paper_A.json` key by key (full
   value, printed value, status, uncertainty) and byte for byte;
3. checks the SHA-256 of every file in `numbers/`, `figures/` and `tables/`
   against the list in `PROVENANCE.md`.

It exits with code 0 only if nothing changed, and lists every difference
otherwise.

The archived output files have two horizontal columns, because WRF allocates a
minimum stencil even for a single-column integration, and every analysis reads
the first one. `python code/check_columns.py` verifies that the two columns are
identical in the variables listed in the script, and exits with a non-zero code
otherwise.

Each entry of `numbers/numbers_paper_A.json` records the full value, the printed
form proposed by the rounding rule, its uncertainty when there is one, the
section of the manuscript where it is used, and the script and function that
produce it. `numbers/numbers_review.md` shows the same content as a table.

## How to cite

Please cite the article and this package:

> Chin-Canche, G. A. Radiative timescale constraints on boundary layer evolution
> in models of Titan. Revised version submitted to *Planetary and Space Science*
> (PSS-D-26-00233).

> Chin-Canche, G. A., 2026. Paper A reproduction package: Radiative Timescale
> Constraints on Boundary Layer Evolution in Models of Titan, version 2.0.1.
> Zenodo, [doi:10.5281/zenodo.23747852](https://doi.org/10.5281/zenodo.23747852).

If you use the data, please also cite Hartwick et al. (2026, *Planet. Sci. J.*
7, 178) and their archive (doi:10.5281/zenodo.17246644), Han and Lora (2025,
*Planet. Sci. J.* 6, 222) and their archive (doi:10.5281/zenodo.16791001), and
the HASI Level 4 product (doi:10.17189/bcm3-4061). Machine-readable metadata are
in `CITATION.cff`.

## License

BSD-3-Clause; see `LICENSE`.
