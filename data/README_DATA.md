# External data

The package reads two published archives that are not redistributed here, and
one observational profile that is. This file says where to get each one, where
to put it, and which parts are used.

## 1. Archive of Hartwick et al. (2026)

Hartwick, V. L., Rafkin, S. C. R., Soto, A. Simulation archive for *Modeling
Titan's planetary boundary layer at the Huygens landing site*, Zenodo,
[doi:10.5281/zenodo.17246644](https://doi.org/10.5281/zenodo.17246644).

The record contains four directories of mtWRF output (`Dry Simulations/`,
`SPS Simulations/`, `NSMT Simulations/`, `MRPS Simulations/`). Only
`Dry Simulations/` is read, and from it three runs:

| Run | Role in the manuscript | Prescribed deep subsurface temperature (TMN) |
|---|---|---|
| `1D_Huygens_TMN91_5sol` | baseline | 91.983 K |
| `1D_Huygens_baseline_5sol` | second baseline | 92.964 K |
| `1D_Huygens_0.1flhc_5sol` | suppressed-flux experiment | 91.865 K |

The archive does not label one dry run as the baseline of the article. The run
`1D_Huygens_TMN91_5sol` is used as the baseline because, apart from the scaling
of the exchange coefficient, its metadata differ from those of the
suppressed-flux run only in TMN; the second baseline is reported alongside it
wherever the archive is read (Tables 6 and 8, Figures 4 and 7).

Each run directory holds one output file, `wrfout_d01_2007-00001_09:30:00`.

**Layout.** Unpack the record so that `Dry Simulations/` is a direct child of a
directory of your choice, and set

```
set TITAN_DATA_DIR=D:\path\to\that\directory
```

Without `TITAN_DATA_DIR`, the scripts look in `data/` of the package.

**File names on Windows.** Colons are not allowed in NTFS file names, so on
Windows the output files unpack as `wrfout_d01_2007-00001_09_30_00`. The scripts
accept both names; nothing needs to be renamed.

## 2. Archive of Han and Lora (2025), version 2

Han, S., Lora, J. M. Archived Titan Atmospheric Model (TAM) output for *Thermal
inertia controls on Titan's surface temperature and planetary boundary layer
structure*, Zenodo, version 2,
[doi:10.5281/zenodo.16791001](https://doi.org/10.5281/zenodo.16791001).

Only the eight diurnal files are read:

```
atmos_diurnal_case1_huygens.nc   ...   atmos_diurnal_case4_huygens.nc
atmos_diurnal_case1_dragonfly.nc ...   atmos_diurnal_case4_dragonfly.nc
```

Put them in one directory and set

```
set TITAN_TAM_DIR=D:\path\to\TAM
```

(default: `$TITAN_DATA_DIR/TAM Simulations_zenodo_16791001`).

**Version check.** Before reading them, the analyses compare the MD5 of each
file with the value published by Zenodo for version 2 of the record, listed in
`code/analysis/inputs/tam_16791001_md5.txt`, and stop if any differs. Version 2
adds the variables `diff_t` and `diff_m`; the variables read here are listed in
the docstring of `code/analysis/r2_t1_tam_boundary.py`.

## 3. HASI Level 4 descent profile (included)

Huygens HASI Level 4 atmospheric descent profile, NASA Planetary Data System,
[doi:10.17189/bcm3-4061](https://doi.org/10.17189/bcm3-4061).

The profile initializes every framework diagnostic and is included in the
package at
`code/titan_pbl_subset/titan_pbl/data/huygens/HASI_L4_ATMO_PROFILE_DESCEN_TAB_v2.txt`.
The PDS distributes it as `HASI_L4_ATMO_PROFILE_DESCEN.TAB`; the file was
renamed so that its extension is not read as `.TAB_v2`. Its content is
unchanged.

## 4. Public WRF source (optional)

Section 5.4 of the manuscript examines the lower boundary of the
Mellor–Yamada–Janjic scheme in the public WRF source
(https://github.com/wrf-model/WRF). `code/analysis/r3_t1_h33_wrf.py` reads
`module_sf_myjsfc.F`, `module_sf_slab.F`, `module_surface_driver.F`,
`module_bl_myjpbl.F` and `module_pbl_driver.F` (directory `phys/` of the
repository) for the tags V3.0, V3.6, V3.9.1.1 and v4.8.0, from

```
set WRF_SRC=D:\path\to\wrf_files
```

with one subdirectory per tag (`%WRF_SRC%\V3.6\module_sf_myjsfc.F`, ...). If
`WRF_SRC` is not set, the script uses the lines it recorded from those files,
with their commit, file and line number, in
`code/analysis/inputs/wrf_myj_transcription.json`.

## Running without the archives

```
set SKIP_ARCHIVE=1
python code/run_all.py
```

runs only the framework diagnostics that need neither archive and ends with a
notice. The single source of numbers, the figures and the tables need both
archives.
