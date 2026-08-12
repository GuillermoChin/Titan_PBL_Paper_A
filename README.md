# Paper A — reproduction package

**Radiative Timescale Constraints on Boundary Layer Evolution in Models of Titan**

Guillermo Adrián Chin Canché
Instituto Tecnológico Superior de Calkiní, Tecnológico Nacional de México, Calkiní, Campeche, Mexico
gcchin@itescam.edu.mx
ORCID:0000-0003-2104-6625

Submitted to *Planetary and Space Science*.

This package reproduces every figure and table in the manuscript. It is a
deliberate subset of the `titan_pbl` framework: only the modules needed for this
paper are included. The full framework is released separately.

## Contents

```
manuscript/     paper_A.docx, captions.md, figure_table_map.md
figures/        six figures, each as PNG (600 dpi), PDF and EPS
tables/         four numbered tables plus two auxiliary ones, as CSV
code/           the scripts, and titan_pbl_subset/ with the modules they need
data/           where the external archive must be placed
environment.yml exact pinned versions of the environment that produced these
```

## Two data sources, kept separate

**Source A** is the framework in `code/titan_pbl_subset/`. It produces Figures 1
to 4 and Tables 1 and 2. Every one of those is a diagnostic evaluated on an
initialised atmospheric state plus the radiative transfer scheme: **no time
integration is performed**. [doi:10.5281/zenodo.21909362](https://doi.org/10.5281/zenodo.21909362)

**Source B** is the published output of Hartwick et al. (2026), Zenodo
[doi:10.5281/zenodo.17246644](https://doi.org/10.5281/zenodo.17246644). It
produces Figures 5 and 6 and Tables 3 and 4. Those are read from their archive,
not reproduced.

## Reproduction

### 1. Environment

```
conda env create -f environment.yml
conda activate titan_pbl_paper_A
```

### 2. External archive (only needed for Source B)

Download the Zenodo archive and unpack it so that `Dry Simulations/` is a direct
child of a directory of your choosing, then point the script at it:

```
set TITAN_DATA_DIR=D:\path\to\archive
```

On Windows the published file names contain colons, which NTFS forbids, so they
unpack as `wrfout_d01_2007-00001_09_30_00`. The script handles both forms.

Without the archive, Source A still runs and Source B is skipped with a notice.

### 3. Generate

```
cd code
python make_paper_A_figures.py
```

Outputs land in `code/paper_A/`. Compare them against the reference copies in
`figures/` and `tables/`.

### 4. Verify

```
python check_columns.py            # the single-column assumption of the archive
python audit_paper_A_numbers.py    # every number in the text against its CSV
```

`audit_paper_A_numbers.py` exits non-zero if any number quoted in the manuscript
has drifted from the output that produced it. It checks 19 quantities.

## The sol length

Rates per Titan solar day use a sol of 1 378 080 s, the value implied by the
archived integration (86 400 steps per sol at the archived time step). This makes
our rates and the archived rates directly comparable.

## Citing

See `CITATION.cff`. The archive of Hartwick et al. (2026) must be cited
separately when Source B results are used.

## License

BSD-3-Clause, see `LICENSE`.
