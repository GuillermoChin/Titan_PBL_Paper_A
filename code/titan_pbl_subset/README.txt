Subset of the titan_pbl framework: only the modules needed to reproduce the
numbers, figures and tables of Paper A (configuration subsystem, vertical grid,
model state and the gray radiation scheme). The full framework is under
development and will be released separately.

The gray scheme offers two surface longwave boundary conditions, selected with
radiation.surface_lw_boundary:

    full_absorption   F_up(0) = eps sigma Ts^4
    kirchhoff         F_up(0) = eps sigma Ts^4 + (1 - eps) F_down(0)

The preset huygens_baseline uses "full_absorption" by default.

The scripts in code/ and code/analysis/ put this directory first on sys.path
and stop if titan_pbl is imported from anywhere else. To use the subset
interactively, from the package root:

    set PYTHONPATH=%CD%\code\titan_pbl_subset
    python -c "import titan_pbl; print(titan_pbl.__file__)"

The Huygens HASI Level 4 descent profile (NASA PDS, doi:10.17189/bcm3-4061) is
included in titan_pbl/data/huygens/. Code comments are partly in Spanish, and
some cite entries H-nn of the author's development log, which is not part of
this package.
