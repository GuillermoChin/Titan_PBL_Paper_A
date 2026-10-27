"""
Tests of the titan_pbl subset shipped with the package (no archive, no full framework).

- the exact heating threshold (Eq. 5) against the numerical zero of the
  scheme over the 384 combinations of Section 2.3;
- the two surface boundary conditions: zero heating of an isothermal column
  with the Kirchhoff-consistent condition, the full-absorption threshold
  Ta (eps^-1/4 - 1) in the same column, and the reflected term (1 - eps) F_down
  between the two upward surface fluxes;
- the configuration rejects an unknown boundary condition.

    python -m pytest code/tests -q
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "code" / "titan_pbl_subset"))
sys.path.insert(0, str(ROOT / "code" / "analysis"))

import r1_threshold_verification as T  # noqa: E402
import titan_pbl  # noqa: E402
from titan_pbl.config import load_config  # noqa: E402
from titan_pbl.config.physical_constants import STEFAN_BOLTZMANN  # noqa: E402
from titan_pbl.core import build_grid, initialize_state  # noqa: E402
from titan_pbl.schemes.radiation import GrayRadiation  # noqa: E402

T_ISO = 93.65


def test_subset_is_the_one_shipped():
    assert (ROOT / "code" / "titan_pbl_subset") in Path(titan_pbl.__file__).resolve().parents


def test_exact_threshold_matches_numerical_zero_384_combinations():
    combos = list(itertools.product(T.PROFILES, T.EMISSIVITIES, T.TAUS, T.WEIGHTS, T.RADIATION, T.BOUNDARIES))
    assert len(combos) == 384
    worst = 0.0
    for (profile, i_ir), eps, tau, f, radiation, boundary in combos:
        r = T.evaluate(profile, i_ir, eps, tau, f, radiation, boundary)
        assert not r["note"].startswith("DISCREPANCIA"), r
        if np.isfinite(r["abs_num_minus_exact_K"]):
            worst = max(worst, r["abs_num_minus_exact_K"])
    # El texto (Seccion 2.3) dice "within about 1e-12 K"; el criterio declarado es 1e-6 K.
    assert worst <= T.TOL_NUM_EXACT
    assert worst < 1e-10


def _isothermal(boundary: str, eps: float = 0.95):
    cfg = load_config("huygens_baseline", overrides={
        "radiation.surface_lw_boundary": boundary,
        "radiation.surface_emissivity": eps,
        "atmosphere.initialization.profile_source": "isothermal",
        "atmosphere.initialization.isothermal_temperature": T_ISO,
        "radiation.ir_top_flux": STEFAN_BOLTZMANN * T_ISO**4,
        "radiation.solar_constant": 0.0,
    })
    state = initialize_state(cfg, build_grid(cfg))
    scheme = GrayRadiation(cfg)
    scheme.initialize(cfg, state)
    return cfg, state, scheme


def test_kirchhoff_isothermal_column_has_zero_heating():
    _, state, scheme = _isothermal("kirchhoff")
    state.skin_temperature = float(state.atmospheric_temperature[0])
    assert float(np.max(np.abs(scheme.heating_rate(state)))) <= 1e-12


@pytest.mark.parametrize("eps", [0.90, 0.95])
def test_full_absorption_isothermal_threshold(eps):
    r = T.evaluate("isothermal", None, eps, 10.0, 0.15, "lw", "full_absorption")
    expected = r["Ta_K"] * (eps ** -0.25 - 1.0)
    assert abs(r["dT_exact_K"] - expected) < 1e-9
    assert abs(r["dT_num_K"] - expected) < 1e-9
    r_k = T.evaluate("isothermal", None, eps, 10.0, 0.15, "lw", "kirchhoff")
    assert abs(r_k["dT_exact_K"]) < T.TOL_ISO_KIRCHHOFF


def test_reflected_term_between_boundary_conditions():
    eps = 0.95
    ups = {}
    for b in ("full_absorption", "kirchhoff"):
        cfg = load_config("huygens_baseline", overrides={"radiation.surface_lw_boundary": b,
                                                         "radiation.surface_emissivity": eps})
        state = initialize_state(cfg, build_grid(cfg))
        scheme = GrayRadiation(cfg)
        scheme.initialize(cfg, state)
        down, up = scheme.longwave_fluxes(state)
        ups[b] = (float(up[0]), float(down[0]))
    reflected = ups["kirchhoff"][0] - ups["full_absorption"][0]
    assert reflected == pytest.approx((1.0 - eps) * ups["kirchhoff"][1], rel=1e-12)


def test_unknown_boundary_condition_is_rejected():
    with pytest.raises(Exception):
        load_config("huygens_baseline", overrides={"radiation.surface_lw_boundary": "reflecting"})
