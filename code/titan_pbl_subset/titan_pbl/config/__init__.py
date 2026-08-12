"""
Configuration subsystem for titan_pbl.

Three tiers of configuration coexist:

Tier 1 -- Universal physical constants
    Defined in :mod:`titan_pbl.config.physical_constants`. Fixed by
    international convention, never tunable.

Tier 2 -- Model parameters
    Declared in :mod:`titan_pbl.config.schema` and supplied through YAML
    presets. Everything that describes Titan, the site, the regolith, the
    numerics or the schemes lives here.

Tier 3 -- Runtime overrides
    Applied by :func:`titan_pbl.config.loader.load_config` at call time.
    This is what enables ensemble sweeps and Bayesian inversion.
"""

from titan_pbl.config.loader import (
    PRESETS_DIR,
    list_presets,
    load_config,
    load_raw,
    save_config,
)
from titan_pbl.config.schema import TitanPBLConfig

__all__ = [
    "TitanPBLConfig",
    "load_config",
    "load_raw",
    "list_presets",
    "save_config",
    "PRESETS_DIR",
]
