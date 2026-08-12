"""
Configuration loader and Tier 3 override machinery.

This module turns a YAML preset (Tier 2) plus optional runtime overrides
(Tier 3) into a fully validated :class:`~titan_pbl.config.schema.TitanPBLConfig`
instance.

The override mechanism is what later allows ensemble runs and Bayesian
inversion to sweep tens of thousands of parameter combinations without ever
editing a source file or a YAML preset.

Examples
--------
Load a bundled preset::

    from titan_pbl.config import load_config
    cfg = load_config("huygens_baseline")

Load a preset and override two nested parameters::

    cfg = load_config(
        "huygens_baseline",
        overrides={
            "methane.scheme": "sps",
            "regolith.porosity": 0.02,
        },
    )

Load a user file from an arbitrary path::

    cfg = load_config("/path/to/my_experiment.yaml")
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Iterable, Mapping

import yaml

from titan_pbl.config.schema import TitanPBLConfig

# Directorio donde viven los presets distribuidos con el paquete.
PRESETS_DIR: Path = Path(__file__).parent / "presets"


# ---------------------------------------------------------------------------
# Utilidades internas
# ---------------------------------------------------------------------------


def _resolve_preset_path(name_or_path: str | Path) -> Path:
    """Resolve a preset name or a filesystem path to an existing YAML file.

    Parameters
    ----------
    name_or_path
        Either a bare preset name (``"huygens_baseline"``), a preset file
        name (``"huygens_baseline.yaml"``) or an absolute/relative path.

    Returns
    -------
    Path
        Path to an existing YAML file.

    Raises
    ------
    FileNotFoundError
        If no matching preset or file can be found.
    """
    candidate = Path(name_or_path)

    # Caso 1: la ruta existe tal cual.
    if candidate.is_file():
        return candidate

    # Caso 2: nombre de preset, con o sin extension.
    stem = candidate.name
    for filename in (stem, f"{stem}.yaml", f"{stem}.yml"):
        preset = PRESETS_DIR / filename
        if preset.is_file():
            return preset

    available = sorted(p.stem for p in PRESETS_DIR.glob("*.y*ml"))
    raise FileNotFoundError(
        f"Configuration '{name_or_path}' was not found.\n"
        f"Searched as a file path and inside {PRESETS_DIR}.\n"
        f"Available presets: {available if available else '(none)'}"
    )


def _deep_merge(base: dict[str, Any], patch: Mapping[str, Any]) -> dict[str, Any]:
    """Recursively merge ``patch`` into ``base`` and return a new dictionary.

    Nested dictionaries are merged key by key; every other type is replaced
    outright. The inputs are never mutated.
    """
    result = copy.deepcopy(base)
    for key, value in patch.items():
        if (
            key in result
            and isinstance(result[key], dict)
            and isinstance(value, Mapping)
        ):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def _expand_dotted_key(dotted: str, value: Any) -> dict[str, Any]:
    """Turn ``"a.b.c"`` and a value into ``{"a": {"b": {"c": value}}}``."""
    parts = dotted.split(".")
    nested: Any = value
    for part in reversed(parts):
        nested = {part: nested}
    return nested


def _expand_overrides(overrides: Mapping[str, Any]) -> dict[str, Any]:
    """Expand a flat mapping of dotted keys into a nested dictionary.

    Both flat and already-nested entries are accepted, so the following are
    equivalent::

        {"regolith.porosity": 0.02}
        {"regolith": {"porosity": 0.02}}
    """
    expanded: dict[str, Any] = {}
    for key, value in overrides.items():
        patch = _expand_dotted_key(key, value) if "." in key else {key: value}
        expanded = _deep_merge(expanded, patch)
    return expanded


def _parse_cli_override(token: str) -> tuple[str, Any]:
    """Parse a ``key=value`` command line token into a typed pair.

    The value is interpreted with the YAML scalar rules, so ``0.02`` becomes a
    float, ``true`` becomes a bool, ``null`` becomes ``None`` and quoted text
    stays a string.
    """
    if "=" not in token:
        raise ValueError(
            f"Malformed override '{token}'. Expected the form key.subkey=value."
        )
    key, raw_value = token.split("=", 1)
    return key.strip(), yaml.safe_load(raw_value)


# ---------------------------------------------------------------------------
# API publica
# ---------------------------------------------------------------------------


def load_raw(name_or_path: str | Path) -> dict[str, Any]:
    """Load a YAML configuration file without validating it.

    Useful for inspection and debugging. Physics code should always use
    :func:`load_config` instead.
    """
    path = _resolve_preset_path(name_or_path)
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"{path} does not contain a YAML mapping at its root.")
    return data


def load_config(
    name_or_path: str | Path = "huygens_baseline",
    overrides: Mapping[str, Any] | None = None,
    cli_overrides: Iterable[str] | None = None,
) -> TitanPBLConfig:
    """Load, override and validate a simulation configuration.

    Parameters
    ----------
    name_or_path
        Preset name (``"huygens_baseline"``) or path to a YAML file.
    overrides
        Optional Tier 3 overrides. Keys may use dot notation
        (``"regolith.porosity"``) or be nested dictionaries.
    cli_overrides
        Optional iterable of ``key=value`` strings, typically coming from a
        command line interface.

    Returns
    -------
    TitanPBLConfig
        A fully validated configuration object.

    Raises
    ------
    FileNotFoundError
        If the preset cannot be located.
    pydantic.ValidationError
        If the resulting configuration violates the schema.
    """
    data = load_raw(name_or_path)

    if overrides:
        data = _deep_merge(data, _expand_overrides(overrides))

    if cli_overrides:
        parsed = dict(_parse_cli_override(token) for token in cli_overrides)
        data = _deep_merge(data, _expand_overrides(parsed))

    return TitanPBLConfig(**data)


def list_presets() -> list[str]:
    """Return the names of every preset bundled with the package."""
    return sorted(p.stem for p in PRESETS_DIR.glob("*.y*ml"))


def save_config(config: TitanPBLConfig, path: str | Path) -> Path:
    """Write a configuration back to disk as YAML.

    This guarantees that every run can be reproduced exactly: the effective
    configuration, including all applied overrides, is archived next to the
    model output.
    """
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    payload = config.model_dump(mode="json")
    with destination.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(payload, handle, sort_keys=False, allow_unicode=True)
    return destination


__all__ = [
    "load_config",
    "load_raw",
    "list_presets",
    "save_config",
    "PRESETS_DIR",
]
