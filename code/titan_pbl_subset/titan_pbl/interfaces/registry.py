"""
Scheme registry and plugin discovery.

The registry is what turns a scheme name in a YAML file into a class. It looks
in two places, in this order:

1. Schemes registered at runtime through :func:`register_scheme`, which is how
   tests and notebooks inject throwaway implementations.
2. Entry points published under the ``titan_pbl.<category>`` groups, which is
   how both the bundled schemes and third party packages are found.

Because the second mechanism is the standard Python packaging one, an external
contributor can publish a package that adds a new scheme and it becomes
available to every installation of the framework without any change to the
core. That is the whole point of the interface layer.

Examples
--------
Discover what is available::

    from titan_pbl.interfaces import available_schemes
    available_schemes("methane")

Resolve the scheme named by a configuration::

    from titan_pbl.interfaces import build_schemes
    schemes = build_schemes(config)
    schemes["radiation"].compute(state, dt)
"""

from __future__ import annotations

from importlib.metadata import entry_points
from typing import TYPE_CHECKING

from titan_pbl.config.schema import TitanPBLConfig
from titan_pbl.interfaces.base import NullScheme, PhysicsScheme

if TYPE_CHECKING:  # pragma: no cover
    pass

# Las cinco categorias que componen una corrida completa.
CATEGORIES: tuple[str, ...] = (
    "radiation",
    "turbulence",
    "surface",
    "methane",
    "subsurface",
)

# Registro en memoria, poblado por register_scheme.
_RUNTIME_REGISTRY: dict[str, dict[str, type[PhysicsScheme]]] = {
    category: {} for category in CATEGORIES
}


def register_scheme(
    category: str, name: str, cls: type[PhysicsScheme], overwrite: bool = False
) -> None:
    """Register a scheme class at runtime.

    Parameters
    ----------
    category
        One of :data:`CATEGORIES`.
    name
        Identifier used in configuration files.
    cls
        The scheme class, which must subclass :class:`PhysicsScheme`.
    overwrite
        Allow replacing an existing runtime registration.

    Raises
    ------
    ValueError
        If the category is unknown, the class is not a scheme, or the name is
        already taken and ``overwrite`` is false.
    """
    if category not in CATEGORIES:
        raise ValueError(
            f"Unknown scheme category '{category}'. Valid categories: {list(CATEGORIES)}"
        )
    if not (isinstance(cls, type) and issubclass(cls, PhysicsScheme)):
        raise ValueError(
            f"{cls!r} must be a subclass of PhysicsScheme to be registered."
        )
    if name in _RUNTIME_REGISTRY[category] and not overwrite:
        raise ValueError(
            f"A scheme named '{name}' is already registered under '{category}'. "
            "Pass overwrite=True to replace it."
        )
    _RUNTIME_REGISTRY[category][name] = cls


def unregister_scheme(category: str, name: str) -> None:
    """Remove a runtime registration, mostly useful for test teardown."""
    _RUNTIME_REGISTRY.get(category, {}).pop(name, None)


def _entry_point_schemes(category: str) -> dict[str, object]:
    """Return the entry points published for a category, keyed by name."""
    group = f"titan_pbl.{category}"
    try:
        points = entry_points(group=group)
    except TypeError:  # pragma: no cover - APIs antiguas de importlib
        points = entry_points().get(group, [])  # type: ignore[union-attr]
    return {point.name: point for point in points}


def available_schemes(category: str) -> list[str]:
    """List every scheme name available for a category.

    Includes both runtime registrations and published entry points, and always
    includes ``"none"``, which resolves to :class:`NullScheme`.
    """
    if category not in CATEGORIES:
        raise ValueError(
            f"Unknown scheme category '{category}'. Valid categories: {list(CATEGORIES)}"
        )
    names = set(_RUNTIME_REGISTRY[category])
    names.update(_entry_point_schemes(category))
    names.add("none")
    return sorted(names)


def get_scheme_class(category: str, name: str) -> type[PhysicsScheme]:
    """Resolve a category and name to a scheme class.

    Raises
    ------
    KeyError
        If no scheme with that name is available, with a message listing what
        is available so the user can see the typo immediately.
    """
    if category not in CATEGORIES:
        raise ValueError(
            f"Unknown scheme category '{category}'. Valid categories: {list(CATEGORIES)}"
        )

    if name == "none":
        return NullScheme

    runtime = _RUNTIME_REGISTRY[category].get(name)
    if runtime is not None:
        return runtime

    point = _entry_point_schemes(category).get(name)
    if point is not None:
        loaded = point.load()  # type: ignore[attr-defined]
        if not (isinstance(loaded, type) and issubclass(loaded, PhysicsScheme)):
            raise TypeError(
                f"Entry point '{name}' in group 'titan_pbl.{category}' resolved to "
                f"{loaded!r}, which is not a PhysicsScheme subclass."
            )
        return loaded

    raise KeyError(
        f"No scheme named '{name}' is available in category '{category}'. "
        f"Available: {available_schemes(category)}"
    )


def build_scheme(category: str, config: TitanPBLConfig) -> PhysicsScheme:
    """Instantiate the scheme a configuration selects for one category."""
    selector = {
        "radiation": config.radiation.scheme,
        "turbulence": config.turbulence.scheme,
        "surface": config.surface_layer.scheme,
        "methane": config.methane.scheme,
        "subsurface": config.subsurface.scheme,
    }
    name = selector[category]
    return get_scheme_class(category, name)(config)


def build_schemes(config: TitanPBLConfig) -> dict[str, PhysicsScheme]:
    """Instantiate every scheme a configuration selects.

    Returns
    -------
    dict
        Mapping from category to an instantiated, not yet initialised, scheme.
        The runner is responsible for calling ``initialize`` on each.
    """
    return {category: build_scheme(category, config) for category in CATEGORIES}


__all__ = [
    "CATEGORIES",
    "register_scheme",
    "unregister_scheme",
    "available_schemes",
    "get_scheme_class",
    "build_scheme",
    "build_schemes",
]
