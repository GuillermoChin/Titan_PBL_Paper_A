"""
Physics scheme contracts.

This subpackage contains no physics at all: it declares the abstract base
classes that every scheme implements, plus the registry that resolves a scheme
name in a configuration file to a class.

Keeping the contracts separate from the implementations is what allows an
external package to add a scheme without touching the core, and what allows
the core solver to compose schemes without knowing anything about them.
"""

from titan_pbl.interfaces.base import NullScheme, PhysicsScheme
from titan_pbl.interfaces.methane import MethaneScheme
from titan_pbl.interfaces.radiation import RadiationScheme
from titan_pbl.interfaces.registry import (
    CATEGORIES,
    available_schemes,
    build_scheme,
    build_schemes,
    get_scheme_class,
    register_scheme,
    unregister_scheme,
)
from titan_pbl.interfaces.subsurface import SubsurfaceScheme
from titan_pbl.interfaces.surface import SurfaceLayerScheme
from titan_pbl.interfaces.turbulence import TurbulenceScheme

__all__ = [
    "PhysicsScheme",
    "NullScheme",
    "RadiationScheme",
    "TurbulenceScheme",
    "SurfaceLayerScheme",
    "MethaneScheme",
    "SubsurfaceScheme",
    "CATEGORIES",
    "register_scheme",
    "unregister_scheme",
    "available_schemes",
    "get_scheme_class",
    "build_scheme",
    "build_schemes",
]
