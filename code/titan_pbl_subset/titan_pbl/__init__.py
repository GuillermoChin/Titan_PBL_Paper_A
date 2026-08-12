"""
titan_pbl -- A modular 1D framework for modelling Titan's Planetary Boundary Layer.

Version 0.1.0 targets a faithful replication of

    Hartwick, V. L., Rafkin, S. C. R. & Soto, A. (2026),
    "Modeling Titan's Planetary Boundary Layer at the Huygens Landing Site:
    Evidence for the Role of Subsurface Methane Evaporation in Limiting
    Convective Growth", The Planetary Science Journal 7:178.
    DOI 10.3847/PSJ/ae809d

Design principle
----------------
No physical or numerical parameter is hardcoded in the source. Everything is
declared in a validated configuration schema and supplied through YAML
presets, so that any parameter can be swept, inverted or overridden at runtime
without touching the code.
"""

__version__ = "0.1.0"
__author__ = "Guillermo Adrian Chin Canche"

__all__ = ["__version__", "__author__"]
