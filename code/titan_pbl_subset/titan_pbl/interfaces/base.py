"""
The physics scheme contract.

Every physics module in the framework implements :class:`PhysicsScheme`, which
declares exactly four methods. That contract is deliberately small and
deliberately fixed: it is what allows an external contributor to publish a new
scheme as a separate package on PyPI and have the core engine pick it up
without a single line of the core changing.

The four methods
----------------
``initialize(config, state)``
    Called once before the first time step. The scheme reads its parameters
    from ``config``, allocates whatever internal storage it needs, and may
    populate diagnostic fields of ``state`` that other schemes will read on
    the first step.

``compute(state, dt)``
    Called every time step. Reads the current state, returns a
    :class:`~titan_pbl.core.state.Tendencies` object. **Must not mutate
    prognostic variables.** Diagnostic fields of the state may be written,
    because they are the scheme's own output.

``update(state, tendencies, dt)``
    Applies tendencies to the state. Most schemes inherit the default
    implementation, which delegates to the solver. Only schemes that need an
    implicit or otherwise non-additive update override it.

``diagnostics(state)``
    Returns a flat dictionary of named scalars for output and validation.

Rules that apply to every scheme
--------------------------------
1. No physical or numerical constant may appear as a literal in the scheme's
   source. Everything is read from ``config``.
2. Every equation implemented must carry a comment naming its reference,
   for example ``# Ec. 12 de Hartwick et al. 2026``.
3. ``compute`` must be free of side effects on prognostic variables, so that
   the solver may call it several times per step for multi stage integrators.
4. A scheme must never import another scheme. Schemes communicate only
   through the state.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, ClassVar

from titan_pbl.config.schema import TitanPBLConfig
from titan_pbl.core.state import ModelState, Tendencies


class PhysicsScheme(ABC):
    """Abstract base class shared by every physics module.

    Attributes
    ----------
    category
        The plugin group this scheme belongs to, for example ``"radiation"``.
        Set by each interface subclass, not by individual implementations.
    scheme_name
        Short identifier of the concrete implementation. Must match the key
        used in the ``entry_points`` declaration so that a configuration
        naming the scheme resolves to this class.
    """

    category: ClassVar[str] = "generic"
    scheme_name: ClassVar[str] = "unnamed"

    def __init__(self, config: TitanPBLConfig) -> None:
        self.config = config
        self._initialized = False

    # -----------------------------------------------------------------------
    # Contrato obligatorio
    # -----------------------------------------------------------------------

    @abstractmethod
    def initialize(self, config: TitanPBLConfig, state: ModelState) -> None:
        """Prepare the scheme before the first time step.

        Implementations must call ``super().initialize(config, state)`` or set
        ``self._initialized = True`` themselves.
        """
        self._initialized = True

    @abstractmethod
    def compute(self, state: ModelState, dt: float) -> Tendencies:
        """Compute this scheme's contribution for one time step.

        Parameters
        ----------
        state
            Current model state. Prognostic fields must be treated as read
            only; diagnostic fields owned by this scheme may be written.
        dt
            Time step [s].

        Returns
        -------
        Tendencies
            The scheme's contribution. Return an empty ``Tendencies()`` when
            the scheme is inactive rather than raising.
        """

    def update(self, state: ModelState, tendencies: Tendencies, dt: float) -> None:
        """Apply tendencies to the state.

        The default implementation does nothing, because the core solver owns
        the update for all additive schemes. Override only when the scheme
        performs its own implicit solve or an adjustment that cannot be
        expressed as a tendency, such as saturation adjustment.
        """
        return None

    @abstractmethod
    def diagnostics(self, state: ModelState) -> dict[str, float]:
        """Return named scalar diagnostics for output and validation."""

    # -----------------------------------------------------------------------
    # Utilidades comunes
    # -----------------------------------------------------------------------

    def require_initialized(self) -> None:
        """Raise if the scheme is used before :meth:`initialize` was called."""
        if not self._initialized:
            raise RuntimeError(
                f"{type(self).__name__} was used before initialize() was called. "
                "The runner must initialise every scheme before the first step."
            )

    def describe(self) -> dict[str, Any]:
        """Return metadata identifying this scheme, for provenance in output."""
        return {
            "category": self.category,
            "scheme_name": self.scheme_name,
            "class": f"{type(self).__module__}.{type(self).__qualname__}",
        }

    def __repr__(self) -> str:  # pragma: no cover - solo presentacion
        return f"<{type(self).__name__} category={self.category!r} name={self.scheme_name!r}>"


class NullScheme(PhysicsScheme):
    """A scheme that contributes nothing.

    Used whenever a category is switched off, most importantly when
    ``methane.scheme`` is ``"none"`` for the dry baseline. Having an explicit
    null implementation means the solver never needs to test whether a
    category is present.
    """

    scheme_name: ClassVar[str] = "none"

    def initialize(self, config: TitanPBLConfig, state: ModelState) -> None:
        self._initialized = True

    def compute(self, state: ModelState, dt: float) -> Tendencies:
        return Tendencies()

    def diagnostics(self, state: ModelState) -> dict[str, float]:
        return {}


__all__ = ["PhysicsScheme", "NullScheme"]
