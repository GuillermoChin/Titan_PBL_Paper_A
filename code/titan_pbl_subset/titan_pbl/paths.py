"""Resolution of data file paths declared in configuration presets.

Presets refer to data files with repository-relative paths such as
``data/huygens/HASI_L4_ATMO_PROFILE_DESCEN_TAB_v2.txt``. Resolving those
against the current working directory alone makes a run depend on where the
interpreter happened to be started, which is exactly the kind of hidden state
that breaks reproducibility: the same preset would load in one shell and fail
in another.

This module resolves such paths against an explicit, ordered list of roots, so
that a preset means the same file regardless of the working directory.

Search order
------------
1. ``$TITAN_PBL_DATA_ROOT``, when set. Escape hatch for machines that keep the
   observational archive outside the source tree.
2. The current working directory. Preserves the historical behaviour, so a
   deliberate local override still wins over the bundled copy.
3. The installed package directory, ``titan_pbl/``. This is the copy shipped
   with the distribution and the one that makes the test suite independent of
   the working directory.
4. The repository root, the parent of the package directory. This is where the
   full external archive lives in a source checkout.

Absolute paths are returned untouched: an author who writes an absolute path
means that file and nothing else.
"""

from __future__ import annotations

import os
from pathlib import Path

# Directorio del paquete instalado y raiz del checkout que lo contiene.
PACKAGE_ROOT: Path = Path(__file__).resolve().parent
REPOSITORY_ROOT: Path = PACKAGE_ROOT.parent

# Nombre de la variable de entorno que antepone una raiz de busqueda propia.
DATA_ROOT_ENV_VAR = "TITAN_PBL_DATA_ROOT"


def data_search_roots() -> list[Path]:
    """Return the ordered roots against which relative data paths are resolved.

    The environment variable is read on every call rather than cached, so that
    a test or a notebook can point the resolver at a different archive without
    reimporting the package.
    """
    roots: list[Path] = []

    override = os.environ.get(DATA_ROOT_ENV_VAR)
    if override:
        roots.append(Path(override))

    roots.extend([Path.cwd(), PACKAGE_ROOT, REPOSITORY_ROOT])

    # Se eliminan duplicados conservando el orden: en una instalacion editable
    # el cwd puede coincidir con la raiz del repositorio.
    unique: list[Path] = []
    for root in roots:
        resolved = root.resolve()
        if resolved not in unique:
            unique.append(resolved)
    return unique


def resolve_data_path(path: str | Path, must_exist: bool = True) -> Path:
    """Resolve a data file path declared in a preset or passed by a caller.

    Parameters
    ----------
    path
        Absolute path, or path relative to one of :func:`data_search_roots`.
    must_exist
        When true, raise if the file is found under no root. When false, return
        the candidate under the first root so that the caller can produce its
        own error message.

    Returns
    -------
    Path
        Absolute path to an existing file when ``must_exist`` is true.

    Raises
    ------
    FileNotFoundError
        If ``must_exist`` is true and no root contains the file. The message
        lists every root searched, because the usual cause is a working
        directory the user did not expect.
    """
    candidate = Path(path)

    if candidate.is_absolute():
        if must_exist and not candidate.is_file():
            raise FileNotFoundError(f"Data file not found: {candidate}")
        return candidate

    roots = data_search_roots()
    for root in roots:
        resolved = root / candidate
        if resolved.is_file():
            return resolved

    if not must_exist:
        return roots[0] / candidate

    listed = "\n".join(f"  - {root}" for root in roots)
    raise FileNotFoundError(
        f"Data file '{candidate}' was not found under any search root:\n{listed}\n"
        f"Set {DATA_ROOT_ENV_VAR} to point at the directory that contains it."
    )


__all__ = [
    "PACKAGE_ROOT",
    "REPOSITORY_ROOT",
    "DATA_ROOT_ENV_VAR",
    "data_search_roots",
    "resolve_data_path",
]
