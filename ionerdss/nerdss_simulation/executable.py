"""Locate the NERDSS executable from a user-supplied ``nerdss_path``.

Every ioNERDSS entry point that launches NERDSS takes a ``nerdss_path`` argument and
passes it through :func:`resolve_nerdss_executable`, so the rule is the same everywhere:

- a path to an executable file is run as-is;
- a directory is searched for ``nerdss`` and then ``nerdss_mpi``, first in the directory
  itself and then in its ``bin/`` subdirectory, so both a NERDSS checkout and its
  ``bin/`` folder work;
- a bare command name such as ``"nerdss"`` is looked up on ``PATH``;
- ``None`` searches ``default_root`` (when the caller has one) and then ``PATH``.
"""

import os
import shutil
import warnings
from pathlib import Path
from typing import List, Optional, Sequence, Union

PathLike = Union[str, os.PathLike]

# Serial build first: an executable launched directly runs as a single process, which is
# all a directly launched nerdss_mpi does too.
NERDSS_EXECUTABLE_NAMES = ("nerdss", "nerdss_mpi")
NERDSS_SEARCH_SUBDIRS = (".", "bin")


def _is_executable_file(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)


def _directory_candidates(directory: Path) -> List[Path]:
    return [directory / subdir / name for name in NERDSS_EXECUTABLE_NAMES for subdir in NERDSS_SEARCH_SUBDIRS]


def _search_directory(directory: Path, tried: List[str]) -> Optional[Path]:
    for candidate in _directory_candidates(directory):
        if _is_executable_file(candidate):
            return candidate
        if candidate.is_file():
            tried.append(f"{candidate} (exists but is not executable)")
        else:
            tried.append(str(candidate))
    return None


def _search_path(names: Sequence[str], tried: List[str]) -> Optional[Path]:
    for name in names:
        found = shutil.which(name)
        if found:
            return Path(found)
        tried.append(f"{name} on PATH")
    return None


def resolve_nerdss_executable(
    nerdss_path: Optional[PathLike] = None,
    *,
    default_root: Optional[PathLike] = None,
) -> Path:
    """Return the absolute path of the NERDSS executable that ``nerdss_path`` refers to.

    Args:
        nerdss_path: The NERDSS executable, a directory holding it (directly or in
            ``bin/``), or a command name to look up on ``PATH``. ``~`` is expanded.
            ``None`` searches ``default_root`` and then ``PATH``.
        default_root: Directory searched when ``nerdss_path`` is ``None``, before ``PATH``;
            :class:`~ionerdss.nerdss_simulation.Simulation` passes ``<work_dir>/NERDSS``,
            where :meth:`~ionerdss.nerdss_simulation.Simulation.install_nerdss` can put it.

    Returns:
        The absolute path of an executable file named by ``nerdss_path`` or, for a
        directory, the first of ``nerdss``, ``bin/nerdss``, ``nerdss_mpi`` and
        ``bin/nerdss_mpi`` that is executable.

    Raises:
        FileNotFoundError: If no executable is found; the message lists every location tried.
        PermissionError: If ``nerdss_path`` names a file that is not executable.
    """
    tried: List[str] = []

    if nerdss_path is None:
        found = None
        if default_root is not None:
            found = _search_directory(Path(default_root).expanduser(), tried)
        if found is None:
            found = _search_path(NERDSS_EXECUTABLE_NAMES, tried)
    else:
        raw = os.fspath(nerdss_path)
        path = Path(raw).expanduser()
        if path.is_dir():
            found = _search_directory(path, tried)
        elif path.is_file():
            if not os.access(path, os.X_OK):
                raise PermissionError(
                    f"The NERDSS executable {path} is not executable. "
                    f"Run `chmod +x {path}` or point nerdss_path at the compiled binary."
                )
            found = path
        elif os.sep not in raw and (os.altsep is None or os.altsep not in raw):
            tried.append(str(path.absolute()))
            found = _search_path([raw], tried)
        else:
            tried.append(str(path))
            found = None

    if found is None:
        locations = "\n".join(f"  - {entry}" for entry in tried)
        raise FileNotFoundError(
            "NERDSS executable not found. Looked for:\n"
            f"{locations}\n"
            "Pass nerdss_path as the NERDSS executable or as a directory containing "
            "nerdss or nerdss_mpi (directly or in bin/), or put nerdss on your PATH."
        )
    return found.absolute()


def _merge_deprecated_nerdss_dir(
    nerdss_path: Optional[PathLike],
    nerdss_dir: Optional[PathLike],
    caller: str,
) -> Optional[PathLike]:
    """Fold the deprecated ``nerdss_dir`` keyword into ``nerdss_path``."""
    if nerdss_dir is None:
        return nerdss_path
    if nerdss_path is not None:
        raise TypeError(f"{caller}() got both nerdss_path and nerdss_dir; pass only nerdss_path.")
    warnings.warn(
        f"{caller}(nerdss_dir=...) is deprecated; use nerdss_path=..., which accepts the "
        "NERDSS executable or a directory containing it (directly or in bin/).",
        DeprecationWarning,
        stacklevel=3,
    )
    return nerdss_dir
