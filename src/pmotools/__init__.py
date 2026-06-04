# pmotools/__init__.py
from __future__ import annotations

from importlib.metadata import version, PackageNotFoundError

#: Version of the PMO schema this package targets. Single source of truth;
#: bump when the schema changes. Independent of the package release version.
__schema_version__ = "1.1.0"

try:
    # Distribution version from installed metadata (matches [project].version)
    __version__ = version("pmotools")
except PackageNotFoundError:
    # Running from a source tree without being installed
    __version__ = "0+local"


def get_pmotools_version() -> str:
    """Return the installed pmotools package version (e.g. '1.0.0')."""
    return __version__


def get_pmo_schema_version() -> str:
    """Return the PMO schema version this package targets (e.g. '1.0.0')."""
    return __schema_version__


__all__ = [
    "__version__",
    "__schema_version__",
    "get_pmotools_version",
    "get_pmo_schema_version",
]
