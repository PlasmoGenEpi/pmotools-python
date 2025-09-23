# pmotools/__init__.py
from __future__ import annotations

try:
    # Python 3.8+
    from importlib.metadata import version, PackageNotFoundError
except Exception:  # pragma: no cover
    # Very old Pythons can fallback to pkg_resources if you ever needed it
    from pkg_resources import get_distribution as _gd  # type: ignore

    class PackageNotFoundError(Exception):
        ...

    def version(pkg: str) -> str:  # type: ignore
        try:
            return _gd(pkg).version
        except Exception as e:  # noqa: BLE001
            raise PackageNotFoundError from e


try:
    # Use the installed distribution name (matches [project].name)
    __version__ = version("pmotools")
except PackageNotFoundError:
    # When running from a source tree without being installed
    __version__ = "0+local"
