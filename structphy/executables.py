import os
import shutil
from functools import cache
from importlib.metadata import PackageNotFoundError, distribution
from pathlib import Path


@cache
def _installed_tools() -> Path | None:
    try:
        # In an editable install __file__ is in the checkout, but the binaries
        # are installed next to the distribution's metadata in site-packages.
        return Path(str(distribution("structphy").locate_file("structphy/_bin")))
    except PackageNotFoundError:
        return None


def find_executable(name: str) -> str:
    """Find an explicitly overridden or bundled native tool."""
    override = os.environ.get("STRUCTPHY_TOOL_DIR") or os.environ.get("STRUCTPHY_CACHE_DIR")
    for directory in (Path(override) if override else None, _installed_tools()):
        if directory is None:
            continue
        executable = shutil.which(str((directory / name).absolute()))
        if executable:
            return executable

    message = f"Could not find StructPhy’s {name} executable. Reinstall StructPhy or set STRUCTPHY_TOOL_DIR."
    raise FileNotFoundError(message)
