import os
import shutil
from pathlib import Path


def find_executable(name: str) -> str:
    """Look in the project cache first, then on PATH."""
    cache_dir = Path(os.environ.get("STRUCTPHY_CACHE_DIR", ".structphy")).resolve()
    executable = shutil.which(str(cache_dir / name)) or shutil.which(name)
    if executable is None:
        message = f"Could not find {name}. Put it in {cache_dir} or on PATH."
        raise FileNotFoundError(message)
    return str(Path(executable).resolve())
