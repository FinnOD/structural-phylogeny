"""Compile the three bundled executables for Hatchling wheels."""

import gzip
import os
import platform
import shlex
import shutil
import subprocess
import sys
import sysconfig
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from hatchling.builders.hooks.plugin.interface import BuildHookInterface

TOOLS = ("TMalign", "fastme", "consense")


def compile_tools(workdir: Path, vendor: Path) -> None:
    if sys.platform not in {"darwin", "linux"}:
        message = "StructPhy builds its native tools on macOS and Linux."
        raise RuntimeError(message)

    cc = shlex.split(os.environ.get("CC", "cc"))
    cxx = shlex.split(os.environ.get("CXX", "c++"))
    cppflags = shlex.split(os.environ.get("CPPFLAGS", ""))
    ldflags = shlex.split(os.environ.get("LDFLAGS", ""))
    for command in (cc[0], cxx[0], "make"):
        if shutil.which(command) is None:
            message = f"Could not find {command}. Install a C/C++ toolchain and make before syncing StructPhy."
            raise RuntimeError(message)

    source = gzip.decompress((vendor / "TMalign-20220412.cpp.gz").read_bytes())
    if sys.platform == "darwin":
        source = source.replace(b"#include <malloc.h>", b"#include <stdlib.h>")
    (workdir / "TMalign.cpp").write_bytes(source)
    cxxflags = shlex.split(os.environ.get("CXXFLAGS", ""))
    subprocess.run(
        [*cxx, "-O3", "-ffast-math", *cppflags, *cxxflags, "-o", "TMalign", "TMalign.cpp", *ldflags],
        cwd=workdir,
        check=True,
    )

    with tarfile.open(vendor / "fastme-2.1.6.4.tar.gz") as archive:
        archive.extractall(workdir, filter="data")
    fastme = workdir / "fastme-2.1.6.4"
    subprocess.run(
        ["./configure", "--disable-OpenMP"], cwd=fastme, env={**os.environ, "CC": shlex.join(cc)}, check=True
    )
    subprocess.run(["make", "-C", "src"], cwd=fastme, check=True)
    shutil.copy2(fastme / "src" / "fastme", workdir / "fastme")

    with tarfile.open(vendor / "phylip-consense-3.697.tar.gz") as archive:
        archive.extractall(workdir, filter="data")
    subprocess.run(
        [
            *cc,
            "-std=gnu17",
            "-O2",
            "-fcommon",
            *cppflags,
            *shlex.split(os.environ.get("CFLAGS", "")),
            "-o",
            str(workdir / "consense"),
            "consense.c",
            "cons.c",
            "phylip.c",
            *ldflags,
            "-lm",
        ],
        cwd=workdir / "phylip-3.697" / "src",
        check=True,
    )


def wheel_tag() -> str:
    if sys.platform == "darwin":
        version = platform.mac_ver()[0].split(".")
        major = int(version[0])
        minor = int(version[1]) if major == 10 else 0
        return f"py3-none-macosx_{major}_{minor}_{platform.machine()}"
    return "py3-none-" + sysconfig.get_platform().replace("-", "_").replace(".", "_")


class NativeToolsHook(BuildHookInterface):
    def initialize(self, version: str, build_data: dict[str, Any]) -> None:
        self._workdir = TemporaryDirectory(prefix="structphy-build-")
        try:
            workdir = Path(self._workdir.name)
            compile_tools(workdir, Path(self.root) / "native" / "vendor")
            binaries = {str(workdir / name): f"structphy/_bin/{name}" for name in TOOLS}
            if version == "editable":
                build_data["force_include_editable"].update(binaries)
            else:
                build_data["force_include"].update(binaries)
                build_data["tag"] = wheel_tag()
            build_data["pure_python"] = False
        except BaseException:
            self._workdir.cleanup()
            raise

    def finalize(self, _version: str, _build_data: dict[str, Any], _artifact_path: str) -> None:
        self._workdir.cleanup()
