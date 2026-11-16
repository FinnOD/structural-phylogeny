import os
import subprocess
import sys
from importlib.metadata import distribution
from pathlib import Path
from zipfile import ZipFile

import pytest

from structphy.executables import find_executable

TOOLS = ("TMalign", "fastme", "consense")


def assert_tools_work(binaries: Path, workdir: Path) -> None:
    tmalign = subprocess.run(
        [binaries / "TMalign", "-h"], cwd=workdir, capture_output=True, text=True, check=True, timeout=10
    )
    assert "TM-align (Version 20220412)" in tmalign.stdout

    (workdir / "matrix").write_text("4\nA 0 2 4 4\nB 2 0 4 4\nC 4 4 0 2\nD 4 4 2 0\n")
    subprocess.run(
        [binaries / "fastme", "-m", "NJ", "-i", "matrix", "-o", "intree"],
        cwd=workdir,
        capture_output=True,
        text=True,
        check=True,
        timeout=10,
    )
    subprocess.run(
        [binaries / "consense"], cwd=workdir, input="Y\n", capture_output=True, text=True, check=True, timeout=10
    )
    assert all(f"{taxon}:" in (workdir / "outtree").read_text() for taxon in "ABCD")


def test_uv_sync_installs_runnable_tools(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("STRUCTPHY_TOOL_DIR", raising=False)
    monkeypatch.delenv("STRUCTPHY_CACHE_DIR", raising=False)
    legacy = tmp_path / ".structphy"
    legacy.mkdir()
    for name in TOOLS:
        decoy = legacy / name
        decoy.write_text("#!/bin/sh\nexit 1\n")
        decoy.chmod(0o755)
    monkeypatch.setenv("PATH", str(legacy))

    binaries = Path(str(distribution("structphy").locate_file("structphy/_bin")))
    for name in TOOLS:
        assert Path(find_executable(name)) == binaries / name
    assert_tools_work(binaries, tmp_path)


def test_explicit_tool_directory_takes_precedence(tmp_path: Path, monkeypatch):
    override = tmp_path / "fastme"
    override.write_text("#!/bin/sh\nexit 0\n")
    override.chmod(0o755)
    monkeypatch.setenv("STRUCTPHY_TOOL_DIR", str(tmp_path))
    assert Path(find_executable("fastme")) == override


@pytest.mark.packaging
def test_uv_build_wheel_installs_runnable_tools(tmp_path: Path):
    project = Path(__file__).resolve().parents[1]
    dist = tmp_path / "dist"
    subprocess.run(["uv", "build", "--no-sources", "--out-dir", dist], cwd=project, check=True, timeout=180)
    (wheel,) = dist.glob("*.whl")
    with ZipFile(wheel) as archive:
        (wheel_metadata,) = (entry for entry in archive.namelist() if entry.endswith(".dist-info/WHEEL"))
        tag = "-".join(wheel.stem.rsplit("-", 3)[1:])
        assert not tag.endswith("-any")
        metadata = archive.read(wheel_metadata).decode().splitlines()
        assert f"Tag: {tag}" in metadata
        assert "Root-Is-Purelib: false" in metadata
        assert any(name.endswith("/extra_metadata/native/fastme-2.1.6.4.tar.gz") for name in archive.namelist())
        assert all(f"structphy/_bin/{name}" in archive.namelist() for name in TOOLS)

    venv = tmp_path / "venv"
    subprocess.run(["uv", "venv", "--python", sys.executable, venv], check=True, timeout=30)
    python = venv / "bin" / "python"
    subprocess.run(
        ["uv", "pip", "install", "--python", python, wheel],
        check=True,
        timeout=180,
    )
    env = {key: value for key, value in os.environ.items() if key not in {"STRUCTPHY_TOOL_DIR", "STRUCTPHY_CACHE_DIR"}}
    paths = subprocess.check_output(
        [
            python,
            "-I",
            "-c",
            "import sys; from structphy.executables import find_executable; "
            "print(*(find_executable(name) for name in sys.argv[1:]), sep='\\n')",
            *TOOLS,
        ],
        cwd=tmp_path,
        env=env,
        text=True,
        timeout=30,
    )
    binaries = Path(paths.splitlines()[0]).parent
    assert binaries.is_relative_to(venv)
    assert [Path(path) for path in paths.splitlines()] == [binaries / name for name in TOOLS]
    assert_tools_work(binaries, tmp_path)
