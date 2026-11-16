import shutil
import subprocess
import tarfile
from contextlib import suppress
from pathlib import Path
from sys import platform
from urllib.parse import urlsplit

import requests


def install_tmalign(CACHE_DIR: Path, TMALIGN_URL: str) -> None:
    source_dir = CACHE_DIR / "TMalign_src"
    with suppress(FileNotFoundError):
        shutil.rmtree(source_dir)
    source_dir.mkdir(parents=True)

    response = requests.get(TMALIGN_URL, timeout=60)
    response.raise_for_status()
    source = response.text
    if platform == "darwin":
        source = source.replace("#include <malloc.h>", "#include <stdlib.h>")
    source_file = source_dir / "TMalign.cpp"
    source_file.write_text(source)

    command = ["g++"]
    if platform != "darwin":
        command.append("-static")
    command.extend(["-O3", "-ffast-math", "-lm", "-o", str(CACHE_DIR / "TMalign"), str(source_file)])
    subprocess.run(command, check=True)


def download_extract_tar(CACHE_DIR: Path, tar_location: Path, extract_location: Path, URL: str) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with suppress(FileNotFoundError):
        shutil.rmtree(extract_location)

    response = requests.get(URL, timeout=60)
    response.raise_for_status()
    tar_location.write_bytes(response.content)
    with tarfile.open(tar_location, "r:gz") as archive:
        extracted_dir = archive.getnames()[0]
        archive.extractall(path=CACHE_DIR, filter="data")
    (CACHE_DIR / extracted_dir).rename(extract_location)


def install_fastme(CACHE_DIR: Path, FASTME_URL: str) -> None:
    source_dir = CACHE_DIR / "FastME-master"
    download_extract_tar(CACHE_DIR, CACHE_DIR / Path(urlsplit(FASTME_URL).path).name, source_dir, FASTME_URL)
    subprocess.run(["./configure"], cwd=source_dir, check=True)
    subprocess.run(["make"], cwd=source_dir, check=True)
    shutil.copy(source_dir / "src" / "fastme", CACHE_DIR / "fastme")


def install_consense(CACHE_DIR: Path, CONSENSE_URL: str) -> None:
    source_dir = CACHE_DIR / "phylip-master"
    download_extract_tar(CACHE_DIR, CACHE_DIR / Path(urlsplit(CONSENSE_URL).path).name, source_dir, CONSENSE_URL)
    makefile = "Makefile.osx" if platform == "darwin" else "Makefile.unx"
    subprocess.run(["make", "-f", makefile, "consense"], cwd=source_dir / "src", check=True)
    shutil.copy(source_dir / "src" / "consense", CACHE_DIR / "consense")
