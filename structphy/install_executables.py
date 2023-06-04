import requests
import subprocess
from sys import platform
from pathlib import Path
import os
import tarfile
import shutil

def install_tmalign(CACHE_DIR, TMALIGN_URL):

    # Define the download location for TMAlign source code.
    # Then remove that directory if it already exists.
    TMalign_dir = CACHE_DIR / 'TMalign_src'
    TMalign_cpp_filename = TMalign_dir / 'TMalign.cpp'
    try:
        shutil.rmtree(TMalign_dir)
    except FileNotFoundError:
        pass
    TMalign_dir.mkdir(parents=False, exist_ok=True)

    # Download the source code
    response = requests.get(TMALIGN_URL)
    response.raise_for_status()
    
    with open(TMalign_cpp_filename, 'w') as file:
        mac_fix = response.content.decode('utf-8')
        # Mac doesnt have malloc.h apparently. We can get it from stdlib.h
        # Got this fix from here: https://github.com/RIOT-OS/RIOT/issues/2361
        if platform == "darwin":
            mac_fix = mac_fix.replace("#include <malloc.h>", 
                """
                #if defined(__MACH__)
                #include <stdlib.h>
                #else 
                #include <malloc.h>
                #endif
                """
            )
        file.write(mac_fix)

    # Compile TMalign.cpp and make sure its executable
    # On Mac don't compile with '-static' argument.
    executable = str(CACHE_DIR / 'TMalign')
    command = ['g++', '-static' if platform != "darwin" else "", '-O3', '-ffast-math', '-lm', '-o', executable, str(TMalign_cpp_filename)]
    completed_process = subprocess.run(command, check=True, capture_output=True, text=True)
    compile_logs = completed_process.stdout + completed_process.stderr
    os.chmod(executable, 0o777)

def download_extract_tar(CACHE_DIR: Path, tar_location: Path, extract_location: Path, URL: str):
    try:
        shutil.rmtree(extract_location)
    except FileNotFoundError:
        pass
    
    # Download the fastme tarball
    response = requests.get(URL)
    response.raise_for_status()
    with open(tar_location, 'wb') as file:
        file.write(response.content)

    # Extract the tarball
    with tarfile.open(tar_location, 'r:gz') as tar:
        extracted_dir = tar.getnames()[0]
        tar.extractall(path=CACHE_DIR)
    (CACHE_DIR / extracted_dir).rename(extract_location)

def install_fastme(CACHE_DIR, FASTME_URL):

    download_extract_tar(
        CACHE_DIR = CACHE_DIR,
        tar_location = CACHE_DIR / os.path.basename(FASTME_URL), 
        extract_location = CACHE_DIR / 'FastME-master',
        URL = FASTME_URL
    )    

    # Run the configure script to generate the makefile
    configure_process = subprocess.run(['./configure'], cwd=str(CACHE_DIR / 'FastME-master'), check=True, capture_output=True, text=True)
    configure_logs = configure_process.stdout + configure_process.stderr

    # make to generate executable
    make_process = subprocess.run(['make'], cwd=str(CACHE_DIR / 'FastME-master'), check=True, capture_output=True, text=True)
    make_logs = make_process.stdout + make_process.stderr

    # Copy the compiled executable and make sure its executable
    subprocess.run(['cp', str(CACHE_DIR / 'FastME-master' / 'src' / 'fastme'), str(CACHE_DIR / 'fastme')], check=True)
    os.chmod(str(CACHE_DIR / 'fastme'), 0o777)

def install_consense(CACHE_DIR, CONSENSE_URL):
    
    download_extract_tar(
        CACHE_DIR = CACHE_DIR,
        tar_location = CACHE_DIR / os.path.basename(CONSENSE_URL), 
        extract_location = CACHE_DIR / 'phylip-master',
        URL = CONSENSE_URL
    )    

    make_process = subprocess.run(['make', '-f', 'Makefile.osx', 'consense'], cwd=str(CACHE_DIR / 'phylip-master' / 'src'), check=False, capture_output=True, text=True)
    make_logs = make_process.stdout + make_process.stderr

    subprocess.run(['cp', str(CACHE_DIR / 'phylip-master' / 'src' / 'consense'), str(CACHE_DIR / 'consense')], check=True)
    os.chmod(str(CACHE_DIR / 'consense'), 0o777)