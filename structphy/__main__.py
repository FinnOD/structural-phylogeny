from pathlib import Path
import click
import os

from structphy.install_executables import install_tmalign, install_fastme

CACHE_DIR = None
TMALIGN_URL = 'https://zhanggroup.org/TM-align/TMalign.cpp'
FASTME_URL = 'http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.6.4.tar.gz'


def setup_working_dir():
    CWD = Path(os.getcwd())
    CACHE_DIR = CWD / ".structphy"
    CACHE_DIR.mkdir(parents=False, exist_ok=True)

    install_tmalign(CACHE_DIR, TMALIGN_URL) #.structphy/TMalign
    install_fastme(CACHE_DIR, FASTME_URL) #.structphy/fastme

    
@click.command()
def main():
    setup_working_dir()
    print('heelleloeo')