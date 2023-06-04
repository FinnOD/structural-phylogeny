from pathlib import Path
import click
import os
import io
import multiprocessing

from structphy.install_executables import install_tmalign, install_fastme, install_consense
from structphy.generate_matrices import generate_bootstrap_matrices_from_structures
from structphy.generate_trees import matrices_to_fastme_newick
from structphy.generate_consensus_tree import bootstrap_trees_to_consensus


TMALIGN_URL = 'https://zhanggroup.org/TM-align/TMalign.cpp'
FASTME_URL = 'http://www.atgc-montpellier.fr/download/sources/fastme/fastme-2.1.6.4.tar.gz'
CONSENSE_URL = 'http://evolution.gs.washington.edu/phylip/download/phylip-3.697.tar.gz'


def setup_working_dir():
    CWD = Path(os.getcwd())
    CACHE_DIR = CWD / ".structphy"
    os.environ["STRUCTPHY_CACHE_DIR"] = str(CACHE_DIR)
    CACHE_DIR.mkdir(parents=False, exist_ok=True)

    # install_tmalign(CACHE_DIR, TMALIGN_URL) #.structphy/TMalign
    # install_fastme(CACHE_DIR, FASTME_URL) #.structphy/fastme
    # install_consense(CACHE_DIR, CONSENSE_URL) #.structphy/consense


    
@click.command()
@click.option('-d', '--structdir', type=click.Path(file_okay=False, path_type=Path))
@click.option('-f', '--fasta', type=click.File(mode='r'))
@click.option('-t', '--threads', type=int, default=multiprocessing.cpu_count())
@click.option('-n', '--n_bootstraps', type=int, default=3)
def main(structdir: Path, fasta: io.TextIOWrapper, threads: int, n_bootstraps: int):
    setup_working_dir()

    if (structdir is None) is (fasta is None): #XOR check, has to be one or the other
        raise click.UsageError('Either a directory of structures (--structdir mydir/), OR a fasta file of sequences (--fasta myseqs.fa) needs to be provided.')
    
    if fasta:
        pass
        # generate structures to .structphy/run_xxx dir
        # structdir = .structphy/run_xxx dir
 
    structure_files = [(structdir / file).resolve() for file in os.listdir(structdir) if file.endswith('.pdb')]
    bootstrap_matrices = generate_bootstrap_matrices_from_structures(structure_files, n_threads=threads, n_bootstraps=n_bootstraps)
    
    # save matrices to csv on flag

    # generate trees from matrices
    bootstrap_trees = matrices_to_fastme_newick(bootstrap_matrices, n_threads=threads)

    # generate consensus tree from bootstrap trees
    consensus_tree = bootstrap_trees_to_consensus(bootstrap_trees)
    