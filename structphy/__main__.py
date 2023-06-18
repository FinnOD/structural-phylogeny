from pathlib import Path
import click
import os
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

    install_tmalign(CACHE_DIR, TMALIGN_URL) #.structphy/TMalign
    install_fastme(CACHE_DIR, FASTME_URL) #.structphy/fastme
    install_consense(CACHE_DIR, CONSENSE_URL) #.structphy/consense


    
@click.command()
@click.option('-d', '--structdir', type=click.Path(file_okay=False, path_type=Path, resolve_path=True))
@click.option('-f', '--fasta', type=click.Path(exists=True,  path_type=Path, resolve_path=True))
@click.option('-t', '--threads', type=int, default=multiprocessing.cpu_count())
@click.option('-n', '--n_bootstraps', type=int, default=10)
@click.option('--n_variants', type=int, default=10)
@click.option('--drop_inserts', is_flag=True, show_default=True, default=False)
@click.option('--fold_dir', type=click.Path(file_okay=False, path_type=Path, resolve_path=True))
@click.option('--dropout', type=str)
def main(structdir: Path, fold_dir: Path, fasta: Path, threads: int, n_bootstraps: int, drop_inserts: bool, dropout: str, n_variants: int):
    setup_working_dir()

    if (structdir is None) is (fasta is None): #XOR check, has to be one or the other
        raise click.UsageError('Either a directory of structures (--structdir mydir/), OR a fasta file of sequences (--fasta myseqs.fa) needs to be provided.')
    
    if fasta:
        from structphy.run_inference_docker import run_esm_dropouts
        from structphy.fasta_loading import fasta_to_dict, fasta_dict_to_bootstrap_string
        from structphy.extract_conserved_pdb import remove_inserts_from_structure

        # Read fasta file into a dict like {id: sequence, ...}
        fasta_dict_full, fasta_dict_no_gaps = fasta_to_dict(fasta)
        
        # If a directory for the folding files isn't given make one at ./fastaname_folddir/
        if not fold_dir:
            fold_dir = Path(os.getcwd()) / (fasta.name.split('.')[0] + '_folddir')
            if not os.path.exists(fold_dir):
                os.mkdir(fold_dir)
                click.echo(f'Folding bootstraps into {fold_dir}')
        # Directory has to be empty to avoid mixing up old and new PDBs
        if os.listdir(fold_dir):
            raise click.UsageError(f'Tried to use {fold_dir}, but it wasn\'t empty.')

        # Write bootstrap fasta temporarily to the fold_dir
        bootstrap_fasta_out = fasta_dict_to_bootstrap_string(fasta_dict_no_gaps, n_variants)
        bootstrap_fasta_path = fold_dir / 'bootstraps.fasta'
        with open(bootstrap_fasta_path, 'w') as f:
            f.write(bootstrap_fasta_out)

        
        # Run inference into 'fold_dir'
        if dropout:
            dropout = [float(x) for x in dropout.split(',')]

        click.echo('Pulling and folding using docker image \'finnod/structphy-esmdropouts-openfold\'')
        structdir = run_esm_dropouts( #fasta_in and output_dir_name must be in share_dir for docker!
            share_dir=fold_dir,
            fasta_in=bootstrap_fasta_path,
            dropout=dropout,
            max_tokens_per_batch=1300, #click option or maybe auto? not sure
        )

        # Remove bootstraps.fa file in fold_dir after structures are made
        os.remove(bootstrap_fasta_path)
        structdir = fold_dir

        # if remove inserts, make a new _conserved directory then strip the inserts from the folded directory
        if drop_inserts:
            full_pdbs = [(fold_dir / file).resolve() for file in os.listdir(fold_dir) if file.endswith('.pdb')]
            out_conserved_dir = fold_dir.parent / (fold_dir.stem + '_conserved')
            out_conserved_dir.mkdir(exist_ok=True)
            for fulL_pdb_file in full_pdbs:
                pdb_name = fulL_pdb_file.name
                id = fulL_pdb_file.stem.split('#')[0]
                pdb_file_out = out_conserved_dir / pdb_name
                remove_inserts_from_structure(fulL_pdb_file,  fasta_dict_full[id], pdb_file_out)
            
            structdir = out_conserved_dir
 
    structure_files = [(structdir / file).resolve() for file in os.listdir(structdir) if file.endswith('.pdb')]

    bootstrap_matrices = generate_bootstrap_matrices_from_structures(structure_files, n_threads=threads, n_bootstraps=n_bootstraps)
    # save matrices to csv on flag

    # generate trees from matrices
    bootstrap_trees = matrices_to_fastme_newick(bootstrap_matrices, n_threads=threads)

    # generate consensus tree from bootstrap trees
    consensus_tree = bootstrap_trees_to_consensus(bootstrap_trees, fake_outgroup=True)

    # bootstrap against the consensus tree
    print(consensus_tree)

    # reweight the consensus branch lengths using distance matrices and optimise routine

    # 