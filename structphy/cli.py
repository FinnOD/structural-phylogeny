from multiprocessing import cpu_count
from pathlib import Path

import click
import pandas as pd

from structphy.bootstrapping import bootstrap_against_tree
from structphy.branch_lengths import get_mean_distance_matrix, get_stacked, get_upgma_tree
from structphy.executables import find_executable
from structphy.folding import FoldError, fold_fasta, parse_dropout
from structphy.generate_consensus_tree import bootstrap_trees_to_consensus
from structphy.generate_matrices import generate_bootstrap_matrices_from_structures, make_fake_outgroups
from structphy.generate_trees import matrices_to_fastme_newick


def _dropout_option(_ctx: click.Context, _param: click.Parameter, value: str | None) -> list[float] | None:
    if value is None:
        return None
    try:
        return parse_dropout(value)
    except FoldError as error:
        raise click.BadParameter(str(error)) from error


def _fold(
    fasta: Path,
    fold_dir: Path,
    n_variants: int,
    dropout: list[float] | None,
    fold_seed: int | None,
    model_dir: Path | None,
    chunk_size: int | None,
    threads: int,
) -> Path:
    click.echo(f"Folding sequences from {fasta}")
    try:
        result = fold_fasta(
            fasta,
            fold_dir,
            n_variants=n_variants,
            dropout=dropout,
            fold_seed=fold_seed,
            model_dir=model_dir,
            chunk_size=chunk_size,
            threads=threads,
        )
    except (FoldError, OSError) as error:
        raise click.ClickException(str(error)) from error
    click.echo(f"Saved PDBs and fold-manifest.json to {result}")
    return result


@click.group(invoke_without_command=True)
@click.option("-d", "--structdir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option("-f", "--fasta", type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("-dm", "--dmdir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option("-o", "--outtree", type=click.Path(dir_okay=False, path_type=Path))
@click.option(
    "--base-tree",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    help="Use a supplied Newick tree for bootstrap support.",
)
@click.option("-t", "--threads", type=click.IntRange(min=1), default=cpu_count(), show_default=True)
@click.option("-n", "--n_bootstraps", "--n-bootstraps", type=click.IntRange(min=1), default=10, show_default=True)
@click.option("--n_variants", "--n-variants", type=click.IntRange(min=1), default=10, show_default=True)
@click.option("--drop_inserts", "--drop-inserts", is_flag=True)
@click.option("--fold_dir", "--fold-dir", type=click.Path(file_okay=False, path_type=Path))
@click.option("--dropout", callback=_dropout_option, help="ESM-2 attention rate or comma-separated per-layer rates.")
@click.option("--fold-seed", type=click.IntRange(min=0, max=2**63 - 1), help="Seed for folding, separate from --seed.")
@click.option("--model-dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option("--chunk-size", type=click.IntRange(min=1), help="Chunk folding-trunk attention to reduce GPU memory.")
@click.option("--seed", type=int, help="Seed for sampling structure variants.")
@click.pass_context
def main(
    ctx: click.Context,
    structdir: Path | None,
    fold_dir: Path | None,
    fasta: Path | None,
    dmdir: Path | None,
    outtree: Path | None,
    base_tree: Path | None,
    threads: int,
    n_bootstraps: int,
    drop_inserts: bool,
    dropout: list[float] | None,
    fold_seed: int | None,
    model_dir: Path | None,
    chunk_size: int | None,
    n_variants: int,
    seed: int | None,
) -> None:
    """Build a tree from structures/matrices/FASTA, or use fold to generate PDBs."""
    if ctx.invoked_subcommand is not None:
        misplaced = [
            f"--{name.replace('_', '-')}"
            for name in ctx.params
            if ctx.get_parameter_source(name) == click.core.ParameterSource.COMMANDLINE
        ]
        if misplaced:
            raise click.UsageError(f"Put options after the fold subcommand: {', '.join(misplaced)}")
        return
    if sum(value is not None for value in (structdir, fasta, dmdir)) != 1:
        raise click.UsageError("Provide one input with --structdir, --fasta, or --dmdir.")

    try:
        find_executable("fastme")
        find_executable("consense")
        if dmdir is None:
            find_executable("TMalign")
    except FileNotFoundError as error:
        raise click.ClickException(str(error)) from error

    if fasta is not None:
        from structphy.extract_conserved_pdb import remove_inserts_from_structure
        from structphy.fasta_loading import fasta_to_dict

        if fold_dir is None:
            fold_dir = Path.cwd() / f"{fasta.stem}_folddir"
        conserved_dir = fold_dir.parent / f"{fold_dir.name}_conserved"
        if drop_inserts and (conserved_dir.exists() or conserved_dir.is_symlink()):
            raise click.UsageError(f"Refusing to reuse an existing conserved directory: {conserved_dir}")
        structdir = _fold(fasta, fold_dir, n_variants, dropout, fold_seed, model_dir, chunk_size, threads)

        if drop_inserts:
            fasta_dict_full, _ = fasta_to_dict(fasta)
            try:
                conserved_dir.mkdir()
            except OSError as error:
                raise click.ClickException(f"Could not create conserved directory {conserved_dir}: {error}") from error
            for pdb_file in fold_dir.glob("*.pdb"):
                sequence_id = pdb_file.stem.split("#")[0]
                remove_inserts_from_structure(pdb_file, fasta_dict_full[sequence_id], conserved_dir / pdb_file.name)
            structdir = conserved_dir

    if dmdir is not None:
        click.echo(f"Reading distance matrices from {dmdir}")
        matrix_files = sorted(dmdir.glob("*.csv"))
        if not matrix_files:
            raise click.UsageError(f"No CSV matrices found in {dmdir}")
        bootstrap_matrices = [pd.read_csv(filename, index_col=0) for filename in matrix_files]
    else:
        assert structdir is not None
        structure_files = sorted(path.resolve() for path in structdir.glob("*.pdb"))
        if len({path.stem.split("#")[0] for path in structure_files}) < 2:
            raise click.UsageError(f"At least two distinct PDB structure IDs are required in {structdir}")
        bootstrap_matrices = generate_bootstrap_matrices_from_structures(
            structure_files, n_threads=threads, n_bootstraps=n_bootstraps, seed=seed
        )

    mean_distance_matrix = get_mean_distance_matrix(get_stacked(bootstrap_matrices))
    mean_distance_matrix = mean_distance_matrix.rename(index=lambda name: name.replace(" ", "_"))
    mean_distance_matrix = mean_distance_matrix.rename(columns=lambda name: name.replace(" ", "_"))
    mean_distance_matrix.to_csv("mean_distance_matrix.csv")

    if dmdir is None:
        matrix_dir = Path("bootstrap_matrices")
        matrix_dir.mkdir(exist_ok=True)
        for index, matrix in enumerate(bootstrap_matrices):
            matrix.to_csv(matrix_dir / f"bootstrap_matrix_{index}.csv", float_format="%.8G")

    fake_outgroup_name = "!_OUTGROUP_!"
    make_fake_outgroups(bootstrap_matrices, fake_outgroup_name)
    bootstrap_trees = matrices_to_fastme_newick(bootstrap_matrices, n_threads=threads)
    Path("bootstrap_trees.newick").write_text("\n".join(bootstrap_trees) + "\n")

    consensus_tree = bootstrap_trees_to_consensus(bootstrap_trees, outgroup_name=fake_outgroup_name)
    Path("consensus_tree.newick").write_text(consensus_tree)

    if base_tree is None:
        base_tree_newick = get_upgma_tree(consensus_tree, mean_distance_matrix)
        Path("upgma_tree.newick").write_text(base_tree_newick)
    else:
        base_tree_newick = base_tree.read_text().strip()

    # Add bootstrap support last because ETE 3 cannot read this output format.
    bootstrapped_tree = bootstrap_against_tree(bootstrap_trees, base_tree_newick)
    (outtree or Path("boostrapped_upgma_tree.newick")).write_text(bootstrapped_tree)


@main.command("fold")
@click.option("--fasta", "-f", required=True, type=click.Path(exists=True, dir_okay=False, path_type=Path))
@click.option("--fold-dir", required=True, type=click.Path(file_okay=False, path_type=Path))
@click.option("--n-variants", type=click.IntRange(min=1), default=10, show_default=True)
@click.option(
    "--dropout",
    required=True,
    callback=_dropout_option,
    help="ESM-2 attention dropout rate or per-layer rates; use 0 for a baseline.",
)
@click.option("--fold-seed", type=click.IntRange(min=0, max=2**63 - 1))
@click.option("--model-dir", type=click.Path(exists=True, file_okay=False, path_type=Path))
@click.option("--chunk-size", type=click.IntRange(min=1), help="Chunk folding-trunk attention to reduce GPU memory.")
@click.option("--threads", type=click.IntRange(min=1), default=4, show_default=True)
def fold_command(
    fasta: Path,
    fold_dir: Path,
    n_variants: int,
    dropout: list[float],
    fold_seed: int | None,
    model_dir: Path | None,
    chunk_size: int | None,
    threads: int,
) -> None:
    """Fold sequences on the allocated CUDA GPU and stop after writing PDBs."""
    _fold(fasta, fold_dir, n_variants, dropout, fold_seed, model_dir, chunk_size, threads)
