from multiprocessing import cpu_count
from pathlib import Path

import click
import pandas as pd

from structphy.bootstrapping import bootstrap_against_tree
from structphy.branch_lengths import get_mean_distance_matrix, get_stacked, get_upgma_tree
from structphy.executables import find_executable
from structphy.generate_consensus_tree import bootstrap_trees_to_consensus
from structphy.generate_matrices import generate_bootstrap_matrices_from_structures, make_fake_outgroups
from structphy.generate_trees import matrices_to_fastme_newick


@click.command()
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
@click.option("--dropout", type=str)
@click.option("--seed", type=int, help="Seed for sampling structure variants.")
def main(
    structdir: Path | None,
    fold_dir: Path | None,
    fasta: Path | None,
    dmdir: Path | None,
    outtree: Path | None,
    base_tree: Path | None,
    threads: int,
    n_bootstraps: int,
    drop_inserts: bool,
    dropout: str | None,
    n_variants: int,
    seed: int | None,
) -> None:
    """Build a phylogenetic tree from PDB structures or distance matrices."""
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
        from structphy.fasta_loading import fasta_dict_to_bootstrap_string, fasta_to_dict

        try:
            from structphy.run_inference_docker import run_esm_dropouts
        except ModuleNotFoundError as error:
            if error.name != "docker":
                raise
            raise click.ClickException("FASTA input requires `uv sync --locked --extra folding`.") from error

        try:
            dropout_rates = [float(rate) for rate in dropout.split(",")] if dropout else None
        except ValueError as error:
            raise click.BadParameter("Expected comma-separated numbers.", param_hint="--dropout") from error

        fasta_dict_full, fasta_dict_no_gaps = fasta_to_dict(fasta)
        if fold_dir is None:
            fold_dir = Path.cwd() / f"{fasta.stem}_folddir"
        fold_dir.mkdir(parents=True, exist_ok=True)
        if any(fold_dir.iterdir()):
            raise click.UsageError(f"Folding directory must be empty: {fold_dir}")

        bootstrap_fasta_path = fold_dir / "bootstraps.fasta"
        bootstrap_fasta_path.write_text(fasta_dict_to_bootstrap_string(fasta_dict_no_gaps, n_variants))
        try:
            structdir = run_esm_dropouts(
                share_dir=fold_dir,
                fasta_in=bootstrap_fasta_path,
                dropout=dropout_rates,
                max_tokens_per_batch=1300,
            )
        finally:
            bootstrap_fasta_path.unlink(missing_ok=True)

        if drop_inserts:
            conserved_dir = fold_dir.with_name(f"{fold_dir.name}_conserved")
            conserved_dir.mkdir(exist_ok=True)
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
