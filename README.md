# Structural phylogeny

StructPhy builds phylogenetic trees from protein structures or saved distance matrices.

## Quick start

To install from a checkout, use [uv](https://docs.astral.sh/uv/getting-started/installation/) 0.13 or newer, a C/C++ toolchain, and `make` on macOS or Linux. Run the example with the kinase structures in the checkout.

```sh
uv sync --locked
uv run --locked structphy --structdir structphy/example_data/kindom_structs \
  --n-bootstraps 2 --threads 2 --seed 42 --outtree example.newick
```

The sync uses Python 3.12 and installs TM-align, FastME, and PHYLIP `consense` automatically. Installing a compatible prebuilt wheel installs the same tools without compiling them. See the [build system docs](docs/build-system.md) for platform setup and release checks.

## Other inputs

- Use `--dmdir` for distance-matrix CSVs. Put row labels in the first column and use matching column labels in the same order.
- Use `--base-tree` to apply bootstrap support to an existing Newick tree with matching taxa.
- For FASTA input on a Linux NVIDIA GPU node, install `uv sync --locked --extra fold`. Use `structphy fold` to write PDBs, or retain `structphy --fasta` for the combined folding and tree workflow. See [HPC folding](docs/folding.md).

StructPhy writes intermediate files to the working directory. Run each analysis in its own directory under the Git-ignored `runs/`, for example `runs/my-analysis/`. From there, use `uv run --project /path/to/checkout --locked structphy …` with paths relative to the run directory. See `uv run --locked structphy --help` for all options.

## Folding on an HPC

Folding runs directly in this repository's Python environment with the optional
`structphy[fold]` dependencies:

```sh
uv sync --locked --extra fold
uv run --locked --extra fold structphy fold \
  --fasta reference.fasta --fold-dir runs/reference-folds \
  --n-variants 10 --dropout 0.1 --fold-seed 42
```

Run in an allocated CUDA GPU job. Choose the dropout rate for your experiment;
use zero for a baseline. Copy the resulting PDB directory to the Mac and use
`--structdir` with the base install. The [folding guide](docs/folding.md)
covers scheduler setup, offline checkpoints, reproducibility, and the existing
FASTA parser limitations.
