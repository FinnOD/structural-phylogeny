# Structural phylogeny

StructPhy builds phylogenetic trees from protein structures or from previously generated distance matrices.

## Set up the Python environment

Install [uv](https://docs.astral.sh/uv/getting-started/installation/) 0.13 or newer, then run these commands from a checkout.

```sh
uv python install 3.12
uv sync --locked
uv run --locked structphy --help
```

Python 3.12 is supported. PhyloDM 3.2.0 uses a native extension whose source build supports Python through 3.12. On a platform without a compatible wheel, install a Rust toolchain with Cargo and a system linker before syncing. ETE 3 also builds from source via uv. The project environment lives in `.venv/`.

## Install the native tools

| Executable | When needed | Source |
| --- | --- | --- |
| `TMalign` | Building distance matrices from PDB files | [TM-align](https://zhanggroup.org/TM-align/) |
| `fastme` | Inferring a tree from each distance matrix | [FastME 2 paper and software availability](https://doi.org/10.1093/molbev/msv150) |
| `consense` | Combining the trees into a consensus | [PHYLIP](https://phylipweb.github.io/phylip/) |

Use binaries built for your OS and CPU. If building from source, you will need a C/C++ toolchain and `make`. StructPhy searches `.structphy/` under the current working directory first, then `PATH`. Set `STRUCTPHY_CACHE_DIR` to an absolute path to share a binary directory across runs. The Python lockfile covers Python packages. Record the versions or checksums of your native tools alongside any results you want to reproduce. The structure and matrix workflows have been verified on Apple Silicon macOS.

## Run an analysis

Run this example from the checkout with the native tools available. The PDB data is in the repository.

```sh
uv run --locked structphy --structdir structphy/example_data/kindom_structs \
  --n-bootstraps 2 --threads 2 --seed 42 --outtree example.newick
```

Use `--dmdir` to reuse a directory of distance-matrix CSVs. Each CSV should have row labels in the first column and matching column labels in the same order. Use `--base-tree path/to/tree.newick` to apply bootstrap support to a supplied tree with matching taxa. By default, StructPhy estimates branch lengths on the consensus tree before applying support. The seed controls PDB variant selection. Keep the input files and native tool versions fixed when comparing runs.

Intermediate files are written to the working directory, including `mean_distance_matrix.csv`, `bootstrap_trees.newick`, and `consensus_tree.newick`. Structure runs also write `bootstrap_matrices/`. To keep results together and out of Git, use a separate working directory for each run under `.local/runs/`. From there, use `uv run --project /absolute/path/to/checkout --locked structphy …`. If your native tools are in the checkout’s `.structphy/` directory, set `STRUCTPHY_CACHE_DIR` to that absolute path.

FASTA input uses the optional Docker folding workflow. Install its Python dependency with `uv sync --locked --extra folding` and run it with `uv run --locked --extra folding structphy --fasta input.fasta`. This mode requires a Docker daemon with NVIDIA GPU access and the configured folding image.

## Develop and build

```sh
uv sync --locked
uv run --locked ruff check .
uv run --locked ruff format --check .
uv run --locked ty check
uv run --locked pytest -q
uv build
```

Apply formatting with `uv run ruff format .`. Use `uv add PACKAGE` and `uv remove PACKAGE` for runtime dependencies, or pass `--group dev` for development tools and `--group research` for notebook dependencies. Run `uv lock --upgrade` when updating dependencies intentionally, then review the lockfile diff and rerun the checks. Keep `uv.lock` in version control.

Start JupyterLab with `uv run --locked --group research jupyter lab`. Install the folding extra as well if your notebook uses Docker. The distributions include the small FASTA examples. PDB datasets and notebooks are kept in the checkout.
