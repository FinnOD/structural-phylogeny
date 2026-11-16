# Folding on an HPC

StructPhy folds proteins directly in Python with the Hugging Face Transformers
port of ESMFold. Install the `fold` extra in this repository's uv environment.
The base install remains suitable for analysing saved PDBs on a Mac.

The tested backend is `facebook/esmfold_v1` at checkpoint revision
`75a3841ee059df2bf4d56688166c8fb459ddd97a`, loaded through Transformers 4.57.2
and PyTorch 2.7.1. This replaces the old fork and container runtime; it does not
change to a different folding model such as ESMFold2.

## Installation

A checkout requires Python 3.12, a C/C++ compiler and make for the bundled tree
tools. uv can use its managed Python; the system Python need not be changed.
On the reviewed cluster, `module load gcc/12.3.0` supplies the compiler.

Start without PyTorch or Transformers modules in the shell: those modules inject
packages through `PYTHONPATH` and can override the uv environment.

```sh
module load gcc/12.3.0
uv sync --locked --extra fold
```

On Linux x86-64, the checkout selects PyTorch's CUDA 12.6 wheel index explicitly.
The NVIDIA driver must support that runtime. A separately loaded CUDA toolkit
does not choose the PyTorch wheel; the folding code does not compile CUDA
extensions. Other CUDA configurations require their own compatibility check.
See [uv with PyTorch](https://docs.astral.sh/uv/guides/integration/pytorch/)
and [PyTorch wheel sources](https://pytorch.org/get-started/previous-versions/).

For a built wheel, the `[fold]` extra travels in package metadata, but uv's
checkout source configuration does not. On Linux x86-64 use:

```sh
uv pip install 'dist/structphy-0.0.1-<platform>.whl[fold]' \
  --torch-backend=cu126
```

Replace the wheel filename with the actual built wheel. uv directs PyTorch to
the CUDA index while resolving ordinary dependencies from PyPI. Use the same CUDA
source when installing `structphy[fold]` from another package source.
Keep `--extra fold` on subsequent `uv run` commands: plain `uv sync` or
`uv run` synchronizes the base environment and can remove optional packages.

## Fold without building a tree

Run within an allocated NVIDIA GPU job:

```sh
uv run --locked --extra fold structphy fold \
  --fasta reference.fasta --fold-dir runs/reference-folds \
  --n-variants 10 --dropout 0.1 --fold-seed 42 --threads 4
```

Here 0.1 illustrates an explicit experimental rate; it is not a validated
scientific default. Use `--dropout 0` for a baseline. The command also accepts
36 comma-separated rates, one per ESM-2 attention layer. Values must be finite
and in `[0, 1)`. A scalar broadcasts to all layers. Wrong schedule lengths
fail before loading model weights.

The model loads once and predicts sequences and variants serially. Inference
uses float32 and eager attention with TF32 disabled. Only the ESM-2
self-attention modules enable dropout; the rest of the model stays in
evaluation mode and gradients remain disabled. This policy differs from the
old runner's whole-model training mode when dropout was supplied, so results
must be treated as a method change.

`--fold-seed` controls folding. If omitted, a generated seed is recorded.
Each sequence/variant gets a stable derived seed, recorded in the manifest.
The tree workflow's `--seed` independently controls variant sampling.
Reproducibility comparisons should use the same device, software and execution
settings and compare parsed coordinates within a stated tolerance.

For longer inputs, `--chunk-size 64` (or a smaller positive value) can reduce
folding-trunk attention memory. It does not guarantee that every sequence fits.
No batching, mixed precision, multi-GPU coordination or resume mode is enabled.

## Scheduler example

This is a Slurm template; adapt the partition, account, memory and time to the
site and workload. Request an H100-capable partition if that hardware is needed.

```sh
#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:30:00

set -euo pipefail
cd /path/to/structural-phylogeny
uv run --locked --extra fold structphy fold \
  --fasta /path/to/reference.fasta \
  --fold-dir /path/to/shared/results/reference-folds \
  --n-variants 10 --dropout 0.1 --fold-seed 42 \
  --threads "$SLURM_CPUS_PER_TASK"
```

Do not override the scheduler's `CUDA_VISIBLE_DEVICES`. StructPhy uses the
current logical CUDA device, normally device 0 within the allocation. It does
not request access to all GPUs on the host. For the existing FASTA-to-tree
workflow, also pass `--threads "$SLURM_CPUS_PER_TASK"` to limit CPU processes.

## Checkpoints and offline nodes

The first run downloads about 7.9 GiB of weights. Set `HF_HOME` to an
appropriate persistent shared cache before downloading; avoid a small home
directory quota. Keep the same cache configuration in compute jobs.

Stage the pinned snapshot on a node with outbound access:

```sh
export HF_HOME=/path/to/shared/scratch/huggingface
uv run --locked --extra fold python -c \
  "from huggingface_hub import snapshot_download; snapshot_download('facebook/esmfold_v1', revision='75a3841ee059df2bf4d56688166c8fb459ddd97a', allow_patterns=['config.json', 'pytorch_model.bin'])"
```

On a compute node, set `HF_HUB_OFFLINE=1` to use that cached revision without
network access. Alternatively, supply `--model-dir /path/to/snapshot`.
A supplied directory is loaded with local-files-only behavior; its revision is
recorded as null because StructPhy cannot verify the provenance of arbitrary
local files. Stage the specified pinned snapshot and retain its provenance.

Transformers may warn about missing contact-head regression weights. The
folding forward path does not call the contact-prediction head; this warning
was also observed during the GPU validation.

## Outputs and tree analysis

PDBs are named `<sequence-id>#<variant>.pdb`, with zero-based variants, at the
top level of the output directory. IDs cannot contain path separators, `#`,
or control characters. The output directory must not already exist, including
an empty directory. Predictions are staged beside it and the complete directory
is published only after validating finite PDB coordinates and the expected
C-alpha count for every prediction. Failed inference identifies the input and
variant and leaves no published destination.

`fold-manifest.json` records the original FASTA hash, model source/revision,
package versions, GPU, effective per-layer dropout, precision, chunking,
fold seed, per-prediction seeds, timings and expected/completed counts.

Copy the complete directory to the Mac and run the base install in a separate
analysis directory:

```sh
uv run --project /path/to/structural-phylogeny --locked structphy \
  --structdir /path/to/reference-folds --n-bootstraps 10 \
  --seed 42 --threads 4 --outtree reference.newick
```

Existing PDB directories still work without a manifest. The existing combined
`structphy --fasta input.fasta` path uses the same backend and continues into
tree analysis. It still supports `--drop-inserts`, `--n-variants` and their
existing underscore aliases. Omitted dropout on that path retains evaluation
mode with zero attention dropout.

Options for the fold-only command must follow the `fold` subcommand; misplaced
top-level options are rejected. For the combined `--fasta --drop-inserts`
workflow, the sibling `<fold-dir>_conserved` directory must also be absent,
including an empty directory. This prevents stale predictions entering analysis.

## FASTA input

Sequence lines can wrap. StructPhy joins them and ignores whitespace within the
sequence. It uses the full header, with surrounding whitespace removed, as
the sequence ID. Empty records and duplicate headers fail before inference.
IDs that differ only by case or equivalent Unicode spelling are also rejected
before inference so the PDB names are distinct when copied to a Mac.

Existing alignment cleanup removes hyphens and uppercases residues before
folding. The combined FASTA-to-tree path retains its lowercase-insert handling
when `--drop-inserts` is supplied.

See [the implementation handoff](hpc-folding-handoff.md) for the measured
configuration, commands and remaining validation limits.
