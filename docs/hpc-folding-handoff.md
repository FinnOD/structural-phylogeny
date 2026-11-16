# HPC folding implementation handoff

Validated on 16 November 2026. Folding now runs in this repository's uv
environment through `structphy[fold]`. The Docker runner and Docker dependency
have been removed. [Usage and installation](folding.md).

## Implementation

- `structphy fold` generates PDBs and stops; the existing top-level
  `--fasta`, `--structdir` and `--dmdir` workflows remain available.
- Both FASTA entry points call `structphy.folding.fold_fasta`.
- `structphy/_esmfold.py` contains the release-specific Transformers adapter.
  Only ESM-2 self-attention dropout is active; gradients are disabled and all
  other modules remain in evaluation mode.
- The model loads once, predicts serially, and writes
  `<id>#<zero-based-variant>.pdb`. A sibling temporary directory holds incomplete
  predictions. The final destination must not already exist.
- Each prediction is checked for finite atom coordinates and the expected
  C-alpha count before publication. The manifest records provenance, settings,
  effective seeds, output counts, timings and observed peak CUDA memory.
- Scalar and 36-layer dropout schedules are supported, along with a local
  checkpoint directory and folding-trunk chunk size.
- The original HPC validation used single-line sequences and unique full headers.
  See the [folding guide](folding.md#fasta-input) for current FASTA behavior.

Omitting dropout on the combined FASTA-to-tree path retains zero attention
dropout and evaluation mode. Supplying dropout now affects ESM-2 attention
only. The old fork enabled whole-model training mode when dropout was supplied,
so this migration is a method change rather than a promise of numerical
equivalence to the old runner.

## Tested environment

| Component | Observed value |
| --- | --- |
| Host | virga-i1, Linux x86-64, SLE HPC 15 SP5, glibc 2.31 |
| Hardware | Four H100s, about 94 GiB each; 2 TiB host RAM; 72 CPU cores |
| Driver | 580.126.20 |
| Project interpreter | Python 3.12.15 in the checkout's .venv |
| uv | 0.13.4 |
| PyTorch | 2.7.1+cu126, CUDA runtime 12.6 |
| Transformers | 4.57.2 |
| Model | facebook/esmfold_v1 |
| Checkpoint revision | 75a3841ee059df2bf4d56688166c8fb459ddd97a |
| Inference | float32, eager attention, TF32 disabled, batch size one, four CPU threads |

No PyTorch or Transformers cluster modules were loaded for these tests.
Installation used `uv sync --locked --extra fold` and the explicit PyTorch
CUDA 12.6 index configured in `pyproject.toml`.

The GPU tests were direct, user-authorized interactive runs on idle physical
GPU 3, isolated with `CUDA_VISIBLE_DEVICES=3`. They were not Slurm jobs.
Do not copy that GPU override into scheduler jobs; preserve the allocation's
device visibility.

## Real CLI validation

Inputs were AURKA, ARAF, BRAF and MET from
`structphy/example_data/kinase_domains.fasta`, with lengths 251, 260, 260 and
264 residues. The selected records and outputs are saved under the Git-ignored
`runs/hpc-folding-validation/`.

The pinned checkpoint had already been downloaded into a temporary cache:

```sh
export HF_HUB_CACHE=/tmp/structphy-fold-smoke-20261116-ThknwB/hf-cache
export HF_HUB_OFFLINE=1
export CUDA_VISIBLE_DEVICES=3

uv run --locked --extra fold structphy fold \
  --fasta runs/hpc-folding-validation/reference.fasta \
  --fold-dir runs/hpc-folding-validation/dropout-42 \
  --n-variants 2 --dropout 0.1 --fold-seed 42 --threads 4
```

The same command was repeated with destination `dropout-42-repeat`.
A third run used `--dropout 0` and destination `zero-42`.
All three directories contain eight readable PDBs and complete manifests.
For future runs, stage a persistent shared cache using the folding guide;
the temporary cache above may be removed by the host.

Parsed coordinates matched exactly across the two separate seed-42 jobs.
Zero-dropout variants also matched exactly despite different derived prediction
seeds. The documented comparison tolerance was 0.001 angstrom. With dropout
0.1, the two variants showed these aligned C-alpha RMSDs:

| Taxon | RMSD (angstrom) |
| --- | ---: |
| AURKA | 0.281 |
| ARAF | 1.068 |
| BRAF | 0.418 |
| MET | 4.397 |

The fold seed is expanded into different stable seeds for each ID/variant;
it is not reused unchanged for every prediction. The observed structural spread
is a mechanism check, not a recommendation of dropout 0.1 as a scientific default.

The repeated run took 1.77–2.33 seconds per prediction after model loading,
with 14.19 GiB peak allocated GPU memory and 14.93 GiB peak reserved memory.
Detailed settings, timings and comparisons are in
`runs/hpc-folding-validation/validation.json` and the fold manifests.

From the validation directory, the generated PDBs were analysed with:

```sh
uv run --project /home/odo065/projects/structural-phylogeny \
  --locked --extra fold structphy \
  --structdir /home/odo065/projects/structural-phylogeny/runs/hpc-folding-validation/dropout-42 \
  --n-bootstraps 10 --threads 4 --seed 42 --outtree reference.newick
```

The Newick output contains exactly ARAF, AURKA, BRAF and MET.

Additional real checks passed:

- The combined `--fasta` path folded the corresponding aligned records from
  `kinase_domains_inserts.fasta`, removed lowercase inserts with
  `--drop-inserts`, and produced a four-taxon tree. It used one variant,
  omitted dropout, fold seed 42 and two tree bootstraps.
- A single AURKA prediction used a local checkpoint, `--chunk-size 64`,
  and 36 alternating rates of zero and 0.1. The manifest recorded local
  provenance with a null revision, the effective schedule and chunk size.
- With `CUDA_VISIBLE_DEVICES=''`, the fold CLI exited nonzero with an
  allocation/PyTorch diagnostic and created no output directory.

## Development and packaging checks

These commands passed in the project's fold-enabled environment:

```sh
uv run --locked --extra fold ruff check .
uv run --locked --extra fold ruff format --check .
uv run --locked --extra fold ty check
uv run --locked --extra fold pytest -q
uv build --no-sources
```

The suite has 45 passing tests. CPU tests cover CLI dispatch, the shared legacy
path and insert cleanup, option validation, filename safety, deterministic seed
derivation, manifest completeness, existing destinations, and failed or malformed
predictions that never publish a destination. Review regression tests also
reject misplaced parent options and existing conserved-output directories before
inference, preventing silently ignored seeds and stale structure reuse.

The packaging test builds from the source distribution, checks that wheel
metadata advertises `fold` with PyTorch and Transformers, excludes Docker,
installs the base wheel into a clean environment, and runs all bundled native
tools. It verifies that the base install contains neither folding dependency
and that attempting inference reports how to install `structphy[fold]`.

A clean temporary Python 3.12 environment also successfully installed the
built wheel with the extra using the documented command:

```sh
uv pip install --python /path/to/clean/venv/bin/python \
  'dist/structphy-0.0.1-py3-none-linux_x86_64.whl[fold]' --torch-backend=cu126
```

ESMFold imports succeeded there with PyTorch 2.7.1+cu126 and Transformers
4.57.2. This check uses wheel metadata rather than checkout uv source settings.
The PyTorch-specific backend option keeps ordinary dependencies on PyPI.

A macOS arm64 base-install dry run resolved successfully and excluded folding
dependencies. Native execution on a Mac was not available during this session.

## Known limits and remaining validation

- Slurm submission and the physical HPC-to-Mac transfer were not exercised.
- The GPU check used two variants per taxon, not a full ten-variant scientific
  comparison. Longer inputs, throughput tuning and larger studies remain untested.
- There is no numerical-equivalence comparison against the old fork, and tree
  support from this small smoke dataset does not validate the dropout method.
- Transformers emits a warning about newly initialized contact-head regression
  weights. The folding forward path does not call that contact-prediction head.
- The existing ETE 3 dependency emits Python deprecation/syntax warnings.
  Tests and tree generation completed successfully.
- Local checkpoint provenance is user-supplied; a revision is not asserted for
  arbitrary directories. Keep staged checkpoints pinned and documented.
