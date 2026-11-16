# Dropout folding on an HPC

StructPhy should provide an optional `structphy[fold]` install for Linux GPU nodes. Run folding in the cluster's Python environment under its scheduler, then copy the PDB directory to the Mac for tree analysis. Keep the install recipe and the folding integration in this repo. The cluster already manages GPU allocation and drivers, so StructPhy can use the GPU allocated to the job directly. Retire the Docker SDK and the hard-coded Docker Hub image once the replacement works on a GPU node.

This plan assumes NVIDIA CUDA for the first tested configuration. Adapt the install and job commands to the actual cluster. The command lines below describe the proposed interface and require implementation.

## 1. Prove the folding backend on the cluster

1. On an allocated GPU node, record the GPU model, driver, `CUDA_VISIBLE_DEVICES`, and the Python available on the cluster. Create a scratch Python 3.12 environment with a driver-compatible PyTorch build and a pinned release of Transformers. Verify that PyTorch can allocate a CUDA tensor. Treat an install that succeeds on the login node and a CUDA run on the compute node as separate checks. [PyTorch installation](https://pytorch.org/get-started/locally/) · [uv with PyTorch](https://docs.astral.sh/uv/guides/integration/pytorch/)
2. Load `facebook/esmfold_v1` through `EsmForProteinFolding` at a fixed checkpoint revision. First fold one short reference sequence with dropout disabled and export its PDB. The Transformers port includes the OpenFold utilities used by its implementation. Start there before adding any of the legacy OpenFold or Conda dependencies. If compute nodes lack outbound access, stage the pinned snapshot from an allowed node and load it from shared scratch or a local checkpoint directory. [ESMFold documentation](https://huggingface.co/docs/transformers/en/model_doc/esm#transformers.EsmForProteinFolding) · [Folding implementation](https://github.com/huggingface/transformers/blob/main/src/transformers/models/esm/modeling_esmfold.py)
3. Inspect the **installed release** and confirm that ESM-2 attention dropout can be enabled selectively during inference. The old fork configured attention dropout per layer and called `train()` on the whole model. Start the new method with ESM-2 attention dropout only, with the rest of ESMFold in evaluation mode and gradients disabled. Check that the configured rate actually reaches the attention forward path. The published checkpoint sets attention dropout to zero, so a global mode switch alone gives no assurance that the intended mechanism is active. Keep any release-specific code in one small module. [Original dropout patch](https://github.com/FinnOD/esm-dropouts/commit/fefc647f25bebdecf862670ec4030a803105b298) · [Transformers ESM attention](https://github.com/huggingface/transformers/blob/main/src/transformers/models/esm/modeling_esm.py) · [Checkpoint configuration](https://huggingface.co/facebook/esmfold_v1/blob/main/config.json)
4. Compare several repeats of that sequence at zero dropout and at an explicit nonzero rate. Compare parsed coordinates within a tolerance, since PDB text differences alone are insufficient. Check that zero-dropout coordinates are stable on the same stack and that positive dropout can produce coordinate variation. Then try a sequence representative of the intended workload and record peak GPU memory. Test single-sequence inference before spending time on batching or precision tuning.

If the installed release makes selective dropout unreliable, try a compatible release and repeat the checks. If a model source patch is genuinely required, keep it small and pinned and record why the stock port was insufficient. A separately installed, pinned checkout of `esm-dropouts` can serve as a temporary reference or backend if the modern stack cannot run on the cluster. Adjust the dependency pins, model loading, batch size, and cluster setup as evidence comes in. Keep the FASTA-to-PDB contract below through those changes.

## 2. Add the local folding command

- Add a `fold` extra in [`pyproject.toml`](../pyproject.toml) for the tested folding dependencies and update `uv.lock`. Verify a fresh `uv sync --locked --extra fold` on the cluster. Also check that the built wheel advertises the extra and can install with `[fold]` in a clean Linux Python 3.12 environment using the documented PyTorch source. Checkout-specific uv settings do not travel with the wheel. If CUDA wheel selection needs cluster-specific configuration, document the tested index or environment setup in the README. Verify plain `uv sync --locked` on macOS with the folding dependencies absent. [uv with PyTorch](https://docs.astral.sh/uv/guides/integration/pytorch/)
- Add `structphy fold` through the existing Click entry point in [`structphy/cli.py`](../structphy/cli.py). Keep existing top-level `--structdir`, `--dmdir`, and `--fasta` invocations working. Keep the tree workflow's executable checks out of the fold-only callback. Import the heavyweight folding dependencies only when folding is requested and report missing dependencies or unavailable CUDA with an actionable error.
- Reuse the same folding function for `structphy fold` and the existing `structphy --fasta` convenience path. Preserve the latter's alignment cleanup and tree generation behavior, including `--drop-inserts`. The GPU job should be able to stop after folding. Once both entry points use the replacement, remove `run_inference_docker.py`, the Docker dependency and type-check exceptions, and the old Docker instructions.
- For `fold`, accept a FASTA, an output directory, `--n-variants` (default 10), an explicit `--dropout` rate, and `--fold-seed`. Use a scalar rate for the usual case and retain the existing comma-separated per-layer form for users with saved settings. Validate finite values in `[0, 1)` and the schedule length against the loaded model. Document that zero dropout is useful for a baseline. Keep the tree `--seed` for bootstrap sampling and keep the legacy `--fasta` behavior when dropout is omitted.

## 3. Keep the PDB directory portable

- Emit `<sequence-id>#<index>.pdb` with indices starting at zero. The existing tree code groups by the prefix before `#`; keep PDBs at the top level of the output directory so it finds them without changes. Reuse the repo's aligned-FASTA cleanup where applicable and handle wrapped records. Validate IDs before inference so duplicate IDs, path separators, or reserved `#` separators cannot create ambiguous or unsafe filenames. [Current FASTA utilities](../structphy/fasta_loading.py) · [Current bootstrap grouping](../structphy/generate_matrices.py)
- Load the model once per job. Start with straightforward serial inference and add batching only if measured throughput warrants it. Keep a predictable sequence and variant order. Record the effective seed and any batching or precision settings that affect reproducibility. Repeatability checks apply to the same software, device and execution settings.
- Stage outputs next to the target and publish the directory after verifying the expected number of nonempty PDBs per sequence. Refuse to overwrite existing results. On a failed prediction, exit nonzero and identify the sequence and variant. Keep incomplete work outside the final directory so the Mac tree workflow cannot accidentally consume it. An out-of-memory error should report the offending input and any relevant supported memory options. Resume support can wait until there is a demonstrated need.
- Include `fold-manifest.json` with the input FASTA hash, model name and revision, package versions, GPU and precision, effective dropout configuration, fold seed, and expected and completed counts. Keep the manifest beside the PDBs. Existing PDB directories remain valid input to `--structdir`.

## 4. Verify the two-machine workflow

From the checkout on the HPC, after selecting a dropout rate for the comparison and obtaining a GPU allocation:

```sh
uv sync --locked --extra fold
uv run --locked --extra fold structphy fold --fasta reference.fasta \
  --fold-dir runs/reference-folds --n-variants 10 \
  --dropout "$STRUCTPHY_DROPOUT_RATE" --fold-seed 42
```

Use several related reference sequences, ideally four or more distinct taxa, and copy the complete `reference-folds` directory to the Mac. Run the analysis from a clean run directory with the base install:

```sh
uv run --project /path/to/structural-phylogeny --locked structphy \
  --structdir /path/to/reference-folds --n-bootstraps 10 \
  --seed 42 --outtree reference.newick
```

- Confirm ten readable PDBs per input, a complete manifest, and a Newick tree with the expected taxa. Repeat the GPU smoke run with the same settings and compare parsed coordinates using a documented tolerance. If the old runtime is available, use the same FASTA to compare within-sequence structural spread and downstream tree support. Record the comparison as a method change because the old runner enabled training mode on the entire model. [Original inference runner](https://github.com/FinnOD/esm-dropouts/blob/95d537f19b97f99918ae05cc075edf3a3948a74b/scripts/esmfold_inference.py#L145-L165)
- Add focused CPU tests using a fake folding function for CLI dispatch, option validation, output naming and manifest completeness, and a failed prediction that never publishes an output directory. Check that the existing FASTA-to-tree path calls the shared folding function. Run the repo's lint, type, test, and packaging checks. Save the real GPU command, versions, observed results and unresolved issues in the implementation handoff. [Development checks](build-system.md#development-checks)
