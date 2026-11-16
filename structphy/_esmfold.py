"""Transformers 4.57.2 ESMFold adapter, including selective ESM-2 attention dropout."""

import os
import platform
import sys
from importlib.metadata import version
from pathlib import Path
from typing import Any

import torch
import torch.version
import transformers
from transformers import EsmConfig, EsmForProteinFolding

from structphy.folding import MODEL_NAME, MODEL_REVISION, FoldError, expand_dropout


class EsmFoldBackend:
    def __init__(
        self, dropout: list[float] | None, model_dir: Path | None, chunk_size: int | None, threads: int
    ) -> None:
        if sys.platform != "linux" or not torch.cuda.is_available():
            message = (
                "ESMFold requires Linux and an available NVIDIA CUDA GPU. "
                "Run inside a GPU allocation with a driver-compatible PyTorch build; "
                "the base install can analyse existing PDBs on macOS."
            )
            raise FoldError(message)
        torch.set_num_threads(threads)
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        source = str(model_dir.resolve()) if model_dir is not None else MODEL_NAME
        load_options: dict[str, Any] = {"local_files_only": model_dir is not None}
        if model_dir is None:
            load_options["revision"] = MODEL_REVISION
        config = EsmConfig.from_pretrained(source, **load_options)
        if config.esmfold_config is None or not config.is_folding_model:
            message = "The checkpoint is not an ESMFold protein-folding model."
            raise FoldError(message)
        rates = expand_dropout(dropout, config.num_hidden_layers)
        self.model: Any = EsmForProteinFolding.from_pretrained(
            source, config=config, dtype=torch.float32, attn_implementation="eager", **load_options
        )
        self.model.eval().requires_grad_(False)
        for layer, rate in zip(self.model.esm.encoder.layer, rates, strict=True):
            attention = layer.attention.self
            # This pinned release uses functional dropout and the attention module's
            # own training flag. Enabling nn.Dropout children does not activate it.
            attention.dropout = rate
            attention.training = rate > 0
        if chunk_size is not None:
            self.model.trunk.set_chunk_size(chunk_size)
        try:
            self.model.cuda()
        except torch.cuda.OutOfMemoryError as error:
            message = "CUDA ran out of memory while loading ESMFold. Use a GPU with more free memory."
            raise FoldError(message) from error
        torch.cuda.reset_peak_memory_stats()
        device = torch.cuda.get_device_properties(torch.cuda.current_device())
        self.metadata: dict[str, Any] = {
            "model": MODEL_NAME,
            # A supplied directory may contain arbitrary weights; do not claim
            # that its revision has been independently verified.
            "model_revision": MODEL_REVISION if model_dir is None else None,
            "model_directory": str(model_dir.resolve()) if model_dir is not None else None,
            "dropout": {"scope": "esm2_self_attention_only", "rates": rates},
            "versions": {
                "structphy": version("structphy"),
                "python": platform.python_version(),
                "torch": torch.__version__,
                "transformers": transformers.__version__,
                "cuda_runtime": torch.version.cuda,
            },
            "device": {
                "name": device.name,
                "total_memory_bytes": device.total_memory,
                "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
                "logical_device": torch.cuda.current_device(),
            },
            "execution": {
                "precision": "float32",
                "attention_backend": "eager",
                "tf32": False,
                "threads": threads,
                "batch_size": 1,
                "chunk_size": chunk_size,
                "max_recycles": config.esmfold_config.trunk.max_recycles,
            },
        }

    def predict(self, sequence: str, seed: int) -> str:
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        try:
            with torch.inference_mode():
                pdb = self.model.infer_pdb(sequence)
        except torch.cuda.OutOfMemoryError as error:
            message = (
                "CUDA ran out of memory. Try --chunk-size 64 (or a smaller positive value) "
                "and check available GPU memory."
            )
            raise FoldError(message) from error
        else:
            self.metadata["execution"]["peak_gpu_memory_allocated_bytes"] = torch.cuda.max_memory_allocated()
            self.metadata["execution"]["peak_gpu_memory_reserved_bytes"] = torch.cuda.max_memory_reserved()
            return pdb
