"""Shared FASTA-to-PDB workflow; heavy dependencies are loaded only for inference."""

import hashlib
import json
import math
import secrets
import time
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Protocol
from unicodedata import normalize

from structphy.fasta_loading import parse_fasta

MODEL_NAME = "facebook/esmfold_v1"
MODEL_REVISION = "75a3841ee059df2bf4d56688166c8fb459ddd97a"


class FoldError(RuntimeError):
    """A folding request could not produce a complete output directory."""


class FoldingBackend(Protocol):
    metadata: dict[str, Any]

    def predict(self, sequence: str, seed: int) -> str: ...


def parse_dropout(value: str) -> list[float]:
    """Parse a scalar or comma-separated ESM-2 attention dropout schedule."""
    try:
        rates = [float(part) for part in value.split(",")]
    except ValueError as error:
        message = "Expected a number or comma-separated numbers for dropout."
        raise FoldError(message) from error
    validate_dropout(rates)
    return rates


def validate_dropout(rates: list[float]) -> None:
    if not rates or any(not math.isfinite(rate) or not 0 <= rate < 1 for rate in rates):
        message = "Dropout rates must be finite numbers in [0, 1)."
        raise FoldError(message)


def expand_dropout(rates: list[float] | None, n_layers: int) -> list[float]:
    """Omitted dropout uses the published checkpoint's zero attention dropout."""
    if rates is None:
        return [0.0] * n_layers
    validate_dropout(rates)
    if len(rates) == 1:
        return rates * n_layers
    if len(rates) != n_layers:
        message = f"Expected one dropout rate or {n_layers} per-layer rates; received {len(rates)}."
        raise FoldError(message)
    return list(rates)


def _load_backend(
    dropout: list[float] | None, model_dir: Path | None, chunk_size: int | None, threads: int
) -> FoldingBackend:
    try:
        from structphy._esmfold import EsmFoldBackend
    except ImportError as error:
        message = (
            "Cannot import folding dependencies. Run `uv sync --locked --extra fold` "
            "or install `structphy[fold]` in a clean Python 3.12 environment. "
            f"Import error: {error}"
        )
        raise FoldError(message) from error
    try:
        return EsmFoldBackend(dropout, model_dir, chunk_size, threads)
    except FoldError:
        raise
    except Exception as error:
        message = (
            f"Could not load ESMFold: {error}. Check the checkpoint/cache and outbound access; "
            "use --model-dir for an already staged checkpoint."
        )
        raise FoldError(message) from error


def _validate_ids(sequences: dict[str, str]) -> None:
    filename_ids: dict[str, str] = {}
    for name, sequence in sequences.items():
        if (
            not name.strip()
            or name in {".", ".."}
            or any(char in name for char in "/\\#")
            or any(ord(char) < 32 or ord(char) == 127 for char in name)
        ):
            message = (
                f"Sequence ID {name!r} is not safe for a PDB filename; avoid path separators, # and control characters."
            )
            raise FoldError(message)
        if not sequence:
            message = f"Sequence {name!r} is empty after removing alignment gaps."
            raise FoldError(message)
        filename_key = normalize("NFD", normalize("NFD", name).casefold())
        if filename_key in filename_ids:
            other = filename_ids[filename_key]
            message = f"Sequence IDs {other!r} and {name!r} would collide as PDB filenames on a Mac."
            raise FoldError(message)
        filename_ids[filename_key] = name


def _validate_pdb(pdb: str, n_residues: int) -> int:
    atoms = 0
    ca_atoms = 0
    for line in pdb.splitlines():
        if not line.startswith("ATOM  "):
            continue
        try:
            xyz = [float(line[start : start + 8]) for start in (30, 38, 46)]
        except ValueError as error:
            message = "The prediction contains unreadable PDB coordinates."
            raise FoldError(message) from error
        if not all(math.isfinite(value) for value in xyz):
            message = "The prediction contains nonfinite PDB coordinates."
            raise FoldError(message)
        atoms += 1
        ca_atoms += line[12:16].strip() == "CA"
    if ca_atoms != n_residues:
        message = f"Expected {n_residues} C-alpha atoms in the prediction; found {ca_atoms}."
        raise FoldError(message)
    return atoms


def prediction_seed(fold_seed: int, sequence_id: str, variant: int) -> int:
    """Keep each prediction's RNG independent of input order and previous predictions."""
    digest = hashlib.sha256(f"{fold_seed}\0{sequence_id}\0{variant}".encode()).digest()
    return int.from_bytes(digest[:8], "big") % (2**63)


def fold_fasta(
    fasta: Path,
    fold_dir: Path,
    *,
    n_variants: int = 10,
    dropout: list[float] | None = None,
    fold_seed: int | None = None,
    model_dir: Path | None = None,
    chunk_size: int | None = None,
    threads: int = 4,
) -> Path:
    """Publish a complete PDB directory, or leave the destination absent on failure."""
    if n_variants < 1 or threads < 1 or (chunk_size is not None and chunk_size < 1):
        message = "Variant count, thread count and chunk size must be positive."
        raise FoldError(message)
    if fold_seed is not None and not 0 <= fold_seed < 2**63:
        message = "The folding seed must be in [0, 2**63)."
        raise FoldError(message)
    if dropout is not None:
        validate_dropout(dropout)
    fold_dir = fold_dir.absolute()
    if fold_dir.exists() or fold_dir.is_symlink():
        message = f"Refusing to overwrite an existing folding directory: {fold_dir}"
        raise FoldError(message)
    try:
        input_bytes = fasta.read_bytes()
        _, sequences = parse_fasta(input_bytes.decode("utf-8-sig"))
    except (OSError, ValueError) as error:
        message = f"Could not read FASTA {fasta}: {error}"
        raise FoldError(message) from error
    _validate_ids(sequences)
    backend = _load_backend(dropout, model_dir, chunk_size, threads)
    effective_seed = secrets.randbits(63) if fold_seed is None else fold_seed
    expected = dict.fromkeys(sorted(sequences), n_variants)
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "input_fasta": {"path": str(fasta.resolve()), "sha256": hashlib.sha256(input_bytes).hexdigest()},
        **backend.metadata,
        "fold_seed": effective_seed,
        "seed_scheme": "sha256(fold_seed\\0sequence_id\\0variant), first 8 bytes modulo 2**63",
        "n_variants": n_variants,
        "sequence_order": sorted(sequences),
        "expected_counts": expected,
        "completed_counts": dict.fromkeys(expected, 0),
        "predictions": [],
    }
    fold_dir.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=f".{fold_dir.name}.", dir=fold_dir.parent) as temporary:
        staging = Path(temporary)
        for sequence_id in sorted(sequences):
            sequence = sequences[sequence_id]
            for variant in range(n_variants):
                seed = prediction_seed(effective_seed, sequence_id, variant)
                filename = f"{sequence_id}#{variant}.pdb"
                started = time.monotonic()
                try:
                    pdb = backend.predict(sequence, seed)
                    atoms = _validate_pdb(pdb, len(sequence))
                    (staging / filename).write_text(pdb)
                except Exception as error:
                    message = f"Folding failed for sequence {sequence_id!r}, variant {variant}: {error}"
                    raise FoldError(message) from error
                manifest["completed_counts"][sequence_id] += 1
                manifest["predictions"].append(
                    {
                        "sequence_id": sequence_id,
                        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                        "residues": len(sequence),
                        "variant": variant,
                        "seed": seed,
                        "filename": filename,
                        "atoms": atoms,
                        "seconds": time.monotonic() - started,
                    }
                )
        (staging / "fold-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
        if fold_dir.exists() or fold_dir.is_symlink():
            message = f"Folding destination appeared during inference; refusing to overwrite: {fold_dir}"
            raise FoldError(message)
        staging.rename(fold_dir)
    return fold_dir
