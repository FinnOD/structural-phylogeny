import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from click.testing import CliRunner

from structphy import cli, folding


def pdb_for(sequence: str) -> str:
    return (
        "".join(
            f"ATOM  {index:5d}  CA  ALA A{index:4d}    {float(index):8.3f}{0.0:8.3f}{0.0:8.3f}"
            "  1.00 90.00           C\n"
            for index in range(1, len(sequence) + 1)
        )
        + "END\n"
    )


class FakeBackend:
    metadata = {
        "model": folding.MODEL_NAME,
        "model_revision": folding.MODEL_REVISION,
        "dropout": {"scope": "esm2_self_attention_only", "rates": [0.1] * 36},
        "versions": {"torch": "fake", "transformers": "fake"},
        "device": {"name": "fake"},
        "execution": {"precision": "float32", "batch_size": 1},
    }

    def __init__(self):
        self.calls = []

    def predict(self, sequence: str, seed: int) -> str:
        self.calls.append((sequence, seed))
        return pdb_for(sequence)


@pytest.fixture
def fake_backend(monkeypatch):
    backend = FakeBackend()
    loads = []

    def load(*args):
        loads.append(args)
        return backend

    monkeypatch.setattr(folding, "_load_backend", load)
    return backend, loads


@pytest.fixture
def fasta(tmp_path: Path) -> Path:
    path = tmp_path / "input.fasta"
    path.write_text(">B\nAA\n>A\nA-a\n")
    return path


@pytest.mark.parametrize("value", ["", "nan", "inf", "-0.1", "1", "0.2,,0.1", "letters"])
def test_invalid_dropout(value):
    with pytest.raises(folding.FoldError):
        folding.parse_dropout(value)


def test_dropout_schedule():
    assert folding.expand_dropout(None, 36) == [0.0] * 36
    assert folding.expand_dropout(folding.parse_dropout("0.1"), 36) == [0.1] * 36
    rates = [index / 100 for index in range(36)]
    assert folding.expand_dropout(rates, 36) == rates
    with pytest.raises(folding.FoldError, match="36 per-layer"):
        folding.expand_dropout([0.1, 0.2], 36)


def test_complete_output_and_manifest(fasta, tmp_path, fake_backend):
    backend, loads = fake_backend
    destination = tmp_path / "folds"
    result = folding.fold_fasta(fasta, destination, n_variants=2, dropout=[0.1], fold_seed=42)

    assert result == destination
    assert {p.name for p in destination.glob("*.pdb")} == {"A#0.pdb", "A#1.pdb", "B#0.pdb", "B#1.pdb"}
    assert len(loads) == 1
    assert loads[0] == ([0.1], None, None, 4)
    manifest = json.loads((destination / "fold-manifest.json").read_text())
    assert manifest["input_fasta"]["sha256"] == hashlib.sha256(fasta.read_bytes()).hexdigest()
    assert manifest["expected_counts"] == manifest["completed_counts"] == {"A": 2, "B": 2}
    assert manifest["sequence_order"] == ["A", "B"]
    assert manifest["fold_seed"] == 42
    assert manifest["model_revision"] == folding.MODEL_REVISION
    assert len(manifest["predictions"]) == 4
    assert all(p["residues"] == 2 for p in manifest["predictions"])
    assert backend.calls == [("AA", p["seed"]) for p in manifest["predictions"]]
    assert len({seed for _, seed in backend.calls}) == 4
    assert not list(tmp_path.glob(".folds.*"))


def test_folding_seed_does_not_depend_on_input_order(fasta, tmp_path, fake_backend):
    backend, _ = fake_backend
    folding.fold_fasta(fasta, tmp_path / "first", n_variants=2, fold_seed=42)
    first = list(backend.calls)
    backend.calls.clear()
    fasta.write_text(">A\nA-a\n>B\nAA\n")
    folding.fold_fasta(fasta, tmp_path / "second", n_variants=2, fold_seed=42)
    assert backend.calls == first


@pytest.mark.usefixtures("fake_backend")
def test_omitted_seed_is_recorded(fasta, tmp_path):
    folding.fold_fasta(fasta, tmp_path / "folds", n_variants=1)
    manifest = json.loads((tmp_path / "folds/fold-manifest.json").read_text())
    assert 0 <= manifest["fold_seed"] < 2**63
    for prediction in manifest["predictions"]:
        assert prediction["seed"] == folding.prediction_seed(
            manifest["fold_seed"], prediction["sequence_id"], prediction["variant"]
        )


def test_failed_prediction_never_publishes(fasta, tmp_path, fake_backend, monkeypatch):
    backend, _ = fake_backend

    def predict(sequence, _seed):
        if len(backend.calls) == 1:
            message = "GPU failed"
            raise RuntimeError(message)
        backend.calls.append(sequence)
        return pdb_for(sequence)

    monkeypatch.setattr(backend, "predict", predict)
    with pytest.raises(folding.FoldError, match="sequence 'A', variant 1: GPU failed"):
        folding.fold_fasta(fasta, tmp_path / "folds", n_variants=2)
    assert not (tmp_path / "folds").exists()
    assert not list(tmp_path.glob(".folds.*"))


@pytest.mark.parametrize("pdb", ["", "not a PDB", pdb_for("A"), pdb_for("AA").replace("   1.000", "     nan")])
def test_invalid_prediction_never_publishes(pdb, fasta, tmp_path, fake_backend, monkeypatch):
    backend, _ = fake_backend
    monkeypatch.setattr(backend, "predict", lambda *_args: pdb)
    with pytest.raises(folding.FoldError, match="Folding failed"):
        folding.fold_fasta(fasta, tmp_path / "folds", n_variants=1)
    assert not (tmp_path / "folds").exists()


def test_existing_directory_is_never_overwritten(fasta, tmp_path, fake_backend):
    _, loads = fake_backend
    destination = tmp_path / "folds"
    destination.mkdir()
    with pytest.raises(folding.FoldError, match="Refusing to overwrite"):
        folding.fold_fasta(fasta, destination)
    assert not loads


@pytest.mark.parametrize("name", ["A#0", "../A", "A/B", "A\\B", ".."])
def test_unsafe_filename_rejected_before_loading(name, fasta, tmp_path, fake_backend):
    _, loads = fake_backend
    fasta.write_text(f">{name}\nAA\n")
    with pytest.raises(folding.FoldError, match="not safe"):
        folding.fold_fasta(fasta, tmp_path / "folds")
    assert not loads


@pytest.mark.usefixtures("fake_backend")
def test_fold_subcommand_bypasses_tree_tools(fasta, tmp_path, monkeypatch):
    def unexpected_tool(_name):
        pytest.fail("Fold-only invocation must not check tree executables")

    monkeypatch.setattr(cli, "find_executable", unexpected_tool)
    result = CliRunner().invoke(
        cli.main,
        [
            "fold",
            "--fasta",
            str(fasta),
            "--fold-dir",
            str(tmp_path / "folds"),
            "--n-variants",
            "2",
            "--dropout",
            "0.1",
            "--fold-seed",
            "42",
        ],
    )
    assert result.exit_code == 0, result.output
    assert (tmp_path / "folds/fold-manifest.json").is_file()


@pytest.mark.parametrize(
    "options",
    [
        [],
        ["--dropout", "nan"],
        ["--dropout", "1"],
        ["--dropout", "0", "--n-variants", "0"],
        ["--dropout", "0", "--fold-seed", "-1"],
        ["--dropout", "0", "--chunk-size", "0"],
    ],
)
def test_fold_cli_rejects_invalid_options(options, fasta, tmp_path, fake_backend):
    _, loads = fake_backend
    result = CliRunner().invoke(
        cli.main, ["fold", "--fasta", str(fasta), "--fold-dir", str(tmp_path / "folds"), *options]
    )
    assert result.exit_code == 2
    assert not loads


def test_cli_reports_backend_failure(fasta, tmp_path, monkeypatch):
    def unavailable(*_args):
        message = "CUDA unavailable; allocate a GPU"
        raise folding.FoldError(message)

    monkeypatch.setattr(folding, "_load_backend", unavailable)
    result = CliRunner().invoke(
        cli.main, ["fold", "--fasta", str(fasta), "--fold-dir", str(tmp_path / "folds"), "--dropout", "0"]
    )
    assert result.exit_code == 1
    assert "CUDA unavailable" in result.output
    assert not (tmp_path / "folds").exists()


def test_base_cli_does_not_import_folding_dependencies():
    subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; import structphy.cli; "
            "assert 'torch' not in sys.modules; assert 'transformers' not in sys.modules",
        ],
        check=True,
        timeout=30,
    )


@pytest.mark.parametrize("dropout", [None, "0.1"])
def test_legacy_fasta_workflow_uses_shared_folding(fasta, tmp_path, monkeypatch, dropout):
    from structphy import extract_conserved_pdb

    monkeypatch.chdir(tmp_path)
    calls = []
    conserved = []
    matrix_inputs = []

    def fake_fold(input_fasta, fold_dir, **kwargs):
        calls.append((input_fasta, kwargs))
        fold_dir.mkdir()
        for name in ("A", "B"):
            (fold_dir / f"{name}#0.pdb").write_text(pdb_for("AA"))
        return fold_dir

    def remove_inserts(source, alignment, output):
        conserved.append(alignment)
        output.write_text(source.read_text())

    def matrices(paths, **kwargs):
        matrix_inputs.extend(paths)
        assert kwargs["seed"] == 19
        return [pd.DataFrame([[0, 2], [2, 0]], index=["A", "B"], columns=["A", "B"])]

    monkeypatch.setattr(cli, "fold_fasta", fake_fold)
    monkeypatch.setattr(cli, "find_executable", lambda _name: None)
    monkeypatch.setattr(extract_conserved_pdb, "remove_inserts_from_structure", remove_inserts)
    monkeypatch.setattr(cli, "generate_bootstrap_matrices_from_structures", matrices)
    monkeypatch.setattr(cli, "matrices_to_fastme_newick", lambda *_args, **_kwargs: ["(A:1,B:1,!_OUTGROUP_!:1);"])
    monkeypatch.setattr(cli, "bootstrap_trees_to_consensus", lambda *_args, **_kwargs: "(A:1,B:1);")
    options = ["--dropout", dropout] if dropout is not None else []
    result = CliRunner().invoke(
        cli.main,
        [
            "--fasta",
            str(fasta),
            "--drop-inserts",
            "--fold-seed",
            "42",
            "--seed",
            "19",
            "--n-variants",
            "2",
            "--outtree",
            "result.newick",
            *options,
        ],
    )

    assert result.exit_code == 0, result.output
    assert len(calls) == 1
    assert calls[0][0] == fasta
    assert calls[0][1]["dropout"] == ([0.1] if dropout is not None else None)
    assert calls[0][1]["fold_seed"] == 42
    assert calls[0][1]["n_variants"] == 2
    assert sorted(conserved) == ["A-a", "AA"]
    assert all(path.parent.name == "input_folddir_conserved" for path in matrix_inputs)
    assert (tmp_path / "result.newick").is_file()


@pytest.mark.parametrize("parent_options", [["--threads", "1"], ["--fold-seed", "123"], ["--n-variants", "2"]])
def test_fold_rejects_misplaced_parent_options(parent_options, fasta, tmp_path, fake_backend):
    _, loads = fake_backend
    result = CliRunner().invoke(
        cli.main,
        [*parent_options, "fold", "--fasta", str(fasta), "--fold-dir", str(tmp_path / "folds"), "--dropout", "0.1"],
    )
    assert result.exit_code == 2
    assert "Put options after the fold subcommand" in result.output
    assert parent_options[0] in result.output
    assert not loads


@pytest.mark.parametrize("stale_file", [False, True])
def test_combined_workflow_refuses_existing_conserved_directory(stale_file, fasta, tmp_path, fake_backend, monkeypatch):
    _, loads = fake_backend
    monkeypatch.setattr(cli, "find_executable", lambda _name: None)
    conserved = tmp_path / "folds_conserved"
    conserved.mkdir()
    sentinel = conserved / "OLD#0.pdb"
    if stale_file:
        sentinel.write_text("previous prediction")
    result = CliRunner().invoke(
        cli.main, ["--fasta", str(fasta), "--fold-dir", str(tmp_path / "folds"), "--drop-inserts"]
    )
    assert result.exit_code == 2
    assert "Refusing to reuse an existing conserved directory" in result.output
    assert not loads
    assert not (tmp_path / "folds").exists()
    if stale_file:
        assert sentinel.read_text() == "previous prediction"
