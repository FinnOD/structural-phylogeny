from pathlib import Path

import pandas as pd
from click.testing import CliRunner

from structphy import cli, generate_matrices
from structphy.generate_trees import convert_matrix_to_phylip


def test_cli_requires_exactly_one_input(tmp_path: Path):
    runner = CliRunner()
    for args in ([], ["--structdir", str(tmp_path), "--dmdir", str(tmp_path)]):
        result = runner.invoke(cli.main, args)
        assert result.exit_code == 2
        assert "Provide one input" in result.output


def test_matrix_workflow_accepts_new_output_and_custom_base_tree(tmp_path: Path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    requested_tools = []
    monkeypatch.setattr(cli, "find_executable", requested_tools.append)
    monkeypatch.setattr(cli, "matrices_to_fastme_newick", lambda *_args, **_kwargs: ["(A:1,B:1,!_OUTGROUP_!:1);"])
    monkeypatch.setattr(cli, "bootstrap_trees_to_consensus", lambda *_args, **_kwargs: "(A:1,B:1);")
    pd.DataFrame([[0, 2], [2, 0]], index=["A", "B"], columns=["A", "B"]).to_csv("matrix.csv")
    Path("base.newick").write_text("(A:2,B:2);\n")

    result = CliRunner().invoke(cli.main, ["--dmdir", ".", "--base-tree", "base.newick", "-o", "result.newick"])

    assert result.exit_code == 0, result.output
    assert requested_tools == ["fastme", "consense"]
    assert "A:2" in Path("result.newick").read_text()
    assert "!_OUTGROUP_!" not in Path("result.newick").read_text()


def test_phylip_labels_and_outgroup_diagonal():
    labels = ["Human kinase#0", "Mouse kinase#1"]
    matrix = pd.DataFrame([[0, 0.5], [0.5, 0]], index=labels, columns=labels, dtype=float)
    generate_matrices.make_fake_outgroups([matrix], "outgroup")

    assert matrix.loc["outgroup", "outgroup"] == 0.0
    assert matrix.loc[labels[0], "outgroup"] == 1_000_000.0
    assert convert_matrix_to_phylip(matrix).splitlines()[1] == "Human_kinase 0.00000000 0.50000000 1000000.00000000"
    assert matrix.index[0] == "Human kinase#0"


def test_seeded_bootstraps_are_independent_of_file_order(monkeypatch):
    def sample_matrix(paths: list[Path], n_threads: int) -> pd.DataFrame:
        return pd.DataFrame({"selection": [path.name for path in paths], "threads": n_threads})

    monkeypatch.setattr(generate_matrices, "generate_matrix_from_bootstraps", sample_matrix)
    variants = [Path(f"{taxon}#{index}.pdb") for taxon in ("a", "b") for index in range(4)]
    first = generate_matrices.generate_bootstrap_matrices_from_structures(variants, 1, 8, seed=42)
    repeated = generate_matrices.generate_bootstrap_matrices_from_structures(variants[::-1], 1, 8, seed=42)

    for matrix, repeated_matrix in zip(first, repeated, strict=True):
        pd.testing.assert_frame_equal(matrix, repeated_matrix)
    assert len({tuple(matrix["selection"]) for matrix in first}) > 1
