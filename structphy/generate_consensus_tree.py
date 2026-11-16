import subprocess
import tempfile
from pathlib import Path

from ete3 import Tree

from structphy.executables import find_executable


def make_command_file(command_path: Path, tree_path: Path, outgroup_position: int | None):
    lines = ["", str(tree_path.resolve())]
    if outgroup_position:
        lines.extend(("O", str(outgroup_position + 1)))
    lines.extend(("Y", ""))
    command_path.write_text("\n".join(lines) + "\n")


def remove_outgroup(tree_newick: str, outgroup_name: str) -> str:
    tree = Tree(tree_newick)
    outgroup_node = tree.search_nodes(name=outgroup_name)[0]
    outgroup_node.delete()

    for node in tree.traverse():
        if len(node.children) == 1:
            return node.children[0].write(format=5)
    return tree.write(format=5)


def bootstrap_trees_to_consensus(bootstrap_trees: list[str], outgroup_name: str) -> str:
    executable = find_executable("consense")
    outgroup_position = None
    if outgroup_name:
        trees = [Tree(newick) for newick in bootstrap_trees]
        positions = {tree.get_leaf_names().index(outgroup_name) for tree in trees}
        if len(positions) != 1:
            message = "Bootstrap trees must have the same outgroup position."
            raise ValueError(message)
        outgroup_position = positions.pop()

    with tempfile.TemporaryDirectory() as directory:
        workdir = Path(directory)
        (workdir / "intree").write_text("\n".join(bootstrap_trees) + "\n")
        commands = ["O", str(outgroup_position + 1)] if outgroup_position else []
        commands.append("Y")
        subprocess.run(
            [executable], input="\n".join(commands) + "\n", cwd=workdir, check=True, capture_output=True, text=True
        )
        consensus_tree = (workdir / "outtree").read_text().replace("\n", "").strip()

    if outgroup_name:
        consensus_tree = remove_outgroup(consensus_tree, outgroup_name)
    return consensus_tree
