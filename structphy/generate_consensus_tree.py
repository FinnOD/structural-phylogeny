from pathlib import Path
from typing import List
import tempfile
import subprocess
import os
from ete3 import Tree

def make_command_file(command_path: Path, tree_path: Path, outgroup_position: int):
    with open(command_path, 'w') as f:
        f.write('\n')
        f.write(str(tree_path.resolve())+'\n')
        if outgroup_position:
            f.write('O\n')
            f.write(f'{outgroup_position+1}\n')
        f.write('Y\n')
        f.write('\n')

def remove_outgroup(tree_newick: str, outgroup_name: str) -> str:
    tree = Tree(tree_newick)

    outgroup_node = tree.search_nodes(name=outgroup_name)[0]
    outgroup_node.delete()

    for node in list(tree.traverse()):
        if len(node.children) == 1:
            return node.children[0].write(format=5)

    return tree.write(format=5)
    
def bootstrap_trees_to_consensus(bootstrap_trees: List[str], outgroup_name: str) -> str:
    CACHE_DIR = Path(os.environ["STRUCTPHY_CACHE_DIR"])

    outgroup_position = None
    if outgroup_name:
        outgrouped_trees = [Tree(newick) for newick in bootstrap_trees]
        leaves_all = [tree.get_leaf_names() for tree in outgrouped_trees]
        outgroup_positions = set(leaves_list.index(outgroup_name) for leaves_list in leaves_all)
        assert len(outgroup_positions) == 1
        outgroup_position = outgroup_positions.pop()
        print(outgroup_position)


    # Set up tempdir for running consense in
    with tempfile.TemporaryDirectory() as tmpdirname:
        temp_dir_path = Path(tmpdirname)
        command_path = temp_dir_path / 'command'
        intrees_path = temp_dir_path / 'intree'
        make_command_file(command_path, intrees_path, outgroup_position=outgroup_position)

        # Have to use a shell script to get around UI of consense
        with open(temp_dir_path / 'run.sh', 'w') as f:
            f.write(str(CACHE_DIR / 'consense') + ' < ' + str(command_path))

        with open(intrees_path, 'w') as f:
            f.writelines(bootstrap_trees)

        # Finally run consense 
        consense_process = subprocess.run(['sh', temp_dir_path / 'run.sh'], cwd=str(temp_dir_path), check=True, capture_output=True, text=True)
        consense_logs = consense_process.stdout + consense_process.stderr

        with open(temp_dir_path / 'outtree', 'r') as f:
            consensus_tree = f.read().replace('\n', '').strip()
    
    if outgroup_name:
        consensus_tree = remove_outgroup(consensus_tree, outgroup_name)

    return consensus_tree