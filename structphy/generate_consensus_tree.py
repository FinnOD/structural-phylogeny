from pathlib import Path
from typing import List
import tempfile
import subprocess
import os
from ete3 import Tree

def make_command_file(command_path: Path, tree_path: Path):
    with open(command_path, 'w') as f:
        f.write('\n')
        f.write(str(tree_path.resolve())+'\n')
        f.write('Y\n')
        f.write('\n')

def add_outgroup(tree_newick: str, outgroup_name: str) -> str:
    # Can't use ete3 because it doesn't place the outgroup as the first group
    # PHYLIP consense goes by leaf number, so here we can assure the outgroup as #1
    return f"({outgroup_name}:1.0, {tree_newick.replace(';', '')});"

def remove_outgroup(tree_newick: str, outgroup_name: str) -> str:
    tree = Tree(tree_newick)

    outgroup_node = tree.search_nodes(name=outgroup_name)[0]
    # We added the outgroup, so there should only be one sister to it
    assert len(outgroup_node.up.children) == 2
    left_root = outgroup_node.up.children[0]
    right_root = outgroup_node.up.children[1]
    
    # Return the root of the tree as the sister of the outgroup.
    # As it was before we added the outgroup.
    if left_root.name == outgroup_name:
        return right_root.write(format=0)
    if right_root.name == outgroup_name:
        return left_root.write(format=0)
    
def bootstrap_trees_to_consensus(bootstrap_trees: List[str], fake_outgroup: bool) -> str:
    CACHE_DIR = Path(os.environ["STRUCTPHY_CACHE_DIR"])

    outgroup_name = '!_OUTGROUP_!'
    if fake_outgroup:
        outgrouped_trees = [add_outgroup(newick, outgroup_name) for newick in bootstrap_trees]
        # Check to make sure that as text, the first group is the outgroup. 
        # This is required by consense as it uses numbered groups not the names.
        first_node = set(x.split(':')[0].replace('(', '').replace(')', '').replace(':', '') for x in outgrouped_trees)
        assert first_node == set([outgroup_name]), (f'{first_node} != {set([outgroup_name])}')
        bootstrap_trees = outgrouped_trees


    # Set up tempdir for running consense in
    with tempfile.TemporaryDirectory() as tmpdirname:
        temp_dir_path = Path(tmpdirname)
        command_path = temp_dir_path / 'command'
        intrees_path = temp_dir_path / 'intree'
        make_command_file(command_path, intrees_path)

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
    
    if fake_outgroup:
        consensus_tree = remove_outgroup(consensus_tree, outgroup_name)

    return consensus_tree