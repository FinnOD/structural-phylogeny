from pathlib import Path
from typing import List
import tempfile
import subprocess
import os

def make_command_file(command_path: Path, tree_path: Path):
    with open(command_path, 'w') as f:
        f.write('\n')
        f.write(str(tree_path.resolve())+'\n')
        f.write('Y\n')
        f.write('\n')


def bootstrap_trees_to_consensus(bootstrap_trees: List[str]) -> str:
    CACHE_DIR = Path(os.environ["STRUCTPHY_CACHE_DIR"])

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
            consensus_tree = f.read()
        
    return consensus_tree