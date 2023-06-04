from typing import List
import pandas as pd
import subprocess
from tqdm.auto import tqdm
from multiprocessing import Pool


def fastme(phylip_matrix: str) -> str:
    # TODO This is a total hack, using stderr as an alternative pipe
    # Works fine if the command never fails :)
    command = ['./.structphy/fastme', '-i', '/dev/stdin', '-I', '/dev/null', '-O', '/dev/null', '-o', '/dev/stderr', '-s', '-m', 'NJ']
    result = subprocess.run(command, input=phylip_matrix.encode(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    newick = result.stderr.decode().strip()
    return newick

def convert_matrix_to_phylip(distance_df: pd.DataFrame) -> str: 

    split_bootstrap_id = lambda x: x.split('#')[0]
    distance_df.index = distance_df.index.map(split_bootstrap_id)
    distance_df.columns = distance_df.columns.map(split_bootstrap_id)

    phylip_out = ""
    phylip_out += (str(len(distance_df))+'\n')
    for i, row in distance_df.iterrows():
        phylip_out += (i+" "+" ".join(f'{x:.8f}' for x in list(row.values))+'\n')

    return phylip_out

# Wrapper required for multiprocessing
def matrix_to_fastme_newick(distance_df: pd.DataFrame) -> str:
    phylip_matrix = convert_matrix_to_phylip(distance_df)
    fastme_newick = fastme(phylip_matrix)
    return fastme_newick

def matrices_to_fastme_newick(distance_dfs: List[pd.DataFrame], n_threads: int) -> List[str]:
    
    fastme_trees = []
    with Pool(n_threads) as pool:
        for fastme_result in tqdm(
            pool.imap_unordered(matrix_to_fastme_newick, distance_dfs),
            total=len(distance_dfs),
            desc='Building trees',
            ascii=True,
            ):
            fastme_trees.append(fastme_result)
    
    return fastme_trees