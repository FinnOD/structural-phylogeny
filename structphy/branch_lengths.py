import pandas as pd
import numpy as np
from ete3 import Tree

def get_stacked(matrices):
  
  idx = None
  vals = []
     
  for df in matrices:
    if idx is None:
      idx = df.stack().index
    vals.append(df.stack().values)

  return pd.DataFrame(
    np.concatenate([vals]).T, # this bs is faster than pandas native
    index=idx,
    columns=None
  )

def get_mean_distance_matrix(stacked):
  mean_df = pd.DataFrame(stacked.mean(axis=1)).reset_index()
  strip_bs_number = lambda id: id.split('#')[0]
  mean_df['level_0'] = mean_df['level_0'].apply(strip_bs_number)
  mean_df['level_1'] = mean_df['level_1'].apply(strip_bs_number)
  return pd.pivot_table(mean_df, index='level_0', columns='level_1', values=0)

def get_upgma_tree(newick_tree: str, distance_matrix):
  tree = Tree(newick_tree)

  for node in list(tree.traverse()):
    if not node.children:
        continue


    if len(node.children) == 2:
        left_leaves = [str(x.name) for x in node.children[0].iter_leaves()]
        right_leaves = [str(x.name) for x in node.children[1].iter_leaves()]

        sibling_distance = distance_matrix.loc[left_leaves][right_leaves].mean().mean()

        node.children[0].dist = sibling_distance/2
        node.children[1].dist = sibling_distance/2


    if len(node.children) == 3:

        left_leaves = [str(x.name) for x in node.children[0].iter_leaves()]
        middle_leaves = [str(x.name) for x in node.children[1].iter_leaves()]
        right_leaves = [str(x.name) for x in node.children[2].iter_leaves()]

        tern_dists = [ 
        distance_matrix.loc[left_leaves][middle_leaves],
        distance_matrix.loc[left_leaves][right_leaves],
        distance_matrix.loc[middle_leaves][right_leaves]
        ]

        sibling_distance = pd.concat(tern_dists).mean().mean()

        node.children[0].dist = sibling_distance/3
        node.children[1].dist = sibling_distance/3
        node.children[2].dist = sibling_distance/3

    return tree.write(format=0)