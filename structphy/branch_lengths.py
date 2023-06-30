from typing import List
import pandas as pd
import numpy as np
from ete3 import Tree
from phylodm import PhyloDM
import dendropy

from structphy.bootstrapping import list_subtree_sets

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

  return tree.write(format=5)

def branch_lengths_on_tree(new_lengths, tree: Tree):
  tree = tree.copy()

  subtrees = sorted(list_subtree_sets(tree))
  subtree_new_length = dict(zip(subtrees, new_lengths, strict=True))

  for node in tree.traverse():
    node_set = frozenset(leaf.name.split(':')[0].strip() for leaf in node.get_leaves())
    branch_length = subtree_new_length[node_set]
    node.dist = branch_length 

  return tree

def get_branch_lengths(tree):
    subtrees = sorted(list_subtree_sets(tree))
    subtree_lengths = dict(zip(subtrees, [None for i in subtrees]))

    for node in tree.traverse():
        node_set = frozenset(leaf.name.split(':')[0].strip() for leaf in node.get_leaves())
        subtree_lengths[node_set] = node.dist

    return list(subtree_lengths.values())

def get_pdm_matrix(new_branch_lengths: List[float], base_tree: Tree):

  new_tree = branch_lengths_on_tree(new_branch_lengths, base_tree)
  dentropy_tree = dendropy.Tree.get(data=new_tree.write(format = 0), schema = 'newick')
  pdm = PhyloDM.load_from_dendropy(dentropy_tree)
  dm = pdm.dm(norm=False)

  labels = pdm.taxa()
  labels = [l.replace(' ', '_') for l in labels]
  dm_df = pd.DataFrame(dm, index=labels, columns=labels)
  dm_df.sort_index(inplace=True)
  dm_df.sort_index(axis=1, inplace=True)

  return dm_df

def optimise_branch_lengths(base_tree_newick: str, mean_distance_matrix: np.array):

    base_tree = Tree(base_tree_newick)

    def residuals(new_branch_lengths): # consensus, mean
        new_distance_matrix = get_pdm_matrix(new_branch_lengths, base_tree)

        residual_matrix = (mean_distance_matrix - new_distance_matrix).to_numpy()
        
        lower_triangle = residual_matrix[np.tril_indices(residual_matrix.shape[0], k = -1)]
        abs_losses = np.abs(lower_triangle)

        return abs_losses
        # return smooth_l1(abs_losses).sum()

    def residuals_lmfit(params, **kwargs):
        new_branch_lengths = np.array(list(params.valuesdict().values()))
        return residuals(new_branch_lengths)
    
    from lmfit import minimize, Parameters, fit_report
    import lmfit
    import math

    init_branch_lengths = np.array(get_branch_lengths(base_tree))
    init_result = residuals(init_branch_lengths).sum()
    print(init_result)

    def scale_model(params, branch_lengths):
        lin = params['m']*branch_lengths + params['c']
        return lin
    
    def scale_fit(params):
        scaled = scale_model(params, init_branch_lengths)
        return residuals(scaled)
    
    scale_params = Parameters()
    scale_params.add(f'm', value=1)
    scale_params.add(f'c', value=0)

    scale_result = minimize(scale_fit, scale_params, method='leastsq', max_nfev=1000)
    print(fit_report(scale_result))

    scaled_branch_lengths = scale_model(scale_result.params, init_branch_lengths)
    print(residuals(scaled_branch_lengths).sum())

    def on_iter(params, iteration, resid, *args, **kwargs):
        if iteration % 100 == 0:
            print(f"Iteration: {iteration: 4d} | Residual: {np.sum(resid):.4f} | {-100*(init_result-resid.sum())/init_result:.4f}%")

    params = Parameters()
    for i, bl in enumerate(scaled_branch_lengths):
        params.add(f'branch_{i}', value=bl, min=0, vary=True)

    leastsq_result = minimize(residuals_lmfit, params, method='leastsq', xtol=1e-10, ftol=1e-10, max_nfev=1_000, iter_cb=on_iter)
    out_tree = branch_lengths_on_tree(np.array(list(leastsq_result.params.valuesdict().values())), base_tree)
    print(fit_report(leastsq_result))

    # Evaluate result from reconstructed distance matrix vs. given mean distance matrix

    def tree_to_matrix(tree):
        dentropy_tree = dendropy.Tree.get(data=tree.write(format = 0), schema = 'newick')
        pdm = PhyloDM.load_from_dendropy(dentropy_tree)
        dm = pdm.dm(norm=False)

        labels = pdm.taxa()
        labels = [l.replace(' ', '_') for l in labels]
        dm_df = pd.DataFrame(dm, index=labels, columns=labels)
        dm_df.sort_index(inplace=True)
        dm_df.sort_index(axis=1, inplace=True)

        return dm_df
    
    mat = tree_to_matrix(out_tree).to_numpy()
    lower_tri_out = mat[np.tril_indices(mat.shape[0], k = -1)]
    mean_out = mean_distance_matrix[np.tril_indices(mean_distance_matrix.shape[0], k = -1)]

    rows = []
    for m, go in zip (mean_out, lower_tri_out):
        rows.append([m, go])
    out_df = pd.DataFrame(rows, columns=['mean target', 'final val'])
    out_df['res'] = out_df['mean target'] - out_df['final val']

    print(f'l1 loss: {out_df["res"].abs().sum():.4f}')
    print(f'l2 loss: {(out_df["res"]*out_df["res"]).sum():.4f}')
    print(f'MAE: {out_df["res"].abs().mean():.6f}')
    print(f'RMSE: {(out_df["res"]*out_df["res"]).mean()**0.5:.6f}')
    print(f'MAPE: {100*(out_df["res"]/out_df["mean target"]).mean():.4f}%')

    return out_tree.write(format=0)