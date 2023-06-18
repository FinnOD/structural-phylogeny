from ete3 import Tree
from tqdm.auto import tqdm

def list_subtrees(node):
    subtrees = []
    if node.is_leaf():
        subtrees.append(node)
    else:
        for child in node.children:
            subtrees.extend(list_subtrees(child))
        subtrees.append(node)

    return subtrees
    
def list_subtree_sets(tree):
  sets = []
  s_t = list_subtrees(tree)
  for s in s_t:
    sets.append(frozenset(leaf.name.split(':')[0].strip() for leaf in s.get_leaves()))
  
  return sets

def bootstrap_against_tree(bootstrap_trees_newick, base_tree_newick):

    base_tree = Tree(base_tree_newick)

    base_subtrees = list_subtree_sets(base_tree)
    bootstrap_counts = {cluster:[] for cluster in base_subtrees}

    for bootstrap_newick in tqdm(bootstrap_trees_newick, desc='Applying bootstraps', ascii=True, position=0):
        bootstrap_tree = Tree(bootstrap_newick)
        bs_tree_clusters = list_subtree_sets(bootstrap_tree)

        for base_cluster in base_subtrees:
            if base_cluster in bs_tree_clusters:
                bootstrap_counts[base_cluster].append(1)
            else:
                bootstrap_counts[base_cluster].append(0)

    bootstrap_ratios = {cluster:sum(hits)/len(hits) for cluster, hits in bootstrap_counts.items()}

    out_tree = base_tree.copy()

    for node in out_tree.traverse():
        if node.is_leaf():
            parent = node.up
            node_set = frozenset(leaf.name.split(':')[0].strip() for leaf in parent.get_leaves())
        else:
            node_set = frozenset(leaf.name.split(':')[0].strip() for leaf in node.get_leaves())
        
        bs = f'{100 * bootstrap_ratios[node_set]:.3f}'
        node.add_features(support=bs)

    nhx_newick = out_tree.write(format=0, features=["support"])
    out_newick = nhx_newick.replace('&&NHX:support=', '')
    
    return out_newick
