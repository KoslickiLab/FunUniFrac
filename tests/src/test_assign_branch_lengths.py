import data
from temp.real_data_branch_assignment import KeggTree, \
    get_subtree, assign_branch_lengths, get_KeggTree_from_edgelist
import networkx as nx


edge_list = 'kegg_ko_edge_df_br_ko00001.txt_AAI_lengths_n_50_f_10_r_100.txt'
G = nx.read_edgelist(data.get_data_abspath(edge_list), delimiter='\t', nodetype=str, create_using=nx.DiGraph, data=(('edge_length', float),))
subgraph_nodes = ['K19768', 'K19773', 'K19765', 'K22878', 'K19767', 'K19764', 'K19774', 'K10141', 'K19766',
                  'K19805', '04212 Longevity regulating pathway - worm', 'K11204', 'K14938', 'K20394',
                  '04213 Longevity regulating pathway - multiple species', 'K19772', 'K13356', 'K17705', '09149 Aging',
                  'K01768', 'K19769', 'K19770', 'K19771', 'K01440', '04211 Longevity regulating pathway']
subtree = G.subgraph(subgraph_nodes)

def test_KeggTree():
    kegg_tree = KeggTree(subtree)
    assert len(kegg_tree.tree.nodes) == len(subtree.nodes)
    assert len(kegg_tree.tree.edges) == len(subtree.edges)
    assert len(kegg_tree.leaf_nodes) == 21
    assert len(kegg_tree.pw_dist) == 210


def test_get_siblings():
    kegg_tree = KeggTree(subtree)
    siblings = kegg_tree.get_siblings('K01768')
    assert len(siblings) == 7
    assert 'K01440' in siblings


def test_multiple_merge_branches():
    edge_list = data.get_data_abspath('test_input/tree_with_multiple_single_descendants.txt')
    kegg_tree = get_KeggTree_from_edgelist(edge_list)
    assert kegg_tree.tree.get_edge_data('A', 'C')['edge_length'] == 3
    assert 'B' not in kegg_tree.tree


def test_assign_branch_with_multiple_single_child():
    edge_list = data.get_data_abspath('test_input/unbalanced_tree.txt')
    kegg_tree = get_KeggTree_from_edgelist(edge_list)
    edge_lengths_solution = {}
    pw_dist = kegg_tree.pw_dist
    assign_branch_lengths(kegg_tree, kegg_tree.leaf_nodes, pw_dist, edge_lengths_solution)
    assert edge_lengths_solution[('root', 'f')] == 8
    assert edge_lengths_solution[('root', 'h')] == 1


def test_make_full_tree():
    edge_list = data.get_data_abspath('test_input/tree_with_multiple_single_descendants.txt')
    kegg_tree = get_KeggTree_from_edgelist(edge_list)
    kegg_tree.make_full_tree()
    assert kegg_tree.size == 7



