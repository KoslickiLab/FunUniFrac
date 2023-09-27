import networkx as nx
from Kegg_tree import KeggTree, write_subgraph_file
import random
import argparse

#original tree
edge_list = '../data/kegg_ko_edge_df_br_ko00001.txt_AAI_lengths_n_50_f_10_r_100.txt'
G = nx.read_edgelist(edge_list, delimiter='\t', nodetype=str, create_using=nx.DiGraph, data=(('edge_length', float),))

def extract_subtree_by_level_random(G, level):
    kegg_tree = KeggTree(G)
    kegg_tree.group_nodes_by_depth()
    root = random.choice(kegg_tree.nodes_by_depth[level])
    subtree_nodes = nx.descendants(G, root)
    subtree_nodes.add(root)
    sub_tree = G.subgraph(subtree_nodes)
    return sub_tree

def extract_subtree_by_root(G, root):
    subtree_nodes = nx.descendants(G, root)
    subtree_nodes.add(root)
    sub_tree = G.subgraph(subtree_nodes)
    return sub_tree

def perturb_edges(subgraph):
    new_tree = nx.DiGraph()
    for edge in subgraph.edges.data():
        new_edge = edge[2]['edge_length'] + random.random()
        new_tree.add_edge(edge[0], edge[1], edge_length=new_edge)
    new_kegg_tree = KeggTree(new_tree)
    new_kegg_tree.group_nodes_by_depth()
    return new_tree


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--level', type=int, default=3)
    args = parser.parse_args()
    sub_tree = extract_subtree_by_level_random(G, args.level)
    new_tree = perturb_edges(sub_tree)
    outfile_name = 'subtree_size' + str(new_tree.size()) + '.txt'
    print(outfile_name)
    write_subgraph_file(new_tree, new_tree.nodes, outfile_name)