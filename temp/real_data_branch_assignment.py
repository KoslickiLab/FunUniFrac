import networkx as nx
import numpy as np
from itertools import combinations

edge_lengths_solutions = {}

class KeggTree():
    #a tree from edge_list with branch lengths
    def __init__(self, tree):
        self.tree = tree #an nx.DiGraph
        try:
            self.tree.remove_edge('parent','child')
            self.tree.remove_node('parent')
            self.tree.remove_node('child')
        except:
            print("'parent' and 'child' are not in the edge list")
        self.leaf_nodes = [node for node in self.tree if self.tree.out_degree(node) == 0]
        self.pw_dist = {}
        self.get_pw_dist()
        self.size = len(self.tree.nodes())

    def get_pw_dist(self):
        undir_tree = self.tree.to_undirected()
        for pair in combinations(self.leaf_nodes, 2):
            distance = nx.shortest_path_length(undir_tree, source=pair[0], target=pair[1], weight='edge_length')
            self.pw_dist[(pair[0], pair[1])] = self.pw_dist[(pair[1], pair[0])] = distance

    def get_siblings(self, node):
        siblings = set()
        parents = self.tree.predecessors(node)
        for p in parents:
            for n in self.tree.successors(p):
                siblings.add(n)
        return siblings


def get_subtree(edge_list, sub_root):
    '''
    For testing purposes. Find a subtree in subgraph induced by nodes rooted at sub_root
    :param edge_list: edge list file
    :param sub_root: a node on the graph
    :return: a directed tree
    '''
    G = nx.read_edgelist(edge_list, delimiter='\t', nodetype=str, create_using=nx.DiGraph, data=(('edge_length', float),))
    Gsub = G.subgraph(nx.descendants(G, sub_root))
    Gsub_undir = Gsub.to_undirected()
    components = nx.connected_components(Gsub_undir)
    Gsub_components = {}
    for i, c in enumerate(components):
        Gsub_components[i] = c
    sub_tree = G.subgraph(Gsub_components[0])
    print(f"check if is tree: {nx.is_tree(G.subgraph(Gsub_components[0]))}")
    print(f"original size {len(G.nodes)}")
    print(f"subgraph size {len(sub_tree)}")
    sub_tree_root = [n for n in sub_tree.nodes if sub_tree.in_degree(n) == 0]
    print(f"subtree root: {sub_tree_root}")
    return sub_tree, sub_tree_root


def assing_branch_lengths(G, leaf_nodes):
    '''
    Given pair wise distances and a graph, assign branch lengths to the edges
    :param G: KeggTree object
    :return:
    '''
    pw_dist = G.pw_dist
    if len(leaf_nodes) < 3:
        return
    #get sibling + 1 nodes from leaf_nodes
    for n in leaf_nodes:
        siblings = G.get_siblings(n)
        if len(siblings) > 3 and len(siblings) < 6:
            print(f"{len(siblings)} siblings for node {n}")
            for pair in combinations(siblings,2):
                print(f"{pair[0]}, {pair[1]}: distance {pw_dist[(pair[0], pair[1])]}")
        break

def solve_system():
    pass


if __name__ == "__main__":
    edge_list_file = 'kegg_ko_edge_df_br_ko00001.txt_AAI_lengths_n_50_f_10_r_100.txt'
    sub_tree, sub_tree_root = get_subtree(edge_list_file, '09150 Organismal Systems') #around 80 nodes
    real_sub_tree = KeggTree(sub_tree)
    #print(real_sub_tree.pw_dist)
    print(real_sub_tree.leaf_nodes)
    assing_branch_lengths(real_sub_tree, real_sub_tree.leaf_nodes)
