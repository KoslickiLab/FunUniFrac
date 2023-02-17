import networkx as nx
import numpy as np
from itertools import combinations

from networkx import DiGraph


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
        for pair in combinations(self.leaf_nodes, 2):
            distance = nx.shortest_path_length(self.tree, 'parent', 'child')
            self.pw_dist[(pair[0], pair[1])] = distance

    def reassign_branch_lengths(self):
        pass

edge_list_file = 'kegg_ko_edge_df.txt'




