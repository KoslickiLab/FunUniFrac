import random
import networkx as nx
import itertools as it
from math import ceil

class BinaryTree(nx.Graph):

    def __init__(self, height):
        nx.Graph.__init__(self)
        nx.balanced_tree(2, height, self)
        edge_lengths = {}
        for e in self.edges():
            edge_lengths[e] = random.randint(1, 10)
        nx.set_edge_attributes(self, values=edge_lengths, name='weight')
        self.leaf_nodes = [node for node in self if nx.degree(self, node) == 1]


    def get_pw_distance(self):
        pw_dist = {}
        for pair in it.combinations(self.leaf_nodes, 2):
            distance = nx.shortest_path_length(self, source=pair[0], target=pair[1], weight='weight')
            pw_dist[(pair[0], pair[1])] = distance
        return pw_dist



tree = BinaryTree(3)
pw_dist = tree.get_pw_distance()
print(tree.edges(data=True))

edge_lengths_solution = {}

def solver(tree, pw_dist, leaf_nodes):
    if len(leaf_nodes) < 3:
        return
    (a,b,c) = leaf_nodes[:3]
    e1 = (pw_dist[(a,c)] - pw_dist[(b,c)] + pw_dist[(a,b)])/2
    e2 = pw_dist[(a,b)] - e1
    edge_lengths_solution[(ceil(a/2)-1,a)] = e1
    edge_lengths_solution[(ceil(b/2)-1,b)] = e2
    for i in range(2, len(leaf_nodes), 2):
        a = leaf_nodes[i]
        b = leaf_nodes[i+1]
        c = leaf_nodes[0]
        e1 = (pw_dist[(c,a)] - pw_dist[(c,b)] + pw_dist[(a, b)]) / 2
        e2 = pw_dist[(a, b)] - e1
        edge_lengths_solution[(ceil(a / 2) - 1, a)] = e1
        edge_lengths_solution[(ceil(b / 2) - 1, b)] = e2
    parent_nodes = []
    for i in range(0, len(leaf_nodes), 2):
        parent_nodes.append(ceil(leaf_nodes[i] / 2) - 1)
    for (a,b) in it.combinations(parent_nodes,2):
        pw_dist[(a, b)] = pw_dist[(a*2+1, b*2+1)] - edge_lengths_solution[(a,a*2+1)] - edge_lengths_solution[(b, b*2+1)]
    solver(tree, pw_dist, parent_nodes)


solver(tree, pw_dist, tree.leaf_nodes)
print(tree.edges(data=True))
print(edge_lengths_solution)