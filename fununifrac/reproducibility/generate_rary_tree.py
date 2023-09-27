import networkx as nx
import random
from Kegg_tree import write_subgraph_file


#parameters: n, r
R = [2, 3, 4]
N = [50, 100, 500, 1000, 2000, 3000]

for n in N:
    for r in R:
        file_name = f'data/size_{n}_{r}ary_tree.txt'
        tree = nx.full_rary_tree(r, n)
        nx.set_edge_attributes(tree, values=1, name="edge_length")
        weighted_tree = nx.DiGraph()
        for edge in tree.edges.data():
            new_edge = random.random()
            weighted_tree.add_edge(edge[0], edge[1], edge_length=new_edge)
            write_subgraph_file(weighted_tree, weighted_tree.nodes, file_name)