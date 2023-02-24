import networkx as nx
import pandas as pd
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

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
            if pair[0] > pair[1]:
                self.pw_dist[(pair[1], pair[0])] = distance
            else:
                self.pw_dist[(pair[0], pair[1])] = distance

    def get_siblings(self, node):
        siblings = set()
        parents = self.tree.predecessors(node)
        for p in parents:
            for n in self.tree.successors(p):
                siblings.add(n)
        siblings.remove(node) #remove the node itself
        return siblings

    def get_parent(self, node):
        #get one parent
        for p in self.tree.predecessors(node):
            return p

    def get_child(self, node):
        #get one child
        for c in self.tree.successors(node):
            return c


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


def assign_branch_lengths(G, leaf_nodes, pw_dist, edge_lengths_solution):
    '''
    Given pair wise distances and a graph, assign branch lengths to the edges
    :param G: KeggTree object
    :return:
    '''
    if len(leaf_nodes) == 1:
        return
    if len(leaf_nodes) == 2:
        #first assume no single child
        parent = G.get_parent(leaf_nodes[0])
        if leaf_nodes[1] < leaf_nodes[0]:
            leaf_nodes[0], leaf_nodes[1] = leaf_nodes[1], leaf_nodes[0]  # swap
            edge_length = pw_dist[(leaf_nodes[0], leaf_nodes[1])]/2
            edge_lengths_solution[(parent, leaf_nodes[0])] = edge_length
            edge_lengths_solution[(parent, leaf_nodes[1])] = edge_length
        return
    else: #at least 2 parents => at least 3 leaf nodes
        leaf_a = leaf_nodes[0]
        siblings = G.get_siblings(leaf_a)
        leaf_b = siblings.pop() #a and b are siblings
        if leaf_b == leaf_nodes[1]:
            leaf_c = leaf_nodes[2]
        else:
            leaf_c = leaf_nodes[1]
        d1 = pw_dist[(leaf_a, leaf_b)] if leaf_a < leaf_b else pw_dist[(leaf_b, leaf_a)]
        d2 = pw_dist[(leaf_a, leaf_c)] if leaf_a < leaf_c else pw_dist[(leaf_c, leaf_a)]
        d3 = pw_dist[(leaf_b, leaf_c)] if leaf_b < leaf_c else pw_dist[(leaf_c, leaf_b)]
        e1 = (d1 + d2 - d3)/2
        e2 = d1 - e1
        parent = G.get_parent(leaf_a)
        edge_lengths_solution[(parent, leaf_a)] = e1
        edge_lengths_solution[(parent, leaf_b)] = e2
        for i, n in enumerate(leaf_nodes[1:]):
            parent = G.get_parent(n)
            if (parent, n) in edge_lengths_solution:
                continue
            else:
                siblings = G.get_siblings(n)
                sibling = siblings.pop()
                if sibling == leaf_a:
                    another_sibling = siblings.pop()
                    siblings.add(sibling)
                    sibling = another_sibling
                d1 = pw_dist[(n, sibling)] if n < sibling else pw_dist[(sibling, n)]
                if (parent, sibling) in edge_lengths_solution:
                    edge_lengths_solution[(parent, n)] = d1 - edge_lengths_solution[(parent, sibling)]
                else:
                    d2 = pw_dist[(n, leaf_a)] if n < leaf_a else pw_dist[(leaf_a, n)]
                    d3 = pw_dist[(sibling, leaf_a)] if sibling < leaf_a else pw_dist[(leaf_a, sibling)]
                    e1 = (d1 + d2 - d3)/2
                    e2 = d1 - e1
                    edge_lengths_solution[(parent, n)] = e1
                    edge_lengths_solution[(parent, sibling)] = e2
    #get all parent nodes
    parent_nodes = set()
    for n in leaf_nodes:
        parent = G.get_parent(n)
        parent_nodes.add(parent)
    #update pw_dist
    for (a, b) in combinations(parent_nodes, 2):
        if a > b:
            a, b = b, a
        a_child = G.get_child(a)
        b_child = G.get_child(b)
        ab_distance = pw_dist[(a_child, b_child)] if a_child < b_child else pw_dist[(b_child, a_child)]
        pw_dist[(a, b)] = ab_distance - \
                          edge_lengths_solution[(a, a_child)] - edge_lengths_solution[(b, b_child)]
    assign_branch_lengths(G, list(parent_nodes), pw_dist, edge_lengths_solution)

def write_dict_to_file(d, filename):
    with open(filename, 'w') as f:
        for k, v in d.items():
            f.write(str(k) + '\t' + str(v) + '\n')
    return


def visualize_diff(edge_lengths_solution, G, outfile_name):
    '''

    :param edge_lengths_solution:
    :param G: KeggTree object
    :return:
    '''
    inferred_edges = [k for k in list(edge_lengths_solution.keys())]
    inferred_edges.sort()
    inferred_lengths = [edge_lengths_solution[e] for e in inferred_edges]
    actual = {(start, end):v  for (start, end, v) in G.tree.edges(data=True)}
    actual_lengths = [actual[k]['edge_length'] for k in inferred_edges]
    df = pd.DataFrame(columns=['edge', 'inferred_length',  'actual_length'])
    df['edge'] = inferred_edges
    df['actual_length'] = actual_lengths
    df['inferred_length'] = inferred_lengths
    print(df)
    sns.scatterplot(data=df, x='inferred_length', y='actual_length')
    #plt.show()
    plt.savefig(outfile_name)
    return

def L1_norm(edge_lengths_solution, G):
    actual_edges = [(start, end) for (start, end, _) in G.tree.edges(data=True)]
    solution_lengths = [edge_lengths_solution[e] for e in actual_edges]
    actual_lengths = [a['edge_length'] for (_, _, a) in G.tree.edges(data=True)]
    print(actual_lengths[:5])
    print(solution_lengths[:5])
    print(np.linalg.norm(np.array(actual_lengths) - np.array(solution_lengths), 1))
    return


if __name__ == "__main__":
    edge_list_file = 'kegg_ko_edge_df_br_ko00001.txt_AAI_lengths_n_50_f_10_r_100.txt'
    G = nx.read_edgelist(edge_list_file, delimiter='\t', nodetype=str, create_using=nx.DiGraph,
                         data=(('edge_length', float),))
    subgraph_nodes = ['K19768', 'K19773', 'K19765', 'K22878', 'K19767', 'K19764', 'K19774', 'K10141', 'K19766',
                      'K19805', '04212 Longevity regulating pathway - worm', 'K11204', 'K14938', 'K20394',
                      '04213 Longevity regulating pathway - multiple species', 'K19772', 'K13356', 'K17705',
                      '09149 Aging',
                      'K01768', 'K19769', 'K19770', 'K19771', 'K01440', '04211 Longevity regulating pathway']
    sub_tree = G.subgraph(subgraph_nodes)
    real_sub_tree = KeggTree(sub_tree)
    edge_lengths_solution = {}
    pw_dist = real_sub_tree.pw_dist
    assign_branch_lengths(real_sub_tree, real_sub_tree.leaf_nodes, pw_dist, edge_lengths_solution)
    #L1_norm(edge_lengths_solution, real_sub_tree)
    visualize_diff(edge_lengths_solution, real_sub_tree, 'scatter_plot_age09149aging.png')
    #write_dict_to_file(edge_lengths_solution, './edge_lengths_solution_09149aging.txt')
    #nx.write_edgelist(real_sub_tree.tree, './sub_tree_09149aging_original_edge_lengths.txt')