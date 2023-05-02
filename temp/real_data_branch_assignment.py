import networkx as nx
import pandas as pd
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../data')
sys.path.append('..')
import argparse


class KeggTree():
    #a tree from edge_list with branch lengths
    def __init__(self, tree, write_merged_file=False, merged_outfile=None):
        self.tree = tree #an nx.DiGraph
        self.root = [node for node in self.tree if self.tree.in_degree(node) == 0][0]
        self.nodes_by_depth = dict()
        self.group_nodes_by_depth()
        self.make_full_tree()
        self.leaf_nodes = self.nodes_by_depth[len(self.nodes_by_depth)-1]
        self.pw_dist = {}
        self.get_pw_dist()
        self.size = len(list(self.tree.nodes()))

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
        if len(siblings) == 0:
            print(f"{node} has no siblings")
        return siblings

    def get_parent(self, node):
        #get one parent
        for p in self.tree.predecessors(node):
            return p

    def get_child(self, node):
        #get one child
        for c in self.tree.successors(node):
            return c

    def merge_single_child_branches(self, write_file=False, outfile=None):
        single_child_parents = [node for node in self.tree if self.tree.out_degree(node) == 1]
        #assume root is not a single parent first
        if len(single_child_parents) == 0:
            return
        print("single_child_parents: ", single_child_parents)
        single_child_parents.reverse()
        for p in single_child_parents:
            single_child = self.get_child(p)
            single_child_grand = self.get_parent(p)
            new_edge_length = self.tree.get_edge_data(p, single_child)['edge_length'] + \
                              self.tree.get_edge_data(single_child_grand, p)['edge_length']
            self.tree.add_edge(single_child_grand, single_child, edge_length=new_edge_length)
            self.tree.remove_node(p)
        if write_file:
            outfile = outfile if outfile is not None else 'tree_with_no_single_child.txt'
            nx.write_weighted_edgelist(self.tree, outfile, delimiter='\t')

    def make_full_tree(self):
        #process tree from root down until the deepest level, if any node has no child, add a dummy node
        dummy_node_count = 0
        for i in range(len(self.nodes_by_depth)):
            for node in self.nodes_by_depth[i]:
                if not self.tree.successors(node):
                    dummy_node = 'dummy' + str(dummy_node_count)
                    self.tree.add_edge(node, dummy_node, edge_length=0)
                    self.nodes_by_depth[i].append(dummy_node)
                    dummy_node_count += 1

    def group_nodes_by_depth(self):
        for node in self.tree.nodes():
            depth = nx.shortest_path_length(self.tree, self.root, node)
            if depth in self.nodes_by_depth:
                self.nodes_by_depth[depth].append(node)
            else:
                self.nodes_by_depth[depth] = [node]



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
        #root
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
        parent_nodes = set()
        leaf_a = leaf_nodes[0]
        siblings = G.get_siblings(leaf_a)
        if len(siblings) == 0:
            parent = G.get_parent(leaf_a)
            parent_nodes.add(parent)
            edge_lengths_solution[(parent, leaf_a)] = 0
        else:
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
            parent_nodes.add(parent)
            edge_lengths_solution[(parent, leaf_a)] = e1
            edge_lengths_solution[(parent, leaf_b)] = e2
            if G.get_parent(leaf_c) == parent:
                e3 = d2 - e1
                edge_lengths_solution[(parent, leaf_c)] = e3
        for i, n in enumerate(leaf_nodes[1:]):
            parent = G.get_parent(n)
            parent_nodes.add(parent)
            if (parent, n) in edge_lengths_solution:
                continue
            else:
                siblings = G.get_siblings(n)
                if len(siblings) == 0:
                    edge_lengths_solution[(parent, n)] = 0
                else:
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
    actual = {(start, end): v for (start, end, v) in G.tree.edges(data=True)}
    actual_lengths = [actual[k]['edge_length'] for k in inferred_edges]
    df = pd.DataFrame(columns=['edge', 'inferred_length',  'actual_length'])
    df['edge'] = inferred_edges
    df['actual_length'] = actual_lengths
    df['inferred_length'] = inferred_lengths
    print(df.to_string())
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

def write_subgraph_file(G, subgraph_nodes, out_file):
    '''
    Not really going to use because most subgraph files are created making jupyter notebook.
    Just include for completeness
    :param G:
    :param subgraph_nodes:
    :param out_file:
    :return:
    '''
    sub_tree = G.subgraph(subgraph_nodes)
    with open(out_file, 'w') as f:
        f.write("#parent\tchild\tedge_length\n")
        for edge in sub_tree.edges(data=True):
            parent = edge[0]
            child = edge[1]
            edge_length = edge[2]['edge_length']
            print(parent, child, edge_length)
            f.write(f"{parent}\t{child}\t{edge_length}\n")
        return

def get_KeggTree_from_edgelist(edge_list_file, write_file=False, outfile=None):
    '''

    :param file:
    :return: KeggTree
    '''
    G = nx.read_edgelist(edge_list_file, delimiter='\t', nodetype=str, create_using=nx.DiGraph,
                         data=(('edge_length', float),))
    keggTree = KeggTree(G, write_merged_file=write_file, merged_outfile=outfile)
    return keggTree


def main():
    parser = argparse.ArgumentParser(description="At this point the script takes in an edge-list file WITH branch lengths"
                                                 " and tries to recompute the edge lengths using the algorithm and compare "
                                                 "the accuracy. With that being said, any sub tree would need to be created"
                                                 " prior to passing into the script. ")
    parser.add_argument('-e', '--edge_list',
                        help='Input edge list file of the KEGG hierarchy. Must have lengths in the '
                             'third column.', required=True)
    parser.add_argument('-s', '--save', help='Path to save the output file.', default='edge_length_solution.png')

    args = parser.parse_args()
    #edge_list = get_data_abspath(args.edge_list)
    edge_list = args.edge_list
    kegg_tree = get_KeggTree_from_edgelist(edge_list)
    edge_lengths_solution = {}
    pw_dist = kegg_tree.pw_dist
    print("preparation done.")
    assign_branch_lengths(kegg_tree, kegg_tree.leaf_nodes, pw_dist, edge_lengths_solution)

    visualize_diff(edge_lengths_solution, kegg_tree, args.save)


if __name__ == "__main__":
    main()