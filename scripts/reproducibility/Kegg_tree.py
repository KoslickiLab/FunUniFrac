from platform import node

import networkx as nx
import pandas as pd
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
import json
sys.path.append('../data')
sys.path.append('..')


class KeggTree:
    #a tree from edge_list with branch lengths
    def __init__(self, tree, write_merged_file=False, merged_outfile=None):
        self.tree = tree #an nx.DiGraph
        self.root = [node for node in self.tree if self.tree.in_degree(node) == 0][0]
        print(f"Identifying root: {self.root}")
        self.nodes_by_depth = dict()
        #self.group_nodes_by_depth()
        #self.make_full_tree()
        self.leaf_nodes = [node for node in self.tree if self.tree.out_degree(node) == 0]
        self.pw_dist = dict()
        self.needed_pairs = dict()
        self.partners = dict()
        #self.get_pw_dist()
        self.size = len(list(self.tree.nodes()))

    def get_pw_dist(self):
        undir_tree = self.tree.to_undirected()
        for pair in combinations(self.leaf_nodes, 2):
            distance = nx.shortest_path_length(undir_tree, source=pair[0], target=pair[1], weight='edge_length')
            if pair[0] > pair[1]:
                self.pw_dist[(pair[1], pair[0])] = distance
            else:
                self.pw_dist[(pair[0], pair[1])] = distance

    def load_pw_dist_from_file(self, json_file):
        with open(json_file, 'r') as f:
            self.pw_dist = json.load(f)


    def get_siblings(self, node):
        siblings = set()
        parents = self.tree.predecessors(node)
        for p in parents:
            for n in self.tree.successors(p):
                siblings.add(n)
        siblings.discard(node) #discard the node itself
        if len(siblings) == 0:
            print(f"{node} has no siblings")
        return siblings

    def get_sibling(self, node):
        # get one sibling
        for parent in self.tree.predecessors(node):
            for c in self.tree.successors(parent):
                if c != node:
                    return c

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
            self.tree.discard_node(p)
        if write_file:
            outfile = outfile if outfile is not None else 'tree_with_no_single_child.txt'
            nx.write_weighted_edgelist(self.tree, outfile, delimiter='\t')

    def make_full_tree(self):
        #process tree from root down until the deepest level, if any node has no child, add a dummy node
        dummy_node_count = 0
        for i in range(len(self.nodes_by_depth)-1):
            for node in self.nodes_by_depth[i]:
                if not self.tree.successors(node):
                    dummy_node = 'dummy' + str(dummy_node_count)
                    self.tree.add_edge(node, dummy_node, edge_length=0)
                    self.nodes_by_depth[i+1].append(dummy_node)
                    dummy_node_count += 1

    def group_nodes_by_depth(self):
        for node in self.tree.nodes():
            depth = nx.shortest_path_length(self.tree, self.root, node)
            if depth in self.nodes_by_depth:
                self.nodes_by_depth[depth].append(node)
            else:
                self.nodes_by_depth[depth] = [node]

    def _get_first_children(self, nodes):
        '''
        Given a set of nodes, return a set of first child nodes
        :return:
        '''
        first_children = set()
        for node in nodes:
            first_children.add(self.get_child(node))
        return first_children

    def get_needed_pairs(self):
        '''
        Loop over each level from the root and get needed parent of the level below
        :return:
        '''
        for i in range(len(self.nodes_by_depth)-1):
            if i == 0: #all siblings
                self.needed_pairs[1] = dict()
                self.partners[1] = dict()
                if len(self.nodes_by_depth[1]) == 1:
                    self.partners[1] = 1
                elif len(self.nodes_by_depth[1]) == 2:
                    self.needed_pairs[1][(self.nodes_by_depth[1][0], self.nodes_by_depth[1][1])] = 0
                elif len(self.nodes_by_depth[1]) == 3:
                    (node1, node2, node3) = (self.nodes_by_depth[1][0], self.nodes_by_depth[1][1], self.nodes_by_depth[1][2])
                    self.needed_pairs[1][(node1, node2)] = 0
                    self.needed_pairs[1][(node2, node3)] = 0
                    self.needed_pairs[1][(node1, node3)] = 0
                    self.needed_pairs[1][(node3, node1)] = 0
                    self.needed_pairs[1][(node3, node2)] = 0
                    self.needed_pairs[1][(node2, node1)] = 0
                    self.partners[1][node1] = [node2, node3]
                    self.partners[1][node3] = [node2, node1]
                else: # >= 4
                    node_set = set(self.nodes_by_depth[1])
                    first_node = self.nodes_by_depth[1][0]
                    last_node = self.nodes_by_depth[1][-1]
                    second_node = self.nodes_by_depth[1][1]
                    self.needed_pairs[1][(first_node, second_node)] = 0
                    self.needed_pairs[1][(first_node, last_node)] = 0
                    self.needed_pairs[1][(second_node, last_node)] = 0
                    self.partners[1][first_node] = [second_node, last_node]
                    node_set.discard(first_node)
                    node_set.discard(second_node)
                    while len(node_set) > 1:
                        cur_node = node_set.pop()
                        sib = node_set.pop()
                        self.needed_pairs[1][(cur_node, sib)] = 0
                        self.needed_pairs[1][(cur_node, first_node)] = 0
                        self.needed_pairs[1][(sib, first_node)] = 0
                        self.partners[1][cur_node] = [sib, first_node]
                    if len(node_set) == 1:
                        cur_node = node_set.pop() #could be last node, but not first nor second
                        self.needed_pairs[1][(cur_node, second_node)] = 0
                        self.needed_pairs[1][(cur_node, first_node)] = 0
                        self.partners[1][cur_node] = [first_node, second_node]
            else: #i>0 not all children nodes are siblings
                first_children = set()
                self.needed_pairs[i+1] = dict()
                self.partners[i+1] = dict()
                first_node = self.nodes_by_depth[i + 1][0]
                last_node = self.nodes_by_depth[i + 1][-1]
                backup_node = self.nodes_by_depth[i + 1][1]
                node_set = set(self.nodes_by_depth[i + 1])
                #handle first node
                sib = self.get_sibling(first_node)
                if not sib:
                    self.partners[i + 1][first_node] = 0
                else:
                    self.needed_pairs[i + 1][(first_node, sib)] = 0
                    if sib == last_node:
                        self.needed_pairs[i + 1][(first_node, backup_node)] = 0
                        self.needed_pairs[i + 1][(sib, backup_node)] = 0
                        self.partners[i + 1][first_node] = [sib, backup_node]
                    else:
                        self.needed_pairs[i + 1][(first_node, last_node)] = 0
                        self.needed_pairs[i + 1][(sib, last_node)] = 0
                        self.partners[i + 1][first_node] = [sib, last_node]
                    node_set.discard(sib)
                node_set.discard(first_node)
                while len(node_set) > 0:
                    cur_node = node_set.pop()
                    sib = self.get_sibling(cur_node)
                    if not sib:
                        self.partners[i+1][cur_node] = 0
                        continue
                    if sib == first_node:
                        if cur_node == backup_node:
                            another_node = last_node
                        else:
                            another_node = backup_node
                    elif sib == last_node:
                        if cur_node == first_node:
                            another_node = backup_node
                        else:
                            another_node = first_node
                    else:
                        if cur_node == first_node:
                            another_node = last_node
                        else:
                            another_node = first_node
                    self.needed_pairs[i+1][(cur_node, sib)] = 0
                    self.needed_pairs[i+1][(cur_node, another_node)] = 0
                    self.needed_pairs[i+1][(sib, another_node)] = 0
                    self.partners[i+1][cur_node] = (sib, another_node)
                    node_set.discard(sib)
                #first child pairs
                for node in self.nodes_by_depth[i]:
                    child = self.get_child(node)
                    first_children.add(child)
                for (a, b) in combinations(first_children, 2):
                    self.needed_pairs[i+1][(a, b)] = self.needed_pairs[i+1][(b, a)] = 0

    def fill_leaf_pairs_distances(self, pw_dist_file, label_file):
        '''
        Can only be run after get_needed_pairs function is run
        :param pw_dist_file: a .npy file
        :param label_file:
        :return:
        '''
        pw_dist = np.load(pw_dist_file)
        labels = [line.strip() for line in open(label_file, 'r')]
        label_pos = {k: v for v, k in enumerate(labels)}
        for (a, b) in self.needed_pairs[len(self.nodes_by_depth)-1]:
            if a and b:
                a_index = label_pos[a]
                b_index = label_pos[b]
                self.needed_pairs[len(self.nodes_by_depth)-1][(a, b)] = pw_dist[a_index][b_index]

    def solve_branch_lengths(self, edge_length_solutions, level):
        '''
        Solves the branch lengths of the given level
        :param edge_length_solutions:
        :param level:
        :return:
        '''
        if level < 1:
            return
        elif level == 1 and len(self.needed_pairs[1]) == 1:
            (node1, node2) = self.needed_pairs[1]
            parent = self.get_parent(node1)
            edge_length_solutions[(parent, node1)] = \
                    edge_length_solutions[(parent, node2)] = self.needed_pairs[1][(node1, node2)]/2
        else:
            for node in self.partners[level]:
                parent = self.get_parent(node)
                if (parent, node) in edge_length_solutions:
                    continue
                if self.partners[level][node] == 0:
                    edge_length_solutions[(parent, node)] = 0
                else:
                    sib = self.partners[level][node][0]
                    another = self.partners[level][node][1]
                    d1 = self.needed_pairs[level][(node, sib)]
                    d2 = self.needed_pairs[level][(node, another)]
                    d3 = self.needed_pairs[level][(sib, another)]
                    e1 = (d1+d2-d3)/2
                    e2 = d1 - e1
                    edge_length_solutions[(parent, node)] = e1
                    edge_length_solutions[(parent, sib)] = e2
        self.update_needed_pairs(level-1, edge_length_solutions)
        self.solve_branch_lengths(edge_length_solutions, level-1)

    def update_needed_pairs(self, level, edge_length_solutions):
        if level >= len(self.nodes_by_depth):
            print(f"level {level} is too big. Max level allowed is {len(self.nodes_by_depth)-1}.")
            return
        if level < 1:
            return
        for (a, b) in self.needed_pairs[level]:
            a_child = self.get_child(a)
            b_child = self.get_child(b)
            a_child_b_child_dist = self.needed_pairs[level+1][(a_child, b_child)]
            self.needed_pairs[level][(a, b)] = a_child_b_child_dist - edge_length_solutions[(a, a_child)] - \
                edge_length_solutions[(b, b_child)]

    def preprocess_pw_dist(self, pw_dist_file, label_file):
        '''
        Create a nested dictionary to store pw distance. Hopefully this can significantly reduce the size
        :param pw_dist_file: a nested dict, each value has 4 keys: sib, this_sib_dist, this_another_dist,
        sib_another_dist. If no sib, value of sib key is 0
        :param label_file:
        :return:
        '''
        pw_dist_dict = dict()
        leaves_set = set(self.leaf_nodes)
        first_node = leaves_set.pop()
        pw_dist = np.load(pw_dist_file)
        labels = [line.strip() for line in open(label_file, 'r')]
        label_pos = {k: v for v, k in enumerate(labels)}
        first_node_index = label_pos[first_node]
        sibs = self.get_siblings(first_node)
        backup_node = leaves_set.pop() #in case first node is the sib of some node
        backup_node_index = label_pos[backup_node]
        leaves_set.add(backup_node)
        if len(sibs) == 0:
            pw_dist_dict[first_node] = dict()
            pw_dist_dict[first_node]['sib'] = 0
        else:
            sib = sibs.pop()
            pw_dist_dict[first_node] = dict()
            pw_dist_dict[first_node]['sib'] = sib
            another_node = leaves_set.pop()
            if another_node == sib:
                yet_another_node = leaves_set.pop()
                leaves_set.add(another_node)
                another_node = yet_another_node
            sib_index = label_pos[sib]
            another_index = label_pos[another_node]
            pw_dist_dict[first_node]['this_sib_dist'] = pw_dist[first_node_index][sib_index]
            pw_dist_dict[first_node]['this_another_dist'] = pw_dist[first_node_index][another_index]
            pw_dist_dict[first_node]['sib_another_dist'] = pw_dist[sib_index][another_index]
            leaves_set.add(another_node) #important
            leaves_set.discard(sib)

        while len(leaves_set) > 0:
            this_node = leaves_set.pop()
            sibs = self.get_siblings(this_node)
            if len(sibs) == 0:
                pw_dist_dict[this_node] = dict()
                pw_dist_dict[this_node]['sib'] = 0
            else:
                sib = sibs.pop()
                pw_dist_dict[this_node] = dict()
                pw_dist_dict[this_node]['sib'] = sib
                sib_index = label_pos[sib]
                this_index = label_pos[this_node]
                if sib == first_node:
                    pw_dist_dict[this_node]['this_sib_dist'] = pw_dist[this_index][sib_index]  # another = first_node
                    pw_dist_dict[this_node]['this_another_dist'] = pw_dist[this_index][backup_node_index]
                    pw_dist_dict[this_node]['sib_another_dist'] = pw_dist[sib_index][backup_node_index]
                else:
                    pw_dist_dict[this_node]['this_sib_dist'] = pw_dist[this_index][sib_index] #another = first_node
                    pw_dist_dict[this_node]['this_another_dist'] = pw_dist[this_index][first_node_index]
                    pw_dist_dict[this_node]['sib_another_dist'] = pw_dist[sib_index][first_node_index]
                leaves_set.discard(sib)
        print(pw_dist_dict)
        self.pw_dist = pw_dist_dict

    def write_leaf_level_pw_dist(self, file_name):
        '''
        Only meaningful after running fill_leaf_level_pw_dist
        Not very useful now.
        :param file_name: .json file
        :return:
        '''
        with open(file_name, 'w') as f:
            json.dump(self.needed_pairs[len(self.nodes_by_depth)-1], f)

def write_edge_list_preserve_order(edge_length_solutions, original_file, file_name):
    '''
    This function writes out the edge list file after computing edge lengths in the format of parent\child\edge_length.
    This function preserves the order of the original file for easy comparison.
    :param edge_length_solutions: a dict. (parent, child):edge_length
    :param original_file: original file with no length. with heading #parent\t child. no edge_length
    :param file_name: file name to save the new file as
    :return: 
    '''
    with open(original_file, 'r') as f:
        f.readline()
        pairs = f.readlines()
    with open(file_name, 'w') as f:
        f.write("#parent\tchild\tedge_length\n")
        for pair in pairs:
            pair = pair.strip()
            (p, c) = pair.split('\t')
            if (p, c) in edge_length_solutions:
                f.write(f"{p}\t{c}\t{edge_length_solutions[(p, c)]}\n")
            else:
                f.write(f"{p}\t{c}\'NA'\n")


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


def visualize_diff(edge_list_inferred, edge_list_actual, outfile_name):
    '''

    :param edge_lengths_solution:
    :param G: KeggTree object
    :return:
    '''
    inferred_df = pd.read_table(edge_list_inferred, header=0)
    reference_df = pd.read_table(edge_list_actual, header=0)
    df = pd.DataFrame(columns=['edge', 'inferred_length',  'actual_length'])
    df['actual_length'] = reference_df['edge_length']
    df['inferred_length'] = inferred_df['edge_length']
    print(df)
    sns.scatterplot(data=df, x='inferred_length', y='actual_length')
    plt.show()
    #plt.savefig(outfile_name)
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
            f.write(f"{parent}\t{child}\t{edge_length}\n")
        return

def get_KeggTree_from_edgelist(edge_list_file, write_file=False, outfile=None, edge_length=True):
    '''

    :param file:
    :return: KeggTree
    '''
    if not edge_length:
        G = nx.read_edgelist(edge_list_file, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
    else:
        G = nx.read_edgelist(edge_list_file, delimiter='\t', nodetype=str, create_using=nx.DiGraph,
                         data=(('edge_length', float),))
    keggTree = KeggTree(G, write_merged_file=write_file, merged_outfile=outfile)
    return keggTree


def post_process(edge_length_solution, kegg_tree):
    '''
    Even out edges with branch length 0 by sharing half of the ancestor
    :param edge_length_solution:
    :param kegg_tree:
    :return:
    '''
    for i in range(len(kegg_tree.nodes_by_depth)):
        for node in kegg_tree.nodes_by_depth[i]:
            if sum(1 for _ in kegg_tree.tree.successors(node)) == 1: #has single child
                parent = kegg_tree.get_parent(node)
                cur_node = node
                child = kegg_tree.get_child(cur_node)
                if child.startswith('dummy') or edge_length_solution[(cur_node, child)] != 0:
                    continue
                nodes_on_the_way = [cur_node]
                while sum(1 for _ in kegg_tree.tree.successors(cur_node)) == 1:
                    cur_node = kegg_tree.get_child(node)
                    if cur_node.startswith('dummy'):
                        break
                    nodes_on_the_way.append(cur_node)
                ave_edge_length = edge_length_solution[(parent, node)]/len(nodes_on_the_way)
                for node in nodes_on_the_way:
                    edge_length_solution[(parent, node)] = ave_edge_length
                    parent = node






