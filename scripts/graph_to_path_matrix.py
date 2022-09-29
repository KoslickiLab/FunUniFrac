#!/usr/bin/env python
import argparse
import os
from os.path import exists
import numpy as np
import networkx as nx
from scipy import sparse
try:
    from blist import blist
except ModuleNotFoundError:
    print("Warning: Could not import blist. Please install blist to speed up the path matrix calculation.")

BRITES = ['br08901', 'br08902', 'br08904', 'ko00001', 'ko00002', 'ko00003', 'br08907',
          'ko01000', 'ko01001', 'ko01009', 'ko01002', 'ko01003', 'ko01005', 'ko01011',
          'ko01004', 'ko01008', 'ko01006', 'ko01007', 'ko00199', 'ko00194', 'ko03000',
          'ko03021', 'ko03019', 'ko03041', 'ko03011', 'ko03009', 'ko03016', 'ko03012',
          'ko03110', 'ko04131', 'ko04121', 'ko03051', 'ko03032', 'ko03036', 'ko03400',
          'ko03029', 'ko02000', 'ko02044', 'ko02042', 'ko02022', 'ko02035', 'ko03037',
          'ko04812', 'ko04147', 'ko02048', 'ko04030', 'ko04050', 'ko04054', 'ko03310',
          'ko04040', 'ko04031', 'ko04052', 'ko04515', 'ko04090', 'ko01504', 'ko00535',
          'ko00536', 'ko00537', 'ko04091', 'ko04990', 'ko03200', 'ko03210', 'ko03100',
          'br08001', 'br08002', 'br08003', 'br08005', 'br08006', 'br08007', 'br08009',
          'br08021', 'br08201', 'br08202', 'br08204', 'br08203', 'br08303', 'br08302',
          'br08301', 'br08313', 'br08312', 'br08304', 'br08305', 'br08331', 'br08330',
          'br08332', 'br08310', 'br08307', 'br08327', 'br08311', 'br08402', 'br08401',
          'br08403', 'br08411', 'br08410', 'br08420', 'br08601', 'br08610', 'br08611',
          'br08612', 'br08613', 'br08614', 'br08615', 'br08620', 'br08621', 'br08605',
          'br03220', 'br03222', 'br01610', 'br01611', 'br01612', 'br01613', 'br01601',
          'br01602', 'br01600', 'br01620', 'br01553', 'br01554', 'br01556', 'br01555',
          'br01557', 'br01800', 'br01810', 'br08020', 'br08120', 'br08319', 'br08329',
          'br08318', 'br08328', 'br08309', 'br08341', 'br08324', 'br08317', 'br08315',
          'br08314', 'br08442', 'br08441', 'br08431']



# get all leaf descendants of a certain node
def get_leaf_descendants(G, node):
    descendants = set()
    for n in nx.descendants(G, node):
        if G.out_degree(n) == 0:
            descendants.add(n)
    return descendants


def get_descendants(G, node):
    """
    Return all descendants of a node, including the node itself.
    :param G: networkx graph
    :param node: name of a node
    :return: set of nodes
    """
    descendants = set()
    descendants.add(node)
    for n in nx.descendants(G, node):
        descendants.add(n)
    return descendants

def get_descendant(graph, v1, v2):
    """
    of the two nodes v1 and v2, ASSUMED TO BE ADJACENT, find out which one is the descendant of the other
    :param graph: networkx graph, directed
    :param v1: node name
    :param v2: node name
    :return: descendant node name
    """
    if v1 in nx.descendants(graph, v2):
        return v1
    elif v2 in nx.descendants(graph, v1):
        return v2
    else:
        raise ValueError("Nodes are not adjacent")

# parse arguments
parser = argparse.ArgumentParser(description='Given the KEGG hierarchy, first this will select the subtree consisting of'
                                             ' all the ancestors of the given brite_id. Then, it will create a matrix '
                                             'where the (i,j) entry will be 1 iff for the ith pair of KOs in the '
                                             '--distances file: (KO1, KO2), edge j is on the shortest path from '
                                             'KO1 to KO2')
parser.add_argument('-e', '--edge_list', help='Input edge list file of the KEGG hierarchy', required=True)
parser.add_argument('-d', '--distances', help='File containing all pairwise distances between KOs. Use sourmash '
                                              'compare', required=True)
parser.add_argument('-o', '--out_dir', help='Output directory', required=True)
parser.add_argument('-b', '--brite_id', help='Brite ID of the KEGG hierarchy you want to focus on. Eg. ko00001', required=True)
args = parser.parse_args()

edge_list = args.edge_list
out_dir = args.out_dir
distances_file = args.distances
brite = args.brite_id
# check that the files exist
if not exists(edge_list):
    raise FileNotFoundError(f"Could not find {edge_list}")
if not exists(distances_file):
    raise FileNotFoundError(f"Could not find {distances_file}")
if not exists(f"{distances_file}.labels.txt"):
    raise FileNotFoundError(f"Could not find {distances_file}.labels.txt in the same directory as {distances_file}")
distances_labels_file = f"{distances_file}.labels.txt"
if not exists(out_dir):
    os.mkdir(out_dir)

# check if brite is legit
if brite not in BRITES:
    raise ValueError(f"{brite} is not a valid BRITE ID. Choices are: {BRITES}")

matrix_name = f"{brite}_{os.path.basename(distances_file)}_A.npz"
basis_name = f"{brite}_{os.path.basename(distances_file)}_column_basis.txt"

########################################################################################################################
# Let's do the following: since I've already computed all pairwise distances, we can just make a large
# least squares problem fitting the tree distances to the pairwise distances
# Let's get the matrix describing which edges are traversed between all pairs of nodes
# This is a sparse matrix, so we'll need to use scipy.sparse

# read in the edge list
G = nx.read_edgelist(edge_list, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
# get the descendants of the brite
descendants = get_descendants(G, brite)
# add the root node back in (this will only have affect if you are using multiple brites. TODO: not yet implemented)
descendants.add('root')
# select the subgraph from the brite to the leaves
G = G.subgraph(descendants)
G_undirected = G.to_undirected()

# set the basis for the tree, which is an ordering of the edges. I'll identify each edge by its terminal node
basis = [x for x in G.nodes()]
basis_index = {node: i for i, node in enumerate(basis)}

# Save the basis
with open(f"{os.path.join(out_dir, basis_name)}", 'w') as f:
    f.write('\n'.join(basis))

# Let's go with a csr_array to store the (pairwise_distance, edge_list) matrix
try:
    data = blist([])
    row_inds = blist([])
    col_inds = blist([])
except NameError:
    data = []
    row_inds = []
    col_inds = []
# import pairwise distances
pairwise_dist = np.load(distances_file)
# import label names
pairwise_dist_KOs = []
with open(distances_labels_file, 'r') as f:
    for line in f.readlines():
        ko = line.strip().split('|')[-1]  # KO number is the last in the pip-delim list
        ko = ko.split(':')[-1]  # remove the ko: prefix
        pairwise_dist_KOs.append(ko)
pairwise_dist_KO_index = {node: i for i, node in enumerate(pairwise_dist_KOs)}

########
# for all pairs of edges, find which of the two connected nodes is the descendant of the other
# then, add the distance to the matrix
edge_2_descendant = {}
edges = list(G.edges())
for i, (v1, v2) in enumerate(edges):
    desc = get_descendant(G, v1, v2)
    edge_2_descendant[(v1, v2)] = desc
    edge_2_descendant[(v2, v1)] = desc

# iterate over the pairs of nodes for which we have pairwise distances
row_ind = -1
for node_i in pairwise_dist_KOs:
    # Some KOs may not be in the subtree selected, so append a row of zeros for those (i.e. don't add to the data list).
    # That won't affect the least squares fit
    if node_i in G_undirected:
        paths = nx.single_source_dijkstra_path(G_undirected, node_i)
    else:
        paths = dict()
    for node_j in pairwise_dist_KOs:
        row_ind += 1
        # if the nodes are the same, skip since this row of the matrix is all zeros
        if node_i != node_j:
            # get the shortest path between the two nodes
            try:
                path = paths[node_j]
                # get the index of each path element in the basis
                # FIXME: there's a problem here since the paths are returned as a list of nodes, but the basis is a
                #  list of edges (edge labeled with the terminal node). So, we need to get the edge that connects
                path_tuples = [(path[i], path[i+1]) for i in range(len(path)-1)]
                # for each tuple, find which is the descendant
                path_edges = [edge_2_descendant[(x[0], x[1])] for x in path_tuples]
                path_indices = [basis_index[node] for node in path_edges]
                # set the corresponding entries in the sparse matrix to 1
                for path_index in path_indices:
                    data.append(1)
                    row_inds.append(row_ind)
                    col_inds.append(path_index)
            except KeyError:
                # if there is no path between the two nodes, skip
                pass
        if row_ind % 100000 == 0:
            print(f"Finished {row_ind}/{len(pairwise_dist_KOs)**2} rows")

A = sparse.csr_matrix((data, (row_inds, col_inds)), shape=(len(pairwise_dist_KOs)**2, len(basis)))
# save the sparse matrix
sparse.save_npz(os.path.join(out_dir, matrix_name), A)
