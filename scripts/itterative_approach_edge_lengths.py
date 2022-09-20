import numpy as np
from scipy import sparse
from scipy.optimize import lsq_linear
import networkx as nx
from scipy import sparse


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

def get_ancestor(G, node):
    """
    Return the ancestor of a node
    :param G: networkx graph
    :param node: name of a node
    :return: set of nodes
    """
    ancestor = set()
    for n in nx.ancestors(G, node):
        ancestor.add(n)
    return ancestor

def KO_pair_to_index(row_KO, column_KO, KO_to_index):
    """
    Given a pair of KOs, return the corresponding row in the A matrix. This is the same as the index of the pair in the
    distances file, when you order the entries of the pairwise distance matrix in row-major order.
    :param row_KO: KO of the row
    :param column_KO: KO of the column
    :param KO_to_index: dictionary mapping KOs to indices
    :return: index
    """
    row_index = KO_to_index[row_KO]
    column_index = KO_to_index[column_KO]
    return row_index * len(KO_to_index) + column_index

def KO_index_pair_to_index(row_index, column_index, KO_to_index):
    """
    Given a pair of KOs, as represented by their index in KO_to_index, return the corresponding row in the A matrix.
    This is the same as the index of the pair in the
    distances file, when you order the entries of the pairwise distance matrix in row-major order.
    :param row_index: row index of the KO
    :param column_KO: column index of the KO
    :return: index
    """
    return row_index * len(KO_to_index) + column_index

brite = "ko00001"
edge_list = "C:\\Users\\dmk333\\PycharmProjects\\FunUniFrac\\real_data\\kegg_ko_edge_df.txt"
distances_file = "C:\\Users\\dmk333\\PycharmProjects\\FunUniFrac\\real_data\\KOs_sketched_scaled_10_compare_5"
distances_labels_file = "C:\\Users\\dmk333\\PycharmProjects\\FunUniFrac\\real_data\\KOs_sketched_scaled_10_compare_5.labels.txt"
A_file = "C:\\Users\\dmk333\\PycharmProjects\\FunUniFrac\\real_data\\ko00001_KOs_sketched_scaled_10_compare_5_A.npz"
A_basis_file = "C:\\Users\\dmk333\\PycharmProjects\\FunUniFrac\\real_data\\ko00001_KOs_sketched_scaled_10_compare_5_column_basis.txt"
pairwise_dist = np.load(distances_file)
A = sparse.load_npz(A_file)

#import the basis file
with open(A_basis_file, 'r') as f:
    basis = f.read().splitlines()
basis_index = {node: i for i, node in enumerate(basis)}

# import label names
pairwise_dist_KOs = []
with open(distances_labels_file, 'r') as f:
    for line in f.readlines():
        ko = line.strip().split('|')[-1]  # KO number is the last in the pip-delim list
        ko = ko.split(':')[-1]  # remove the ko: prefix
        pairwise_dist_KOs.append(ko)
pairwise_dist_KO_index = {node: i for i, node in enumerate(pairwise_dist_KOs)}

# create the y vector of all pairwise distances
y = []
for ko1 in pairwise_dist_KOs:
    for ko2 in pairwise_dist_KOs:
        y.append(pairwise_dist[pairwise_dist_KO_index[ko1], pairwise_dist_KO_index[ko2]])
y = np.array(y)
# by default, the values are: 0 = most dissimilar, 1 = most similar, so to convert to a distance, we subtract from 1
y = 1 - y

# import the graph, just take the part of the graph containing the KOs and all paths to the given brite
# read in the edge list
G = nx.read_edgelist(edge_list, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
# get the descendants of the brite
descendants = get_descendants(G, brite)
# add the root node back in
descendants.add('root')
# select the subgraph from the brite to the leaves
G = G.subgraph(descendants)
G_undirected = G.to_undirected()

# take a random KO
leaf_nodes = get_leaf_descendants(G, brite)
ko = list(leaf_nodes)[0]
# get the ancestor
ancestors = list(G.predecessors(ko))
if len(ancestors) > 1:
    print("Warning: more than one ancestor, just taking the first")
ancestor = ancestors[0]
# get the descendants
descendants = get_leaf_descendants(G, ancestor)
# subset the descendants to the ones that are in the pairwise distance matrix
descendants = descendants.intersection(set(pairwise_dist_KOs))
# get the portion of the matrix that corresponds to the descendants, but only the ones we have pairwise distances for
descendants_index = [pairwise_dist_KO_index.get(node) for node in descendants]
# each pair of kos in the descendants is a row in the matrix
# so to subset the matrix, I will need to find out which of the rows in the matrix correspond to the descendants indexes
# I can do this by finding the rows that have a 1 in the descendants index column
row_indices = []
for row_index in descendants_index:
    for column_index in descendants_index:
        row_indices.append(KO_index_pair_to_index(row_index, column_index, pairwise_dist_KOs))
col_indices = []
for ko1 in descendants:
    col_indices.append(basis_index[ko1])

A_sub = A[row_indices, :][:, col_indices]
y_sub = y[row_indices]
# solve the least squares problem
res = lsq_linear(A_sub, y_sub, bounds=(0, 1))
# get the edge lengths
edge_lengths = res.x
