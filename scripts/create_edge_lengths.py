import argparse
import os
from os.path import exists
import numpy as np
import networkx as nx
from scipy import sparse
import pandas as pd

# get all leaf descendants of a certain node
def get_leaf_descendants(G, node):
    descendants = set()
    for n in nx.descendants(G, node):
        if G.out_degree(n) == 0:
            descendants.add(n)
    return descendants


def get_descendants(G, node):
    descendants = set()
    descendants.add(node)
    for n in nx.descendants(G, node):
        descendants.add(n)
    return descendants

# parse arguments
parser = argparse.ArgumentParser(description='This will take the matrix made by graph_to_matrix.py and the all '
                                             'pairwise distance matrix and solve the least squares problem of '
                                             'inferring the edge lengths of the graph.')
parser.add_argument('-e', '--edge_list', help='Input edge list file of the KEGG hierarchy', required=True)
parser.add_argument('-d', '--distances', help='File containing all pairwise distances between KOs. Use sourmash '
                                              'compare', required=True)
parser.add_argument('-o', '--out_dir', help='Output directory')
parser.add_argument('-A', '--A_matrix', help='A matrix file created by graph_to_matrix.py', required=True)
args = parser.parse_args()

edge_list = args.edge_list
out_dir = args.out_dir
distances_file = args.distances
A_matrix_file = args.A_matrix
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
if not exists(A_matrix_file):
    raise FileNotFoundError(f"Could not find {A_matrix_file}")
basis_name = f"{A_matrix_file}_column_basis.txt"
if not exists(basis_name):
    raise FileNotFoundError(f"Could not find the basis file {basis_name} that should accompany {A_matrix_file}.")



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

# create the y vector of all pairwise distances
y = []
for ko1 in pairwise_dist_KOs:
    for ko2 in pairwise_dist_KOs:
        y.append(pairwise_dist[pairwise_dist_KO_index[ko1], pairwise_dist_KO_index[ko2]])
y = np.array(y)

# import the A matrix
A = sparse.load_npz(A_matrix_file)

# solve the least squares problem
x, residuals, rank, s = np.linalg.lstsq(A, y, rcond=None)

# import the basis
with open(basis_name, 'r') as f:
    basis = f.readlines()
    basis = [line.strip() for line in basis]

# import the edge list
edge_list = os.path.join("test_data", "kegg_ko_edge_df.txt")
df = pd.read_csv(edge_list, sep='\t', header=0)
# add a new column for the edge lengths
df['edge_length'] = np.nan
# iterate through the basis and add the edge lengths to the dataframe
for i, tail in enumerate(basis):
    df.loc[df['child'] == tail, 'edge_length'] = x[i]
