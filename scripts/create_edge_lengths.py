#!/usr/bin/env python
import argparse
import os
from os.path import exists
import numpy as np
import networkx as nx
from scipy import sparse
import pandas as pd


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
parser.add_argument('-b', '--brite_id', help='Brite ID of the KEGG hierarchy you want to focus on. Eg. ko00001',
                    required=True)
args = parser.parse_args()
brite = args.brite_id
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
basis_name = f"{A_matrix_file.removesuffix('_A.npz')}_column_basis.txt"
if not exists(basis_name):
    raise FileNotFoundError(f"Could not find the basis file {basis_name} that should accompany {A_matrix_file}.")

if brite not in BRITES:
    raise ValueError(f"{brite} is not a valid BRITE ID. Choices are: {BRITES}")


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
# remove all the rows that aren't in the basis, so only focusing on the subtree defined by the given brite
df = df[df['child'].isin(basis)]
# write the new edge list to file
df.to_csv(os.path.join(out_dir, f"kegg_ko_edge_df_with_edge_lengths_{brite}_{os.path.basename(distances_file)}.txt"), sep='\t',
          index=False)
