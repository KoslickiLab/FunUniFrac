import sys, os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import src.LP_EMD_helper as LH
import src.EMDU as EMDU
import multiprocessing
from itertools import combinations, repeat
import glob
from src.CONSTANTS import BRITES
from src.LP_EMD_helper import get_descendants
import os
import pandas as pd
import json
import numpy as np
import pydot
from networkx.drawing.nx_pydot import graphviz_layout
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib
from src.KEGG_helpers import make_nodes_readable
#matplotlib.use('MacOSX')
matplotlib.use('Agg')

# The code above is a bit messy, but it should be fairly easy to follow. The main thing to note is that the diffabs
# are stored in a 3D array, where the first two dimensions are the pairwise distances between the samples in the two
# clusters, and the third dimension is the diffabs between the samples in the first two dimensions. The diffabs are
# stored in the same order as the pairwise distances, so the ith diffab corresponds to the ith row and ith column of
# the pairwise distances. The diffabs are also stored in the same order as the nodes in the tree, so the ith diffab
# corresponds to the ith node in the tree. The code above uses the diffabs to set the edge lengths of the tree, and
# then draws the tree using the edge lengths. The edge lengths are set to be the absolute value of the diffabs, and
# the edge color is set to be red if the diffab is positive, and blue if the diffab is negative. The edge width is set
# to be proportional to the absolute value of the diffab.

edge_list_file = "experiments/KOtree/kegg_ko_edge_df_br_ko00001.txt_lengths_n_50_f_10_r_100.txt"
brite = "ko00001"
clusters_file = 'experiments/QIITA_study/dendro_uniform_pw_fu_ko00001.npy_clusters.json'
pw_dist_file = 'experiments/QIITA_study/AAI_pw_fu_ko00001_all.npy'
pw_dist_file_basis = pw_dist_file + ".basis.txt"
diffabs_file = pw_dist_file + ".diffab.npy"
diffabs_basis_file = pw_dist_file + ".diffab.nodes.txt"
std_dev_factor = 0  # number of standard deviations above the mean to plot
# import the clusters
with open(clusters_file) as f:
    clustered_samples = json.load(f)
# import the diffabs
diffabs = np.load(diffabs_file)
# import the file basis (rows and columns of the distances, as well as the indices of the first two dimensions of the
# diffabs)
file_basis = []
with open(pw_dist_file_basis) as f:
    for line in f.readlines():
        file_basis.append(os.path.basename(line.strip()))
diffab_3rd_dim_basis = []
with open(diffabs_basis_file) as f:
    for line in f.readlines():
        diffab_3rd_dim_basis.append(line.strip())
# dictionary mapping files to their index in the basis
file_to_loc_in_basis = {file: i for i, file in enumerate(file_basis)}
cut_level = 0
# get the samples at the cut level
clusters_at_cut_level = clustered_samples[str(cut_level)]
cluster_numbers = list(clusters_at_cut_level.keys())
cluster_1 = clusters_at_cut_level[cluster_numbers[0]]
cluster_2 = clusters_at_cut_level[cluster_numbers[1]]
# get the subset of the diffabs corresponding to the samples in the clusters
diffabs_subset = np.zeros((len(cluster_1), len(cluster_2), len(diffab_3rd_dim_basis)))
for i, sample_1 in enumerate(cluster_1):
    for j, sample_2 in enumerate(cluster_2):
        diffabs_subset[i, j, :] = diffabs[file_to_loc_in_basis[sample_1], file_to_loc_in_basis[sample_2], :]

assert diffabs_subset.shape == (len(cluster_1), len(cluster_2), len(diffab_3rd_dim_basis))
mean_diffab_bet_clusters = np.mean(diffabs_subset, axis=(0, 1))
var_diffab_bet_clusters = np.var(diffabs_subset, axis=(0, 1))

# import the graph
Gdir = LH.import_graph(edge_list_file, directed=True)
# Select the subtree of the KEGG hierarchy rooted at the given BRITE ID
descendants = get_descendants(Gdir, brite)
# add the root node back in (this will only have affect if you are using multiple brites. TODO: not yet implemented)
descendants.add('root')
# select the subgraph from the brite to the leaves
Gdir = Gdir.subgraph(descendants)

# set all edge lengths to zero
for u, v in Gdir.edges():
    Gdir[u][v]['edge_length'] = 0
    Gdir[u][v]['weight'] = .00000000000000000000001
    Gdir[u][v]['color'] = 'k'
# set edge properties for the diffab values. Note that the ith diffab entry should be the edge length between the ith
# and ancestor(i)th nodes in the tree
thresh = np.mean(np.abs(mean_diffab_bet_clusters)) + std_dev_factor * np.std(np.abs(mean_diffab_bet_clusters))
#thresh = 0
important_vertices = []
for i, mean_val in enumerate(mean_diffab_bet_clusters):
    if i != len(mean_diffab_bet_clusters) - 1 and np.abs(mean_val) > thresh:
        u = diffab_3rd_dim_basis[i]
        v = list(Gdir.predecessors(u))[0]
        if mean_val > 0:
            Gdir[v][u]['A'] = mean_val
            Gdir[v][u]['Avar'] = var_diffab_bet_clusters[i]
            Gdir[v][u]['color'] = 'r'
            Gdir[v][u]['weight'] = np.abs(mean_val)
            important_vertices.append(u)
            important_vertices.append(v)
        if mean_val < 0:
            Gdir[v][u]['B'] = -mean_val
            Gdir[v][u]['Bvar'] = var_diffab_bet_clusters[i]
            Gdir[v][u]['color'] = 'b'
            Gdir[v][u]['weight'] = np.abs(mean_val)
            important_vertices.append(u)
            important_vertices.append(v)

T = Gdir.subgraph(important_vertices)
# relabel the nodes to make them human readable
T = make_nodes_readable(T)
# rename nodes to escape : in the names
T = nx.relabel_nodes(T, {node: node.replace(':', '_') for node in T.nodes()})
#pos = graphviz_layout(T, prog="dot")
pos = graphviz_layout(T, prog="twopi")
plt.figure(figsize=(50, 50))
widths = [3000*T[u][v]['weight'] for u, v in T.edges()]
colors = [T[u][v]['color'] for u, v in T.edges()]

new_labels = {}
for u in T.nodes():
    if T.degree(u) <= 1:
        new_labels[u] = ""
    else:
        new_labels[u] = u


node_size_by_degree = [ T.degree(node)*15 for node in T.nodes() ]

nx.draw(T, pos, node_size=node_size_by_degree, with_labels=True, arrows=False, arrowsize=0, width=widths, edge_color=colors, labels=new_labels)
plt.savefig('test_mrh.png')
