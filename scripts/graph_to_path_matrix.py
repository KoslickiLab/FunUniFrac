#!/usr/bin/env python
import argparse
import os
import sys
from os.path import exists
import numpy as np
import networkx as nx
from scipy import sparse
import multiprocessing
try:
    from blist import blist
except ModuleNotFoundError:
    print("Warning: Could not import blist. Please install blist to speed up the path matrix calculation.")
# relative imports
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    sys.path.append(os.path.dirname(SCRIPT_DIR))
except:
    pass
from src.CONSTANTS import BRITES
from src.LP_EMD_helper import get_descendants, get_descendant, get_KO_labels_and_index
from itertools import repeat


def map_func(node_i, G):
    if node_i in G:
        return {node_i: nx.single_source_dijkstra_path(G, node_i)}
    else:
        return {}


def map_star(args):
    return map_func(*args)


def argument_parser():
    parser = argparse.ArgumentParser(
        description='Given the KEGG hierarchy, first this will select the subtree consisting of'
                    ' all the ancestors of the given brite_id. Then, it will create a matrix '
                    'where the (i,j) entry will be 1 iff for the ith pair of KOs in the '
                    '--distances file: (KO1, KO2), edge j is on the shortest path from '
                    'KO1 to KO2')
    parser.add_argument('-e', '--edge_list', help='Input edge list file of the KEGG hierarchy', required=True)
    parser.add_argument('-d', '--distances', help='File containing all pairwise distances between KOs. Use sourmash '
                                                  'compare', required=True)
    parser.add_argument('-o', '--out_dir', help='Output directory', required=True)
    parser.add_argument('-b', '--brite_id', help='Brite ID of the KEGG hierarchy you want to focus on. Eg. ko00001',
                        required=True)
    return parser


def main():
    # parse arguments
    parser = argument_parser()
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
    pairwise_dist_KOs, pairwise_dist_KO_index = get_KO_labels_and_index(distances_labels_file, basis=basis)

    ########
    # for all pairs of edges, find which of the two connected nodes is the descendant of the other
    edge_2_descendant = {}
    edges = list(G.edges())
    num_edges = len(edges)
    for i, (v1, v2) in enumerate(edges):
        if i % 100 == 0:
            print(f"descendant iteration: {i}/{num_edges}")
        desc = get_descendant(G, v1, v2)
        edge_2_descendant[(v1, v2)] = desc
        edge_2_descendant[(v2, v1)] = desc

    # iterate over the pairs of nodes for which we have pairwise distances
    # In parallel, get all shortest paths
    print("Getting all shortest paths...")
    num_processes = multiprocessing.cpu_count() // 2
    pool = multiprocessing.Pool(num_processes)
    paths_list = pool.imap(map_star, zip(pairwise_dist_KOs, repeat(G_undirected)), chunksize=max(1, len(pairwise_dist_KOs) // num_processes))
    # The results should be ordered the same as the pairwise_dist_KOs
    print("Done getting all shortest paths")


    row_ind = -1
    for i, node_i in enumerate(pairwise_dist_KOs):
        # Some KOs may not be in the subtree selected, so append a row of zeros for those (i.e. don't add to the data list).
        # That won't affect the least squares fit
        if node_i in G_undirected:
            paths = next(paths_list)
        else:
            paths = dict()
        for node_j in pairwise_dist_KOs:
            row_ind += 1
            # if the nodes are the same, skip since this row of the matrix is all zeros
            if node_i != node_j:
                # get the shortest path between the two nodes
                try:
                    path = paths[node_i][node_j]
                    # get the index of each path element in the basis
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
    print(f"Saving sparse matrix to {os.path.join(out_dir, matrix_name)}")
    sparse.save_npz(os.path.join(out_dir, matrix_name), A)


if __name__ == '__main__':
    main()
