#!/usr/bin/env python
import os
from os.path import exists
import numpy as np
import networkx as nx
from scipy import sparse
import pandas as pd
from scipy.optimize import lsq_linear
import multiprocessing
from itertools import repeat
import src.utils.kegg_db as kegg_db
from src.algorithms.lp_edge_length import get_descendants, get_descendant, get_KO_labels_and_index
import data

def map_func(itr, A, y, factor, reg_factor):
    num_rows = int(factor * A.shape[1])
    row_indices = np.random.choice(A.shape[0], num_rows, replace=False)
    A_small = A[row_indices, :]
    y_small = y[row_indices]
    # append a row of 1's to A_small
    A_small = sparse.vstack([reg_factor * A_small, sparse.csr_matrix(np.ones(A_small.shape[1]))])
    # append a 0 to y_small
    y_small = np.append(reg_factor * y_small, 0)
    # Use lsq_linear to solve the NNLS problem
    res = lsq_linear(A_small, y_small, bounds=(0, 1), verbose=2)
    x = res.x
    return x


def map_star(args):
    return map_func(*args)


def main(args):
    brite = args.brite_id
    edge_list = args.edge_list
    out_file = args.out_file
    distances_file = args.distances
    A_matrix_file = args.A_matrix
    force = args.force
    num_iter = int(args.num_iter)
    factor = int(args.factor)
    reg_factor = float(args.reg_factor)
    isdistance = args.distance
    if num_iter < 1:
        raise ValueError('Number of iterations must be at least 1')
    if factor < 1:
        raise ValueError('Factor must be at least 1')

    # check that the files exist
    edge_list = data.get_data_abspath(edge_list)
    distances_file = data.get_data_abspath(distances_file)
    distances_labels_file = data.get_data_abspath(f"{distances_file}.labels.txt")
    A_matrix_file = data.get_data_abspath(A_matrix_file)
    if A_matrix_file.endswith('_A.npz'):
        basis_name = f"{A_matrix_file[:-len('_A.npz')]}_column_basis.txt"
    else:
        raise FileNotFoundError(f"Could not find the basis file that should accompany {A_matrix_file}. It appears the "
                                    f"matrix file does not end in '_A.npz'. Was it created with graph_to_path_matrix.py?")    
    out_file = data.get_data_abspath(out_file, raise_if_not_found=False)
    if os.path.exists(out_file) and not force:
        raise FileExistsError(f"{out_file} already exists. Please delete it or choose another name, or use --force.")
    if brite not in kegg_db.instance.brites:
        raise ValueError(f"{brite} is not a valid BRITE ID. Choices are: {kegg_db.instance.brites}")

    # import the basis of the A matrix (the shortest path matrix)
    with open(basis_name, 'r') as f:
        basis = f.readlines()
        basis = [line.strip() for line in basis]

    # import pairwise distances
    pairwise_dist = np.load(distances_file)
    # import label names
    pairwise_dist_KOs, pairwise_dist_KO_index = get_KO_labels_and_index(distances_labels_file, basis=basis)

    # create the y vector of all pairwise distances
    y = []
    for ko1 in pairwise_dist_KOs:
        for ko2 in pairwise_dist_KOs:
            y.append(pairwise_dist[pairwise_dist_KO_index[ko1], pairwise_dist_KO_index[ko2]])
    y = np.array(y)
    # by default, the values are: 0 = most dissimilar, 1 = most similar, so to convert to a distance, we subtract from 1
    if not isdistance:
        y = 1 - y
    # import the A matrix
    A = sparse.load_npz(A_matrix_file)
    if A.shape[0] != y.shape[0]:
        raise ValueError(f"The A matrix has {A.shape[0]} rows, but the y vector has {y.shape[0]} elements. "
                         f"Something is wrong.")

    num_threads = 1  # numpy apparently uses all PHYSICAL cores, twice that for hyperthreading
    pool = multiprocessing.Pool(num_threads)

    xs = np.array(pool.map(map_star, zip(range(num_iter), repeat(A), repeat(y), repeat(factor), repeat(reg_factor)), chunksize=num_iter // num_threads))

    # take the average of the solutions
    x = np.mean(xs, axis=0)

    # import the edge list
    df = pd.read_csv(edge_list, sep='\t', header=0)
    # add a new column for the edge lengths
    df['edge_length'] = np.nan
    # iterate through the basis and add the edge lengths to the dataframe
    for i, tail in enumerate(basis):
        df.loc[df['child'] == tail, 'edge_length'] = x[i]
    # remove all the rows that aren't in the basis, so only focusing on the subtree defined by the given brite
    df = df[df['child'].isin(basis)]
    # write the new edge list to file
    df.to_csv(out_file, sep='\t', index=False)
