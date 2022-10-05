#!/usr/bin/env python
import argparse
import os
import sys
from os.path import exists
import numpy as np
import networkx as nx
from scipy import sparse
import pandas as pd
from scipy.optimize import lsq_linear
import multiprocessing
from itertools import repeat
# relative imports
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from src.CONSTANTS import BRITES
from src.LP_EMD_helper import get_descendants, get_descendant, get_KO_labels_and_index


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


def argument_parser():
    parser = argparse.ArgumentParser(
        description='This will take the matrix made by graph_to_path_matrix.py and the all '
                    'pairwise distance matrix and solve the least squares problem of '
                    'inferring the edge lengths of the graph.')
    parser.add_argument('-e', '--edge_list', help='Input edge list file of the KEGG hierarchy', required=True)
    parser.add_argument('-d', '--distances', help='File containing all pairwise distances between KOs. Use sourmash '
                                                  'compare', required=True)
    parser.add_argument('-o', '--out_file', help='Output file name: edge list with lengths in the last column',
                        required=True)
    parser.add_argument('-A', '--A_matrix', help='A matrix file created by graph_to_path_matrix.py', required=True)
    parser.add_argument('-b', '--brite_id', help='Brite ID of the KEGG hierarchy you want to focus on. Eg. ko00001',
                        required=True)
    parser.add_argument('-n', '--num_iter', help='Number of random selections on which to perform the NNLS',
                        default=100)
    parser.add_argument('-f', '--factor', help='Selects <--factor>*(A.shape[1]) rows for which to do the NNLS',
                        default=5)
    parser.add_argument('-r', '--reg_factor', help='Regularization factor for the NNLS', default=1)
    parser.add_argument('--force', help='Overwrite the output file if it exists', action='store_true')
    parser.add_argument('--distance',
                        help='Flag indicating that the input matrix is a distance (0=identical). If not, it is assumed to be a similarity (1=identical).',
                        action='store_true')
    return parser

def main():
    # parse arguments
    parser = argument_parser()
    args = parser.parse_args()
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
    if not exists(edge_list):
        raise FileNotFoundError(f"Could not find {edge_list}")
    if not exists(distances_file):
        raise FileNotFoundError(f"Could not find {distances_file}")
    if not exists(f"{distances_file}.labels.txt"):
        raise FileNotFoundError(f"Could not find {distances_file}.labels.txt in the same directory as {distances_file}")
    distances_labels_file = f"{distances_file}.labels.txt"
    if os.path.exists(out_file) and not force:
        raise FileExistsError(f"{out_file} already exists. Please delete it or choose another name, or use --force.")
    if not exists(A_matrix_file):
        raise FileNotFoundError(f"Could not find {A_matrix_file}")
    try:
        basis_name = f"{A_matrix_file.removesuffix('_A.npz')}_column_basis.txt"
    except AttributeError:
        if A_matrix_file.endswith('_A.npz'):
            basis_name = f"{A_matrix_file[:-6]}_column_basis.txt"
        else:
            raise FileNotFoundError(f"Could not find the basis file that should accompany {A_matrix_file}. It appears the "
                                    f"matrix file does not end in '_A.npz'. Was it created with graph_to_path_matrix.py?")
    if not exists(basis_name):
        raise FileNotFoundError(f"Could not find the basis file {basis_name} that should accompany {A_matrix_file}.")

    if brite not in BRITES:
        raise ValueError(f"{brite} is not a valid BRITE ID. Choices are: {BRITES}")


    # import pairwise distances
    pairwise_dist = np.load(distances_file)
    # import label names
    pairwise_dist_KOs, pairwise_dist_KO_index = get_KO_labels_and_index(distances_labels_file)

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

    num_threads = min(multiprocessing.cpu_count() // 2, num_iter)
    pool = multiprocessing.Pool(num_threads)

    xs = np.array(pool.map(map_star, zip(range(num_iter), repeat(A), repeat(y), repeat(factor), repeat(reg_factor)), chunksize=num_iter // num_threads))

    # take the average of the solutions
    x = np.mean(xs, axis=0)

    # import the basis
    with open(basis_name, 'r') as f:
        basis = f.readlines()
        basis = [line.strip() for line in basis]

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


if __name__ == '__main__':
    main()
