#!/usr/bin/env python
# Script to compute all pairwise functional unifrac from a directory of functional profiles
import argparse
import os
import sys
from os.path import exists
import numpy as np
import networkx as nx
from scipy import sparse
import pandas as pd
from scipy.optimize import lsq_linear
# relative imports
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    sys.path.append(os.path.dirname(SCRIPT_DIR))
except:
    pass
import src.CONSTANTS as CONSTANTS
import src.LP_EMD_helper as LH
import src.EMDU as EMDU
import multiprocessing
from itertools import combinations, repeat
import glob
from src.CONSTANTS import BRITES
from src.LP_EMD_helper import get_descendants
import logging


def map_func(file, Tint, lint, nodes_in_order, EMDU_index_2_node, abundance_key):
    P = EMDU.functional_profile_to_EMDU_vector(file, EMDU_index_2_node,
                                               abundance_key=abundance_key, normalize=True)
    P_pushed = EMDU.push_up_L1(P, Tint, lint, nodes_in_order)
    return file, P_pushed


def map_star(args):
    return map_func(*args)


def argument_parser():
    # parse arguments
    parser = argparse.ArgumentParser(description="This script will take a directory of sourmash gather results and "
                                                 "a weighted edge list representing the KEGG hierarchy and compute "
                                                 "all pairwise functional unifrac distances.")
    parser.add_argument('-e', '--edge_list',
                        help='Input edge list file of the KEGG hierarchy. Must have lengths in the '
                             'third column.', required=True)
    parser.add_argument('-d', '--directory', help='Directory containing sourmash gather results.', required=True)
    parser.add_argument('-fp', '--file_pattern', help="Pattern to match files in the directory. Default is "
                                                      "'*_gather.csv'", default='*_gather.csv')
    parser.add_argument('-o', '--out_file', help='Output file name: will be a numpy array.', required=True)
    parser.add_argument('-f', '--force', help='Overwrite the output file if it exists', action='store_true')
    parser.add_argument('-a', '--abundance_key',
                        help='Key in the gather results to use for abundance. Default is `median_abund`',
                        default='median_abund')
    parser.add_argument('-t', '--threads', help='Number of threads to use. Default is half the cores available.',
                        default=int(multiprocessing.cpu_count() / 2), type=int)
    parser.add_argument('-b', '--brite', help='Use the subtree of the KEGG hierarchy rooted at the given BRITE ID. '
                                              'eg. brite:ko00001', default="ko00001", type=str, required=False)
    parser.add_argument('--diffab', action='store_true', help='Also return the difference abundance vectors.')
    parser.add_argument('-v', help="Be verbose", action="store_const", dest="loglevel", const=logging.INFO,
                        default=logging.WARNING)
    return parser


def main():
    parser = argument_parser()
    # parse the args and error check
    args = parser.parse_args()
    logging.basicConfig(level=args.loglevel, format='%(asctime)s %(levelname)s: %(message)s')
    edge_list_file = args.edge_list
    directory = args.directory
    out_file = args.out_file
    file_pattern = args.file_pattern
    force = args.force
    abundance_key = args.abundance_key
    num_threads = args.threads
    brite = args.brite
    make_diffab = args.diffab
    if brite not in BRITES:
        raise ValueError(f'Invalid BRITE ID: {brite}. Must be one of {BRITES}')
    if not exists(edge_list_file):
        raise FileNotFoundError(f"Could not find {edge_list_file}")
    if not exists(directory):
        raise FileNotFoundError(f"Could not find {directory}")
    if exists(out_file) and not force:
        raise FileExistsError(f"{out_file} already exists. Please delete or rename it, or use --force.")
    # look for files in that directory with that file pattern
    fun_files = glob.glob(os.path.join(directory, file_pattern))
    if len(fun_files) == 0:
        raise FileNotFoundError(f"No files in {directory} match the pattern {file_pattern}.")
    fun_files = sorted(fun_files)
    logging.info(f"Found {len(fun_files)} files in {directory} matching {file_pattern}")
    logging.info(f"Parsing graph")
    # Parse the graph and get the FunUniFrac objects Tint, length, and edge_list
    Gdir = LH.import_graph(edge_list_file, directed=True)
    # Select the subtree of the KEGG hierarchy rooted at the given BRITE ID
    descendants = get_descendants(Gdir, brite)
    # add the root node back in (this will only have affect if you are using multiple brites. TODO: not yet implemented)
    descendants.add('root')
    # select the subgraph from the brite to the leaves
    Gdir = Gdir.subgraph(descendants)

    logging.info(f"Converting graph into EMDUniFrac format")
    # Then create the inputs for EMDUniFrac
    Tint, lint, nodes_in_order, EMDU_index_2_node = LH.weighted_tree_to_EMDU_input(Gdir)
    # switch the order for convenience
    node_2_EMDU_index = {val: key for key, val in EMDU_index_2_node.items()}

    logging.info(f"Computing pairwise distances in parallel")
    # convert the functional profiles to vectors and push them up in parallel
    Ps_pushed = {}
    pool = multiprocessing.Pool(args.threads)
    results = pool.imap(map_star, zip(fun_files, repeat(Tint), repeat(lint), repeat(nodes_in_order),
                                      repeat(EMDU_index_2_node),
                                      repeat(abundance_key)),
                        chunksize=max(2, len(fun_files) // num_threads))
    pool.close()
    pool.join()
    #results = map(map_star, zip(fun_files))
    for file, P_pushed in results:
        Ps_pushed[file] = P_pushed

    logging.info(f"Computing pairwise distances")
    # Then compute the pairwise distances
    dists = np.zeros((len(fun_files), len(fun_files)))
    # TODO: convert diffabs to something sparse. Eg. https://sparse.pydata.org/en/stable/construct.html
    diffabs = np.zeros((len(fun_files), len(fun_files), len(nodes_in_order)))
    for i, j in combinations(range(len(fun_files)), 2):
        if not make_diffab:
            dists[i, j] = dists[j, i] = EMDU.EMD_L1_on_pushed(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
        else:
            Z, diffab = EMDU.EMD_L1_and_diffab_on_pushed(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
            dists[i, j] = dists[j, i] = Z
            diffabs[i, j, :] = diffab
            diffabs[j, i, :] = -diffab

    logging.info(f"Saving results")
    # save the distances
    np.save(out_file, dists)
    logging.info(f"Saved distances to {out_file}")
    if make_diffab:
        logging.info("Converting to an npy file")
        np.save(out_file + '.diffab.npy', diffabs)
        logging.info(f"Saved difference abundance diffabs to {out_file}.diffab.npy")
        # save the nodes in order, but using the original node names
        with open(out_file + '.diffab.nodes.txt', 'w') as f:
            for node in nodes_in_order:
                f.write(f'{EMDU_index_2_node[node]}\n')
    # save the basis
    with open(out_file + '.basis.txt', 'w') as f:
        for file in fun_files:
            f.write(f"{file}\n")


if __name__ == '__main__':
    main()



