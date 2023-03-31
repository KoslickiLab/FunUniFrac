#!/usr/bin/env python
# Script to compute all pairwise functional unifrac from a directory of functional profiles
import argparse
import os
import numpy as np
import networkx as nx
from scipy import sparse
import pandas as pd
from scipy.optimize import lsq_linear
import src.LP_EMD_helper as LH
import src.EMDU as EMDU
import multiprocessing
from itertools import combinations, repeat
from src.CONSTANTS import BRITES
from src.LP_EMD_helper import get_descendants
import logging
from blist import blist
import sparse
import pickle
import data


def map_func(file, Tint, lint, nodes_in_order, EMDU_index_2_node, abundance_key, unweighted, is_L2):
    P = EMDU.functional_profile_to_EMDU_vector(file, EMDU_index_2_node,
                                               abundance_key=abundance_key, normalize=True)
    if unweighted:
        # if entries are positive, set to 1
        P[P > 0] = 1
    if is_L2:
        P_pushed = EMDU.push_up_L2(P, Tint, lint, nodes_in_order)
    else:
        P_pushed = EMDU.push_up_L1(P, Tint, lint, nodes_in_order)
    return file, P_pushed


def map_star(args):
    return map_func(*args)


def construct_graph():
    """
    logging.info(f"Parsing graph")
    # Parse the graph and get the FunUniFrac objects Tint, length, and edge_list
    Gdir = LH.import_graph(edge_list_file, directed=True)
    # Select the subtree of the KEGG hierarchy rooted at the given BRITE ID
    descendants = get_descendants(Gdir, brite)
    # add the root node back in (this will only have affect if you are using multiple brites. TODO: not yet implemented)
    descendants.add('root')
    # select the subgraph from the brite to the leaves
    Gdir = Gdir.subgraph(descendants)
    """
    return None


def interface_emd():
    """
    logging.info(f"Converting graph into EMDUniFrac format")
    # Then create the inputs for EMDUniFrac
    Tint, lint, nodes_in_order, EMDU_index_2_node = LH.weighted_tree_to_EMDU_input(Gdir)
    # switch the order for convenience
    node_2_EMDU_index = {val: key for key, val in EMDU_index_2_node.items()}
    """

def run_in_parallel():
    """
    logging.info(f"Computing pairwise distances in parallel")
    # convert the functional profiles to vectors and push them up in parallel
    Ps_pushed = {}
    pool = multiprocessing.Pool(args.threads)
    results = pool.imap(map_star, zip(fun_files, repeat(Tint), repeat(lint), repeat(nodes_in_order),
                                      repeat(EMDU_index_2_node),
                                      repeat(abundance_key), repeat(unweighted), repeat(is_L2)),
                        chunksize=max(2, len(fun_files) // num_threads))
    pool.close()
    pool.join()
    #results = map(map_star, zip(fun_files))
    for file, P_pushed in results:
        Ps_pushed[file] = P_pushed
    """

def pairwise_dist_cal():
    """
    logging.info(f"Computing pairwise distances")
    # Then compute the pairwise distances
    dists = np.zeros((len(fun_files), len(fun_files)))
    i_coords = blist([])
    j_coords = blist([])
    k_coords = blist([])
    data_vals = blist([])
    diffab_dims = (len(fun_files), len(fun_files), len(nodes_in_order))
    for i, j in combinations(range(len(fun_files)), 2):
        if not make_diffab:
            if not is_L2:
                dists[i, j] = dists[j, i] = EMDU.EMD_L1_on_pushed(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
            else:
                dists[i, j] = dists[j, i] = EMDU.EMD_L2_on_pushed(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
        else:
            if not is_L2:
                Z, diffab = EMDU.EMD_L1_and_diffab_on_pushed(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
            else:
                Z, diffab = EMDU.EMD_L2_and_diffab_on_pushed(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
            dists[i, j] = dists[j, i] = Z
            nonzero_diffab_locs = np.nonzero(diffab)[0]
            i_coords.extend([i] * len(nonzero_diffab_locs))
            j_coords.extend([j] * len(nonzero_diffab_locs))
            k_coords.extend(nonzero_diffab_locs)
            data_vals.extend(diffab[nonzero_diffab_locs])
    if make_diffab:
        logging.info(f"Converting diffabs to sparse matrix")
        coords = [i_coords, j_coords, k_coords]
        diffabs_sparse = sparse.COO(coords, data_vals, shape=diffab_dims)
        # add the negative transpose ([[0, 1], [0, 0]] -> [[0, 1], [-1, 0]])
        diffabs_sparse = diffabs_sparse - diffabs_sparse.transpose(axes=(1, 0, 2))
    """

def pre_processing():
    # I will move this part to scripts
    """
    if brite not in BRITES:
        raise ValueError(f'Invalid BRITE ID: {brite}. Must be one of {BRITES}')
    edge_list_file = data.get_data_abspath(edge_list_file)
    out_file = data.get_data_abspath(out_file, raise_if_not_found=False)
    if os.path.exists(out_file) and not force:
        raise FileExistsError(f"{out_file} already exists. Please delete it or choose another name, or use --force.")
    # look for files in that directory with that file pattern
    fun_files = data.get_data_abspaths(file_pattern)
    fun_files = sorted(fun_files)
    """

def post_processing():
    # I will move this part to scripts
    """
    logging.info(f"Saving results")
    # save the distances
    np.save(out_file, dists)
    logging.info(f"Saved distances to {out_file}")
    if make_diffab:
        logging.info("Saving diffabs as a sparse npz file")
        diffabs_out_file = out_file + '.diffab.npz'
        #np.save(out_file + '.diffab.npy', diffabs)
        sparse.save_npz(diffabs_out_file, diffabs_sparse)
        logging.info(f"Saved difference abundance diffabs to {diffabs_out_file}")
        # save the nodes in order, but using the original node names
        with open(diffabs_out_file + '.nodes.txt', 'w') as f:
            for node in nodes_in_order:
                f.write(f'{EMDU_index_2_node[node]}\n')
    # save the basis
    with open(out_file + '.basis.txt', 'w') as f:
        for file in fun_files:
            f.write(f"{file}\n")
    # if asked, save the pushed profiles
    if save_Ppushed:
        logging.info(f"Saving pushed profiles")
        ppushed_out_file = out_file + '.ppushed.pkl'
        with open(ppushed_out_file, 'wb') as f:
            pickle.dump(Ps_pushed, f)
        logging.info(f"Saved pushed profiles to {ppushed_out_file}")
    """


def main(args):
    logging.basicConfig(level=args.loglevel, format='%(asctime)s %(levelname)s: %(message)s')
    edge_list_file = args.edge_list
    out_file = args.out_file
    file_pattern = args.file_pattern
    force = args.force
    abundance_key = args.abundance_key
    num_threads = args.threads
    brite = args.brite
    make_diffab = args.diffab
    unweighted = args.unweighted
    is_L2 = args.L2
    save_Ppushed = args.Ppushed
    if brite not in BRITES:
        raise ValueError(f'Invalid BRITE ID: {brite}. Must be one of {BRITES}')
    edge_list_file = data.get_data_abspath(edge_list_file)
    out_file = data.get_data_abspath(out_file, raise_if_not_found=False)
    if os.path.exists(out_file) and not force:
        raise FileExistsError(f"{out_file} already exists. Please delete it or choose another name, or use --force.")
    # look for files in that directory with that file pattern
    fun_files = data.get_data_abspaths(file_pattern)
    fun_files = sorted(fun_files)
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
                                      repeat(abundance_key), repeat(unweighted), repeat(is_L2)),
                        chunksize=max(2, len(fun_files) // num_threads))
    pool.close()
    pool.join()
    #results = map(map_star, zip(fun_files))
    for file, P_pushed in results:
        Ps_pushed[file] = P_pushed

    logging.info(f"Computing pairwise distances")
    # Then compute the pairwise distances
    dists = np.zeros((len(fun_files), len(fun_files)))
    i_coords = blist([])
    j_coords = blist([])
    k_coords = blist([])
    data_vals = blist([])
    diffab_dims = (len(fun_files), len(fun_files), len(nodes_in_order))
    for i, j in combinations(range(len(fun_files)), 2):
        if not make_diffab:
            if not is_L2:
                dists[i, j] = dists[j, i] = EMDU.EMD_L1_on_pushed(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
            else:
                dists[i, j] = dists[j, i] = EMDU.EMD_L2_on_pushed(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
        else:
            if not is_L2:
                Z, diffab = EMDU.EMD_L1_and_diffab_on_pushed(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
            else:
                Z, diffab = EMDU.EMD_L2_and_diffab_on_pushed(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
            dists[i, j] = dists[j, i] = Z
            nonzero_diffab_locs = np.nonzero(diffab)[0]
            i_coords.extend([i] * len(nonzero_diffab_locs))
            j_coords.extend([j] * len(nonzero_diffab_locs))
            k_coords.extend(nonzero_diffab_locs)
            data_vals.extend(diffab[nonzero_diffab_locs])
    if make_diffab:
        logging.info(f"Converting diffabs to sparse matrix")
        coords = [i_coords, j_coords, k_coords]
        diffabs_sparse = sparse.COO(coords, data_vals, shape=diffab_dims)
        # add the negative transpose ([[0, 1], [0, 0]] -> [[0, 1], [-1, 0]])
        diffabs_sparse = diffabs_sparse - diffabs_sparse.transpose(axes=(1, 0, 2))

    logging.info(f"Saving results")
    # save the distances
    np.save(out_file, dists)
    logging.info(f"Saved distances to {out_file}")
    if make_diffab:
        logging.info("Saving diffabs as a sparse npz file")
        diffabs_out_file = out_file + '.diffab.npz'
        #np.save(out_file + '.diffab.npy', diffabs)
        sparse.save_npz(diffabs_out_file, diffabs_sparse)
        logging.info(f"Saved difference abundance diffabs to {diffabs_out_file}")
        # save the nodes in order, but using the original node names
        with open(diffabs_out_file + '.nodes.txt', 'w') as f:
            for node in nodes_in_order:
                f.write(f'{EMDU_index_2_node[node]}\n')
    # save the basis
    with open(out_file + '.basis.txt', 'w') as f:
        for file in fun_files:
            f.write(f"{file}\n")
    # if asked, save the pushed profiles
    if save_Ppushed:
        logging.info(f"Saving pushed profiles")
        ppushed_out_file = out_file + '.ppushed.pkl'
        with open(ppushed_out_file, 'wb') as f:
            pickle.dump(Ps_pushed, f)
        logging.info(f"Saved pushed profiles to {ppushed_out_file}")


if __name__ == '__main__':
    main()



