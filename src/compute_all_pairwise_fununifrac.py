#!/usr/bin/env python
# Script to compute all pairwise functional unifrac from a directory of functional profiles
import os
import numpy as np
import pandas as pd
from scipy import sparse
from src.algorithms.emd_unifrac import EarthMoverDistanceUniFracAbstract, EarthMoverDistanceUniFracSolver
import multiprocessing
from itertools import combinations, repeat
import src.utility.kegg_db as kegg_db
import logging
from blist import blist
import sparse
import pickle
import data
from src.objects.func_tree import FuncTreeEmduInput
from src.objects.profile_vector import get_L2_diffab, get_L2, get_L1, get_L1_diffab
import src.factory.make_tree as make_tree
import src.factory.make_emd_input as make_emd_input

solver: EarthMoverDistanceUniFracAbstract = EarthMoverDistanceUniFracSolver() 


def get_func_profiles_parallel(args):
    def map_func(file, input: FuncTreeEmduInput, abundance_key, unweighted, is_L2):
        df = pd.read_csv(file)
        P = make_emd_input.functional_profile_to_vector(df, input, abundance_key=abundance_key, normalize=True)
        if unweighted:
            # if entries are positive, set to 1
            P[P > 0] = 1
        if is_L2:
            P_pushed = solver.push_up_L2(P, input)
        else:
            P_pushed = solver.push_up_L1(P, input)
        return file, P_pushed
    return map_func(*args)


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
    if brite not in kegg_db.instance.brites:
        raise ValueError(f"{brite} is not a valid BRITE ID. Choices are: {kegg_db.instance.brites}")
    edge_list_file = data.get_data_abspath(edge_list_file)
    out_file = data.get_data_abspath(out_file, raise_if_not_found=False)
    if os.path.exists(out_file) and not force:
        raise FileExistsError(f"{out_file} already exists. Please delete it or choose another name, or use --force.")
    # look for files in that directory with that file pattern
    fun_files = data.get_data_abspaths(file_pattern)
    fun_files = sorted(fun_files)
    logging.info(f"Parsing graph")
    # Parse the graph and get the FunUniFrac objects Tint, length, and edge_list
    tree = make_tree.import_graph(edge_list_file, directed=True)

    logging.info(f"Converting graph into EMDUniFrac format")
    # Then create the inputs for EMDUniFrac
    input: FuncTreeEmduInput = make_emd_input.tree_to_EMDU_input(tree, brite)

    logging.info(f"Computing pairwise distances in parallel")
    # convert the functional profiles to vectors and push them up in parallel
    Ps_pushed = {}
    pool = multiprocessing.Pool(args.threads)
    results = pool.imap(get_func_profiles_parallel, zip(fun_files, repeat(input), repeat(abundance_key), repeat(unweighted), repeat(is_L2)),
                        chunksize=max(2, len(fun_files) // num_threads))
    pool.close()
    pool.join()
    
    for file, P_pushed in results:
        Ps_pushed[file] = P_pushed

    logging.info(f"Computing pairwise distances")
    # Then compute the pairwise distances
    dists = np.zeros((len(fun_files), len(fun_files)))
    i_coords = blist([])
    j_coords = blist([])
    k_coords = blist([])
    data_vals = blist([])
    diffab_dims = (len(fun_files), len(fun_files), len(input.basis))
    for i, j in combinations(range(len(fun_files)), 2):
        if not make_diffab:
            if not is_L2:
                dists[i, j] = dists[j, i] = get_L1(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
            else:
                dists[i, j] = dists[j, i] = get_L2(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
        else:
            if not is_L2:
                Z = get_L1(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
                diffab = get_L1_diffab(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
            else:
                Z = get_L2(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
                diffab = get_L2_diffab(Ps_pushed[fun_files[i]], Ps_pushed[fun_files[j]])
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
            for node in input.basis:
                f.write(f'{input.idx_to_node[node]}\n')
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



