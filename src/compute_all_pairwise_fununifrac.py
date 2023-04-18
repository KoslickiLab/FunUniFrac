#!/usr/bin/env python
# Script to compute all pairwise functional unifrac from a directory of functional profiles
import os, glob
import numpy as np
from scipy import sparse
from src.algorithms.emd_unifrac import EarthMoverDistanceUniFracSolver
import src.utility.kegg_db as kegg_db
import logging
import sparse
import pickle
import data
from src.objects.func_tree import FuncTreeEmduInput
import src.factory.make_tree as make_tree
import src.factory.make_emd_input as make_emd_input
import src.utility.constant as constant
import datetime
from pathlib import Path  


def main(args):
    logging.basicConfig(level=args.loglevel, format='%(asctime)s %(levelname)s: %(message)s')
    edge_list_file = args.edge_list
    out_dir = args.out_dir
    out_id = args.out_id
    file_dir = args.file_dir
    file_pattern = args.file_pattern
    abundance_key = args.abundance_key
    num_threads = args.threads
    brite = args.brite
    make_diffab = args.diffab
    unweighted = args.unweighted
    is_L2 = args.L2
    save_Ppushed = args.Ppushed
    ##############################################################################
    # Validations
    ##############################################################################
    if brite not in kegg_db.instance.brites:
        raise ValueError(f"{brite} is not a valid BRITE ID. Choices are: {kegg_db.instance.brites}")
    data.check_data_abspath(edge_list_file)
    # if data.check_data_abspath(out_file, raise_if_not_found=False) and not force:
    #     raise FileExistsError(f"{out_file} already exists. Please delete it or choose another name, or use --force.")
    # look for files in that directory with that file pattern
    file_pattern = os.path.join(file_dir, file_pattern)
    data.check_data_abspaths(file_pattern)
    fun_files = glob.glob(file_pattern)
    fun_files = sorted(fun_files)
    logging.info(f"Parsing graph")
    ##############################################################################
    # Main objects
    ##############################################################################
    solver = EarthMoverDistanceUniFracSolver() 
    # Parse the graph and get the FunUniFrac objects Tint, length, and edge_list
    tree = make_tree.import_graph(edge_list_file, directed=True)
    logging.info(f"Converting graph into EMDUniFrac format")
    # Then create the inputs for EMDUniFrac
    input: FuncTreeEmduInput = make_emd_input.tree_to_EMDU_input(tree, brite)
    logging.info(f"Computing pairwise distances in parallel")
    ##############################################################################
    # Computations
    ##############################################################################
    # convert the functional profiles to vectors
    results = make_emd_input.parallel_functional_profile_to_vector(num_threads, fun_files, input, abundance_key, unweighted)
    # push vectors up in parallel
    Ps_pushed = {}
    for file, P in results:
        if is_L2:
            P_pushed = solver.push_up_L2(P, input)
        else:
            P_pushed = solver.push_up_L1(P, input)
        Ps_pushed[file] = P_pushed
    logging.info(f"Computing pairwise distances")
    # Then compute the pairwise distances
    dists, diffabs_sparse = solver.pairwise_computation(Ps_pushed, fun_files, input, make_diffab, is_L2)
    ##############################################################################
    # Save Outputs
    ##############################################################################
    logging.info(f"Saving results")
    # save the distances
    fid = out_id if out_id else datetime.datetime.now().strftime('%Y%m%d-%H%M%S')
    out_file = Path(out_dir, constant.FUNUNIFRAC_OUT__MAIN_FILE_NAME.format(fid))
    out_file.parent.mkdir(parents=True, exist_ok=True)  

    np.save(out_file, dists)
    # save the basis
    out_basis_file = os.path.join(out_dir, constant.FUNUNIFRAC_OUT__MAIN_BASIS__FILE_NAME.format(fid))
    with open(out_basis_file, 'w') as f:
        for file in fun_files:
            f.write(f"{file}\n")
    logging.info(f"Saved distances to {out_file}")
    if make_diffab:
        logging.info("Saving diffabs as a sparse npz file")
        out_diffab_file = os.path.join(out_dir, constant.FUNUNIFRAC_OUT__DIFFAB_FILE_NAME.format(fid))
        #np.save(out_file + '.diffab.npy', diffabs)
        sparse.save_npz(out_diffab_file, diffabs_sparse)
        logging.info(f"Saved difference abundance diffabs to {out_file}")
        # save the nodes in order, but using the original node names
        out_diffab_node_file = os.path.join(out_dir, constant.FUNUNIFRAC_OUT__DIFFAB_NODES_FILE_NAME.format(fid))
        with open(out_diffab_node_file, 'w') as f:
            for node in input.basis:
                f.write(f'{input.idx_to_node[node]}\n')
    
    # if asked, save the pushed profiles
    if save_Ppushed:
        logging.info(f"Saving pushed profiles")
        out_pushed_file = os.path.join(out_dir, constant.FUNUNIFRAC_OUT__PUSH_FILE_NAME.format(fid))
        with open(out_pushed_file, 'wb') as f:
            pickle.dump(Ps_pushed, f)
        logging.info(f"Saved pushed profiles to {out_pushed_file}")
