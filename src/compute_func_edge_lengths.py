#!/usr/bin/env python
import os
import datetime
import pandas as pd
import networkx as nx
import src.utility.kegg_db as kegg_db
import src.objects.func_tree as func_tree
import src.utility.kegg_process as kegg_process
import src.factory.make_pairwise_distance as make_pairwise_distance
import src.algorithms.edge_length_computation as edge_length_computation
import data
import src.utility.constant as constant
from pathlib import Path  


def main(args):
    brite = args.brite_id
    edge_file = args.edge_file
    out_dir = args.out_dir
    out_id = args.out_id
    distance_file = args.distance_file
    distance_label_file=f"{distance_file}.labels.txt"
    num_iter = int(args.num_iter)
    factor = int(args.factor)
    reg_factor = float(args.reg_factor)
    isdistance = args.distance
    ##############################################################################
    # Validations
    ##############################################################################
    if num_iter < 1:
        raise ValueError('Number of iterations must be at least 1')
    if factor < 1:
        raise ValueError('Factor must be at least 1')
    edge_file = data.get_data_abspath(edge_file)
    distance_file = data.get_data_abspath(distance_file)
    distance_label_file = data.get_data_abspath(distance_label_file)
    out_dir = data.get_data_abspath(out_dir, raise_if_not_found=False)
    if brite not in kegg_db.instance.brites:
        raise ValueError(f"{brite} is not a valid BRITE ID. Choices are: {kegg_db.instance.brites}")
    ##############################################################################
    # Main objects
    ##############################################################################
    edge_list = pd.read_csv(edge_file, sep='\t', header=0)
    G = nx.read_edgelist(edge_file, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
    tree = func_tree.FuncTree(G)
    tree.make_subtree(brite)
    pairwise_distances = make_pairwise_distance.get_KO_pairwise_dist(distance_file, distance_label_file)
    ##############################################################################
    # First, get A matrix following the process: "create_func_edge_matrix"
    ##############################################################################
    solver = edge_length_computation.EdgeLengthSolver()
    A = solver.get_A_matrix(tree, pairwise_distances)
    ##############################################################################
    # Compute edge length
    ##############################################################################
    edge_list = solver.compute_edges(A, tree.basis, edge_list, pairwise_distances, num_iter, factor, reg_factor, isdistance)
    ##############################################################################
    # Output
    ##############################################################################
    if out_id:
        out_file = f"{constant.EDGE_LIST_OUT__FILE_NAME.format(out_id)}"
    else:
        out_file = f"{constant.EDGE_LIST_OUT__FILE_NAME.format(datetime.datetime.now().strftime('%Y%m%d-%H%M%S'))}"
    filepath = Path(f"{out_dir}/{out_file}")  
    filepath.parent.mkdir(parents=True, exist_ok=True)  
    edge_list.to_csv(filepath, sep='\t', index=False)  


