#!/usr/bin/env python
import os
from os.path import exists
import networkx as nx
from scipy import sparse
import src.utility.kegg_db as kegg_db
import src.objects.func_tree as func_tree
import src.factory.make_leaf_distance as make_leaf_distance
import src.algorithms.edge_length_computation as edge_length_computation
import src.utility.constant as constant
import data


def main(args):
    edge_file           =args.edge_file
    distance_file       =args.distance_file
    distance_label_file =f"{args.distance_file}.labels.txt"
    out_dir             =args.out_dir
    brite_id            =args.brite_id
    ##############################################################################
    # Validations
    ##############################################################################
    # check if brite is legit
    if brite_id not in kegg_db.instance.brites:
        raise ValueError(f"{brite_id} is not a valid BRITE ID. Choices are: {kegg_db.instance.brites}")
    edge_file           =data.get_data_abspath(edge_file)
    distance_file       =data.get_data_abspath(distance_file)
    distance_label_file =data.get_data_abspath(distance_label_file)
    out_dir             =data.get_data_abspath(out_dir)
    ##############################################################################
    # Main Objects
    ##############################################################################
    G = nx.read_edgelist(edge_file, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
    tree = func_tree.FuncTree(G)
    tree.set_subtree(brite_id)
    pairwise_distances = make_leaf_distance.get_KO_pairwise_dist(distance_file, distance_label_file)
    ##############################################################################
    # Get A matrix
    ##############################################################################
    solver = edge_length_computation.EdgeLengthSolver()
    A = solver.get_A_matrix(tree, pairwise_distances)
    ##############################################################################
    # Output
    ##############################################################################
    if not exists(out_dir):
        os.mkdir(out_dir)
    # save the basis
    with open(f"{os.path.join(out_dir, constant.A_MATRIX_BASIS__FILE_NAME)}", 'w') as f:
        f.write('\n'.join(tree.basis))
    # save the sparse matrix
    print(f"Saving sparse matrix")
    sparse.save_npz(os.path.join(out_dir, constant.A_MATRIX__FILE_NAME), A)
