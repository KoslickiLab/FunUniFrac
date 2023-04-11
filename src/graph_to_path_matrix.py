#!/usr/bin/env python
import os
from os.path import exists
import networkx as nx
from scipy import sparse
import src.utils.kegg_db as kegg_db
import src.base.func_tree as func_tree
import src.algorithms.kegg_process as kegg_process
import src.algorithms.edge_length_computation as edge_length_computation
import src.base.constant as constant
import data



def main(args):
    edge_file       =data.get_data_abspath(args.edge_file)
    distance_file   =data.get_data_abspath(args.distance_file)
    out_dir         =data.get_data_abspath(args.out_dir)
    brite_id        =args.brite_id

    # check if brite is legit
    if brite_id not in kegg_db.instance.brites:
        raise ValueError(f"{brite_id} is not a valid BRITE ID. Choices are: {kegg_db.instance.brites}")

    KO_dist_labels = []
    with open(distance_file, 'r') as f:
        for line in f.readlines():
            ko = line.strip().split('|')[-1]  # KO number is the last in the pip-delim list
            ko = ko.split(':')[-1]  # remove the ko: prefix
            KO_dist_labels.append(ko)
    G = nx.read_edgelist(edge_file, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
    tree = func_tree.FuncTree(G)
    tree.apply_classification(brite_id)
    pairwise_distances, _ = kegg_process.get_KO_indices(KO_dist_labels, basis=tree.basis)

    solver = edge_length_computation.EdgeLengthSolver()
    basis, A = solver.get_A_matrix(tree, pairwise_distances)

    if not exists(out_dir):
        os.mkdir(out_dir)
    # Save the basis
    with open(f"{os.path.join(out_dir, constant.A_MATRIX_BASIS__FILE_NAME)}", 'w') as f:
        f.write('\n'.join(basis))
    # save the sparse matrix
    print(f"Saving sparse matrix")
    sparse.save_npz(os.path.join(out_dir, constant.A_MATRIX__FILE_NAME), A)
