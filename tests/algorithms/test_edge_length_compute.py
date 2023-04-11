import os
import data
import subprocess
import numpy as np
from scipy import sparse
import networkx as nx
from src.algorithms.edge_length_computation import EdgeLengthSolver
from src.utils.func_tree import FuncTree
from src.algorithms.kegg_process import get_KO_indices


def test_small_edge_list():
    """
    Uses a complete binary tree with 4 leaf nodes, all branch lengths set to 1
    :return: None
    """
    dist_name = "small_pairwise_distances.npy"
    edge_list = data.get_data_abspath("small_edge_list.txt")
    brite = "ko00001"
    distances_labels_file = data.get_data_abspath(f"sourmash/{dist_name}.labels.txt")
    KO_dist_labels = []
    with open(distances_labels_file, 'r') as f:
        for line in f.readlines():
            ko = line.strip().split('|')[-1]  # KO number is the last in the pip-delim list
            ko = ko.split(':')[-1]  # remove the ko: prefix
            KO_dist_labels.append(ko)
    G = nx.read_edgelist(edge_list, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
    tree = FuncTree(G)
    tree.apply_classification(brite)
    pairwise_distances, _ = get_KO_indices(KO_dist_labels, basis=tree.basis)

    solver = EdgeLengthSolver()
    _, A = solver.get_A_matrix(tree, pairwise_distances)

    assert A.shape == (4**2, 6+1)
    # use the hand-calculated A matrix to check that the output is correct
    A_correct = np.array([[0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, 1, 0, 0],
                          [0, 1, 1, 1, 0, 1, 0],
                          [0, 1, 1, 1, 0, 0, 1],
                          [0, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0],
                          [0, 1, 1, 0, 1, 1, 0],
                          [0, 1, 1, 0, 1, 0, 1],
                          [0, 1, 1, 1, 0, 1, 0],
                          [0, 1, 1, 0, 1, 1, 0],
                          [0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 1, 1],
                          [0, 1, 1, 1, 0, 0, 1],
                          [0, 1, 1, 0, 1, 0, 1],
                          [0, 0, 0, 0, 0, 1, 1],
                          [0, 0, 0, 0, 0, 0, 0]])
    assert np.allclose(A.toarray(), A_correct)
