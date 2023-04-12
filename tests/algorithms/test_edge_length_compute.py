import data
import numpy as np
from scipy import sparse
import networkx as nx
import pandas as pd
from src.algorithms.edge_length_computation import EdgeLengthSolver
from src.base.func_tree import FuncTree
from src.base.pairwise_dist import PairwiseDistance
from src.algorithms.kegg_process import get_KO_indices


def test__create_A_matrix__with_small_edge_list():
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



def test__edge_length_computation__with_small_edge_lengths():
    """
    Uses a complete binary tree with 4 leaf nodes, all branch lengths set to 1
    :return: None
    """
    edge_list_file = data.get_data_abspath("small_edge_list.txt")
    distances_file = data.get_data_abspath(f"sourmash/small_pairwise_distances.npy")
    distances_labels_file = data.get_data_abspath(f"sourmash/small_pairwise_distances.npy.labels.txt")
    brite = "ko00001"
    A_file = data.get_data_abspath(f"{brite}_small_pairwise_distances.npy_A.npz")
    basis_name = data.get_data_abspath(f"{brite}_small_pairwise_distances.npy_column_basis.txt")
    
    edge_list = pd.read_csv(edge_list_file, sep='\t', header=0)

    with open(basis_name, 'r') as f:
        basis = f.readlines()
        basis = [line.strip() for line in basis]

    pw_dist = np.load(distances_file)
    KO_dist_labels = []
    with open(distances_labels_file, 'r') as f:
        for line in f.readlines():
            ko = line.strip().split('|')[-1]  # KO number is the last in the pip-delim list
            ko = ko.split(':')[-1]  # remove the ko: prefix
            KO_dist_labels.append(ko)
    pw_dist_labels, pw_dist_indices = get_KO_indices(KO_dist_labels, basis=basis)
    pairwise_distances = PairwiseDistance(pw_dist, pw_dist_labels, pw_dist_indices)

    A = sparse.load_npz(A_file)    

    solver = EdgeLengthSolver()
    df = solver.compute_edges(A=A, 
                              basis=basis, 
                              pairwise_distances=pairwise_distances, 
                              edge_list=edge_list, 
                              num_iter=10, 
                              factor=2, 
                              reg_factor=100, 
                              isdistance=True)

    # get the edge lengths
    edge_lengths = df["edge_length"].values
    # check if the edge lengths are all close to 1
    assert np.allclose(edge_lengths, np.ones_like(edge_lengths), atol=1e-2)
