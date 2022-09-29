# unit tests for LP_EMD_helper.py
import networkx as nx
import numpy as np

from src.LP_EMD_helper import get_EMD_pyemd,\
get_EMDUniFrac_from_functional_profiles,\
get_EMD_from_edge_file,\
get_ID_index_dict,\
get_distance_matrix_from_edge_list,\
get_leaf_nodes,\
get_leaf_nodes_only_graph,\
get_distance_matrix_from_edge_list,\
make_edge_list_file_len_1_tmp,\
parse_edge_list,\
simulate_leaf_supported_vector,\
get_graphs_and_index,\
LeafDistributionSimulator



def test_get_distance_matrix_from_edge_list():
    test_edge_file ='test_data/small_edge_list_with_lengths.txt'
    distance_matrix, node_list = get_distance_matrix_from_edge_list(test_edge_file)
    # test that the distance matrix is the right size
    assert distance_matrix.shape == (len(node_list), len(node_list))
    # test against precomputed distance matrix
    # distance matrix for binary with 4 leaves and edge lengths of 1
    known_D = np.array([[0, 1, 1, 2, 2, 2, 2],
                        [1, 0, 2, 1, 1, 3, 3],
                        [1, 2, 0, 3, 3, 1, 1],
                        [2, 1, 3, 0, 2, 4, 4],
                        [2, 1, 3, 2, 0, 4, 4],
                        [2, 3, 1, 4, 4, 0, 2],
                        [2, 3, 1, 4, 4, 2, 0]])
    assert np.allclose(distance_matrix, known_D, atol=1e-2)


def test_get_EMD():
    D = np.array([[0, 1, 1, 2, 2, 2, 2],
                  [1, 0, 2, 1, 1, 3, 3],
                  [1, 2, 0, 3, 3, 1, 1],
                  [2, 1, 3, 0, 2, 4, 4],
                  [2, 1, 3, 2, 0, 4, 4],
                  [2, 3, 1, 4, 4, 0, 2],
                  [2, 3, 1, 4, 4, 2, 0]])
    P = np.array([0, 0, 0, 0, 0, 0, 1])
    Q = np.array([0, 0, 0, 0, 0, 1, 0])
    emd = get_EMD_pyemd(P, Q, D)
    assert np.isclose(emd, 2, atol=1e-2)
    P = np.array([0, 0, 0, 0, 0, 0, 1])
    Q = np.array([0, 0, 0, 0, 1, 0, 0])
    emd = get_EMD_pyemd(P, Q, D)
    assert np.isclose(emd, 4, atol=1e-2)
    P = np.array([0, 0, 0, 0, 0, 0, 2])
    Q = np.array([0, 0, 0, 0, 0, 2, 0])
    # the following should fail since P and Q are not normalized
    try:
        emd = get_EMD_pyemd(P, Q, D)
        assert False
    except:
        assert True
    P = np.array([1, 0, 0, 0, 0, 0, 0])
    Q = np.array([0, 0, 0, 0, 0, 1, 0])
    emd = get_EMD_pyemd(P, Q, D)
    assert np.isclose(emd, 2, atol=1e-2)

def test_get_graphs_and_index():
    test_edge_file ='test_data/small_edge_list.txt'
    Gdir, Gundir, basis, index = get_graphs_and_index(test_edge_file, "ko00001")
    assert len(Gdir.edges()) == 6
    assert len(Gundir.edges()) == 6
    assert len(basis) == 7
    assert len(index) == 7
    assert index["ko00001"] == 0
    assert nx.is_directed(Gdir)
    assert not nx.is_directed(Gundir)


def test_leaf_node_simulator():
    edge_list_file = 'test_data/small_edge_list.txt'
    brite = 'ko00001'
    l = LeafDistributionSimulator(edge_list_file, brite)
    assert len(l.basis) == 7
    P = l.get_random_dist_on_leaves()
    # check that the simulated distribution is a valid probability distribution
    assert np.allclose(np.sum(P), 1, atol=1e-2)
    assert np.all(P >= 0)
    assert np.all(P <= 1)
    # check if the distribution is supported only on the leaves
    for node in l.basis:
        if node not in l.leaf_nodes:
            assert np.isclose(P[l.basis_index[node]], 0, atol=1e-8)
        else:
            assert P[l.basis_index[node]] > 0