import src.LP_EMD_helper as LH
import src.EMDU as EMDU
import numpy as np
import pytest

def test_emdu_vs_pyemd_simple():
    """
    Compare PyEMD with EMDUniFrac
    :return: None
    """
    test_edge_file = 'test_data/small_edge_list_with_lengths.txt'
    distance_matrix, node_list = LH.get_distance_matrix_from_edge_list(test_edge_file)
    # simple known distributions
    P = np.array([0, 0, 0, 0, 0, 0, 1])
    Q = np.array([0, 0, 0, 0, 0, 1, 0])
    # pyemd version
    pyemd_val = LH.get_EMD_pyemd(P, Q, distance_matrix)
    # EMDUniFrac version
    Gdir = LH.import_graph(test_edge_file, directed=True)
    Tint, lint, nodes_in_order, EMDU_index_2_node = LH.weighted_tree_to_EMDU_input(Gdir)
    # switch the order for convenience
    node_2_EMDU_index = {val: key for key, val in EMDU_index_2_node.items()}
    # convert from the distance matrix ordering to the EMDU ordering
    PU = np.zeros_like(P)
    QU = np.zeros_like(Q)
    for i, node in enumerate(node_list):
        PU[node_2_EMDU_index[node]] = P[i]
        QU[node_2_EMDU_index[node]] = Q[i]
    emdu_val1 = EMDU.EMDUnifrac_weighted_plain(Tint, lint, nodes_in_order, PU, QU)
    emdu_val2, _ = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, PU, QU)
    emdu_val3, _, _ = EMDU.EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, PU, QU)
    assert np.isclose(pyemd_val, emdu_val1, atol=1e-8)
    assert np.isclose(pyemd_val, emdu_val2, atol=1e-8)
    assert np.isclose(pyemd_val, emdu_val3, atol=1e-8)


def test_emdu_vs_pyemd_random():
    """
    Compare PyEMD with EMDUniFrac with a bunch of random Ps and Qs
    :return: None
    """
    test_edge_file = 'test_data/small_edge_list_with_lengths.txt'
    distance_matrix, node_list = LH.get_distance_matrix_from_edge_list(test_edge_file)
    Gdir = LH.import_graph(test_edge_file, directed=True)
    Tint, lint, nodes_in_order, EMDU_index_2_node = LH.weighted_tree_to_EMDU_input(Gdir)
    # switch the order for convenience
    node_2_EMDU_index = {val: key for key, val in EMDU_index_2_node.items()}

    num_its = 10
    for i in range(num_its):
        # random dirichlet distributions
        P = np.random.dirichlet(np.ones(len(node_list)))
        Q = np.random.dirichlet(np.ones(len(node_list)))
        # pyemd version
        pyemd_val = LH.get_EMD_pyemd(P, Q, distance_matrix)
        PU = np.zeros_like(P)
        QU = np.zeros_like(Q)
        for i, node in enumerate(node_list):
            PU[node_2_EMDU_index[node]] = P[i]
            QU[node_2_EMDU_index[node]] = Q[i]
        emdu_val1 = EMDU.EMDUnifrac_weighted_plain(Tint, lint, nodes_in_order, PU, QU)
        emdu_val2, _ = EMDU.EMDUnifrac_weighted(Tint, lint, nodes_in_order, PU, QU)
        emdu_val3, _, _ = EMDU.EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, PU, QU)
        assert np.isclose(pyemd_val, emdu_val1, atol=1e-8)
        assert np.isclose(pyemd_val, emdu_val2, atol=1e-8)
        assert np.isclose(pyemd_val, emdu_val3, atol=1e-8)


# Test functional profile conversion
def test_func_profile_convert():
    test_edge_file = 'test_data/small_edge_list_with_lengths.txt'
    functional_profile_file = 'test_data/small_sim_10_KOs_gather.csv'
    Gdir = LH.import_graph(test_edge_file, directed=True)
    Tint, lint, nodes_in_order, EMDU_index_2_node = LH.weighted_tree_to_EMDU_input(Gdir)
    node_2_EMDU_index = {val: key for key, val in EMDU_index_2_node.items()}
    P = EMDU.functional_profile_to_EMDU_vector(functional_profile_file, EMDU_index_2_node,
                                               abundance_key='median_abund', normalize=True)
    assert np.isclose(np.sum(P), 1.0, atol=1e-8)
    assert np.all(P >= 0)
    assert np.all(P <= 1)
    assert np.isclose(P[node_2_EMDU_index['e']], 2 / 6, atol=1e-8)
    assert np.isclose(P[node_2_EMDU_index['d']], 3 / 6, atol=1e-8)
    assert np.isclose(P[node_2_EMDU_index['f']], 1 / 6, atol=1e-8)
    # Try a bad abundance key
    with pytest.raises(ValueError):
        P = EMDU.functional_profile_to_EMDU_vector(functional_profile_file, EMDU_index_2_node,
                                                   abundance_key='name', normalize=True)
    # Try a good, but different abundance_key
    P = EMDU.functional_profile_to_EMDU_vector(functional_profile_file, EMDU_index_2_node,
                                               abundance_key='f_orig_query', normalize=True)
    assert np.isclose(np.sum(P), 1.0, atol=1e-8)
    P = EMDU.functional_profile_to_EMDU_vector(functional_profile_file, EMDU_index_2_node,
                                               abundance_key='f_orig_query', normalize=False)
    assert np.isclose(np.sum(P), 1.258921, atol=1e-4)