import src.factory.make_tree as make_tree
import src.factory.make_emd_input as make_emd_input
import src.utility.pyemd_simulation as pyemd_simulation
from src.algorithms.emd_unifrac import EarthMoverDistanceUniFracAbstract, EarthMoverDistanceUniFracSolver
import src.objects.emdu_diffab as emdu_diffab
import src.objects.func_tree as func_tree
import src.objects.emdu_vector as emdu_vector
import numpy as np
import pytest
import pandas as pd
import data


def test_emdu_vs_pyemd_simple():
    """
    Compare PyEMD with EMDUniFrac
    :return: None
    """
    solver: EarthMoverDistanceUniFracAbstract = EarthMoverDistanceUniFracSolver()
    test_edge_file = data.get_data_abspath('small_edge_list_with_lengths.txt')
    distance_matrix, node_list = pyemd_simulation.get_distance_matrix_from_edge_list(test_edge_file)
    # simple known distributions
    P = np.array([0, 0, 0, 0, 0, 0, 1])
    Q = np.array([0, 0, 0, 0, 0, 1, 0])
    # pyemd version
    pyemd_val = pyemd_simulation.get_EMD_pyemd(P, Q, distance_matrix)
    # EMDUniFrac version
    tree = make_tree.import_graph(test_edge_file, directed=True)
    Gdir = tree.main_tree
    input = make_emd_input.tree_to_EMDU_input(tree)
    Tint, lint, nodes_in_order, EMDU_index_2_node = input.Tint, input.lint, input.basis, input.idx_to_node
    # switch the order for convenience
    node_2_EMDU_index = {val: key for key, val in EMDU_index_2_node.items()}
    # convert from the distance matrix ordering to the EMDU ordering
    PU = np.zeros_like(P)
    QU = np.zeros_like(Q)
    for i, node in enumerate(node_list):
        PU[node_2_EMDU_index[node]] = P[i]
        QU[node_2_EMDU_index[node]] = Q[i]

    input = func_tree.EmdInput(Tint, lint, nodes_in_order, node_2_EMDU_index)
    # emdu_val1 = solver.solve_plain(input, PU, QU, weighted=True)
    emdu_val2, _ = solver.solve(input, PU, QU, weighted=True)
    emdu_val3, _, _ = solver.solve_with_flow(input, PU, QU, weighted=True)
    # assert np.isclose(pyemd_val, emdu_val1, atol=1e-8)
    assert np.isclose(pyemd_val, emdu_val2, atol=1e-8)
    assert np.isclose(pyemd_val, emdu_val3, atol=1e-8)


def test_emdu_vs_pyemd_random():
    """
    Compare PyEMD with EMDUniFrac with a bunch of random Ps and Qs
    :return: None
    """
    solver: EarthMoverDistanceUniFracAbstract = EarthMoverDistanceUniFracSolver()
    test_edge_file = data.get_data_abspath('small_edge_list_with_lengths.txt')
    distance_matrix, node_list = pyemd_simulation.get_distance_matrix_from_edge_list(test_edge_file)
    tree = make_tree.import_graph(test_edge_file, directed=True)
    input = make_emd_input.tree_to_EMDU_input(tree)
    Tint, lint, nodes_in_order, EMDU_index_2_node = input.Tint, input.lint, input.basis, input.idx_to_node
    # switch the order for convenience
    node_2_EMDU_index = {val: key for key, val in EMDU_index_2_node.items()}

    num_its = 10
    for i in range(num_its):
        # random dirichlet distributions
        # These are distributed over the nodes in the order of the distance matrix
        P = np.random.dirichlet(np.ones(len(node_list)))
        Q = np.random.dirichlet(np.ones(len(node_list)))
        # pyemd version
        pyemd_val = pyemd_simulation.get_EMD_pyemd(P, Q, distance_matrix)
        # The EMDU version needs the distributions in the order of the nodes_in_order list
        PU = np.zeros_like(P)
        QU = np.zeros_like(Q)
        for i, node in enumerate(node_list):
            PU[node_2_EMDU_index[node]] = P[i]
            QU[node_2_EMDU_index[node]] = Q[i]

        input = func_tree.EmdInput(Tint, lint, nodes_in_order, node_2_EMDU_index)
        # emdu_val1 = solver.solve_plain(input, PU, QU, weighted=True)
        emdu_val2, _ = solver.solve(input, PU, QU, weighted=True)
        emdu_val3, _, _ = solver.solve_with_flow(input, PU, QU, weighted=True)
        # assert np.isclose(pyemd_val, emdu_val1, atol=1e-8)
        assert np.isclose(pyemd_val, emdu_val2, atol=1e-8)
        assert np.isclose(pyemd_val, emdu_val3, atol=1e-8)


# Test functional profile conversion
def test_func_profile_convert():
    solver: EarthMoverDistanceUniFracAbstract = EarthMoverDistanceUniFracSolver()
    test_edge_file = data.get_data_abspath('small_edge_list_with_lengths.txt')
    functional_profile_file = data.get_data_abspath('small_sim_10_KOs_gather.csv')
    tree = make_tree.import_graph(test_edge_file, directed=True)
    input = make_emd_input.tree_to_EMDU_input(tree)
    functional_profile = pd.read_csv(functional_profile_file)
    P = make_emd_input.functional_profile_to_vector(functional_profile, input, abundance_key='median_abund', normalize=True)
    assert np.isclose(np.sum(P), 1.0, atol=1e-8)
    assert np.all(P >= 0)
    assert np.all(P <= 1)
    node_2_EMDU_index = {val: key for key, val in input.idx_to_node.items()}
    assert np.isclose(P[node_2_EMDU_index['e']], 2 / 6, atol=1e-8)
    assert np.isclose(P[node_2_EMDU_index['d']], 3 / 6, atol=1e-8)
    assert np.isclose(P[node_2_EMDU_index['f']], 1 / 6, atol=1e-8)
    # Try a bad abundance key
    with pytest.raises(ValueError):
        P = make_emd_input.functional_profile_to_vector(functional_profile, input, abundance_key='name', normalize=True)
    # Try a good, but different abundance_key
    P = make_emd_input.functional_profile_to_vector(functional_profile, input, abundance_key='f_orig_query', normalize=True)
    assert np.isclose(np.sum(P), 1.0, atol=1e-8)
    P = make_emd_input.functional_profile_to_vector(functional_profile, input, abundance_key='f_orig_query', normalize=False)
    assert np.isclose(np.sum(P), 1.258921, atol=1e-4)


def test_push_up_L1():
    solver: EarthMoverDistanceUniFracAbstract = EarthMoverDistanceUniFracSolver()
    test_edge_file = data.get_data_abspath('small_edge_list_with_lengths.txt')
    distance_matrix, node_list = pyemd_simulation.get_distance_matrix_from_edge_list(test_edge_file)
    tree = make_tree.import_graph(test_edge_file, directed=True)
    input = make_emd_input.tree_to_EMDU_input(tree)
    Tint, lint, nodes_in_order, EMDU_index_2_node = input.Tint, input.lint, input.basis, input.idx_to_node
    # switch the order for convenience
    node_2_EMDU_index = {val: key for key, val in EMDU_index_2_node.items()}

    num_its = 10
    for i in range(num_its):
        # random dirichlet distributions
        P = np.random.dirichlet(np.ones(len(node_list)))
        Q = np.random.dirichlet(np.ones(len(node_list)))
        # pyemd version
        pyemd_val = pyemd_simulation.get_EMD_pyemd(P, Q, distance_matrix)
        PU = np.zeros_like(P)
        QU = np.zeros_like(Q)
        for i, node in enumerate(node_list):
            PU[node_2_EMDU_index[node]] = P[i]
            QU[node_2_EMDU_index[node]] = Q[i]
        P_pushed = solver.push_up_L1(PU, input)
        Q_pushed = solver.push_up_L1(QU, input)
        pushed_emd = np.sum(np.abs(P_pushed - Q_pushed))
        assert np.isclose(pyemd_val, pushed_emd, atol=1e-5)


def test_diffab_indexer():
    solver: EarthMoverDistanceUniFracAbstract = EarthMoverDistanceUniFracSolver()
    edge_list_file = data.get_data_abspath("small_edge_list_with_lengths.txt")
    brite = "ko00001"
    file_pattern = "*_gather.csv"
    fun_files = data.get_data_abspaths(file_pattern)
    fun_files = sorted(fun_files)
    tree = make_tree.import_graph(edge_list_file, directed=True)
    input = make_emd_input.tree_to_EMDU_input(tree)
    Tint, lint, nodes_in_order, EMDU_index_2_node = input.Tint, input.lint, input.basis, input.idx_to_node
    # compute the diffabs
    profile_P = pd.read_csv(fun_files[0])
    P = make_emd_input.functional_profile_to_vector(profile_P, input, abundance_key="median_abund", normalize=True)
    P_pushed = solver.push_up_L1(P, input)
    profile_Q = pd.read_csv(fun_files[1])
    Q = make_emd_input.functional_profile_to_vector(profile_Q, input, abundance_key="median_abund", normalize=True)
    Q_pushed = solver.push_up_L1(Q, input)
    diffabs = np.zeros((len(fun_files), len(fun_files), len(nodes_in_order)))
    diffabs[0, 0, :] = emdu_vector.get_L1_diffab(P_pushed, P_pushed)
    diffabs[0, 1, :] = emdu_vector.get_L1_diffab(P_pushed, Q_pushed)
    diffabs[1, 0, :] = emdu_vector.get_L1_diffab(Q_pushed, P_pushed)
    diffabs[1, 1, :] = emdu_vector.get_L1_diffab(Q_pushed, Q_pushed)
    # instantiate the indexer
    indexer = emdu_diffab.DiffabArrayIndexer(diffabs, nodes_in_order, fun_files, EMDU_index_2_node)
    # test the indexer
    assert np.allclose(indexer.get_diffab(fun_files[0], fun_files[0]), np.zeros_like(nodes_in_order), atol=1e-8)
    assert np.allclose(indexer.get_diffab(fun_files[1], fun_files[1]), np.zeros_like(nodes_in_order), atol=1e-8)
    assert np.allclose(indexer.get_diffab(fun_files[1], fun_files[0]), diffabs[1, 0], atol=1e-8)
    assert np.allclose(indexer.get_diffab(fun_files[0], fun_files[1]), diffabs[0, 1], atol=1e-8)
    assert np.allclose(indexer.get_diffab(fun_files[0], fun_files), diffabs[[0], :, :], atol=1e-8)
    #print(f"first: {indexer.get_diffab(fun_files, fun_files[1])}")
    #print(f"second: {diffabs[:, 1, :]}")
    # Problem is that the first one is transposed-ish: I have an extra set of brackets between the two
    # solution is to enclose the index with brackets
    assert np.allclose(indexer.get_diffab(fun_files, fun_files[1]), diffabs[:, [1], :], atol=1e-8)
    assert np.allclose(indexer.get_diffab(fun_files, fun_files), diffabs, atol=1e-8)
    assert np.allclose(indexer.get_diffab_for_node(fun_files[0], fun_files[1], ['d', 'e', 'b', 'f', 'g', 'c', 'ko00001']), diffabs[0, 1, :], atol=1e-8)
    assert np.allclose(indexer.get_diffab_for_node(fun_files[0], fun_files[1], 'ko00001'), diffabs[0, 1, -1],
                      atol=1e-8)
    assert np.allclose(indexer.get_diffab_for_node(fun_files[0], fun_files, 'd'), diffabs[0][:, [0]],
                      atol=1e-8)
    assert np.allclose(indexer.get_diffab_for_node(fun_files[0], fun_files, ['ko00001', 'g']), diffabs[0][:, [-1, 4]],
                       atol=1e-8)


def test_EMDUnifrac_weighted_flow():
    solver: EarthMoverDistanceUniFracAbstract = EarthMoverDistanceUniFracSolver()
    path = data.get_data_abspath('small_edge_list_with_lengths_emdu.txt')
    tree = make_tree.import_graph(path)
    input = make_emd_input.tree_to_EMDU_input(tree)
    Tint, lint, nodes_in_order, EMDU_index_2_node = input.Tint, input.lint, input.basis, input.idx_to_node
    nodes_samples = {
        'C': {'sample1': 1, 'sample2': 0},
        'B': {'sample1': 1, 'sample2': 1},
        'A': {'sample1': 0, 'sample2': 0},
        'temp0': {'sample1': 0,
                  'sample2': 1}}  # temp0 is the root node, not named in Newick format, but included in nodes_in_order
    P = np.zeros(4)
    Q = np.zeros(4)
    for i in EMDU_index_2_node:
        P[i] = nodes_samples[EMDU_index_2_node[i]]['sample1']
        Q[i] = nodes_samples[EMDU_index_2_node[i]]['sample2']
    P = P/2
    Q = Q/2 #normalize
    (Z, F, diffab) = solver.solve_with_flow(input, P, Q, weighted=True)  # Run the weighted version of EMDUnifrac that returns the flow
    # Check to make sure results make sense
    print(EMDU_index_2_node)
    print(nodes_in_order)
    print(P, Q)
    print(F)
    assert Z == 0.25  # This is the Unifrac distance
    #assert F[(1, 1)] == 0.5
    #assert F[(0, 3)] == 0.5  # F is the flow and is in a sparse matrix format: a dictionary with tuple keys using elements of Tint and values T[(i, j)] equal to amount of abundance moved from organism nodes_in_order(i) to nodes_in_order(j)
    assert F[(0, 0)] == 0.5
    assert F[(1, 3)] == 0.5
    assert sum(F.values()) == 1
    #assert diffab == {(2, 3): 0.14999999999999999, (0, 2): 0.10000000000000001}  # diffab is the differential abundance vector, also in a sparse matrix format: a dictionary with tuple keys using elements of Tint and values diffab[(i, j)] equal to the signed difference of abundance between the two samples restricted to the sub-tree defined by nodes_in_order(i) and weighted by the edge length lint[(i,j)].
    assert np.isclose(diffab[(1, 2)], 0.1)
    assert np.isclose(diffab[(2, 3)], 0.15)


def test_EMDUnifrac_weighted():
    solver: EarthMoverDistanceUniFracAbstract = EarthMoverDistanceUniFracSolver()
    path = data.get_data_abspath('small_edge_list_with_lengths_emdu.txt')
    tree = make_tree.import_graph(path)
    input = make_emd_input.tree_to_EMDU_input(tree)
    Tint, lint, nodes_in_order, EMDU_index_2_node = input.Tint, input.lint, input.basis, input.idx_to_node
    G = tree.main_tree
    nodes_samples = {
        'C': {'sample1': 1, 'sample2': 0},
        'B': {'sample1': 1, 'sample2': 1},
        'A': {'sample1': 0, 'sample2': 0},
        'temp0': {'sample1': 0, 'sample2': 1}}  # temp0 is the root node
    P = np.zeros(4)
    Q = np.zeros(4)
    for i in EMDU_index_2_node:
        P[i] = nodes_samples[EMDU_index_2_node[i]]['sample1']
        Q[i] = nodes_samples[EMDU_index_2_node[i]]['sample2']
    P = P/2
    Q = Q/2 #normalize
    print(EMDU_index_2_node)
    print(nodes_in_order)
    print(P, Q)
    (Z, diffab) = solver.solve(input, P, Q, weighted=True)
    assert Z == 0.25
    print(diffab)
    #assert diffab == {(2, 3): 0.14999999999999999, (0, 2): 0.10000000000000001}
    assert np.isclose(diffab[(2,3)], 0.15)
    #assert np.isclose(diffab[(0,2)], 0.1)
    assert np.isclose(diffab[(1, 2)], 0.1)

def test_EMDUnifrac_unweighted():
    solver: EarthMoverDistanceUniFracAbstract = EarthMoverDistanceUniFracSolver()
    path = data.get_data_abspath('small_edge_list_with_lengths_emdu.txt')
    tree = make_tree.import_graph(path)
    input = make_emd_input.tree_to_EMDU_input(tree)
    Tint, lint, nodes_in_order, EMDU_index_2_node = input.Tint, input.lint, input.basis, input.idx_to_node
    G = tree.main_tree
    nodes_samples = {
        'C': {'sample1': 1, 'sample2': 0},
        'B': {'sample1': 1, 'sample2': 1},
        'A': {'sample1': 0, 'sample2': 0},
        'temp0': {'sample1': 0, 'sample2': 1}}
    P = np.zeros(4)
    Q = np.zeros(4)
    for i in EMDU_index_2_node:
        P[i] = nodes_samples[EMDU_index_2_node[i]]['sample1']
        Q[i] = nodes_samples[EMDU_index_2_node[i]]['sample2']
    (Z, diffab) = solver.solve(input, P, Q, weighted=False)
    assert Z == 0.5
    #assert diffab == {(2, 3): 0.29999999999999999, (0, 2): 0.20000000000000001}
    print(diffab)
    assert np.isclose(diffab[(2, 3)], 0.3)
    #assert np.isclose(diffab[(0, 2)], 0.2)
    assert np.isclose(diffab[(1, 2)], 0.2)



def test_EMDUnifrac_unweighted_flow():
    solver: EarthMoverDistanceUniFracAbstract = EarthMoverDistanceUniFracSolver()
    path = data.get_data_abspath('small_edge_list_with_lengths_emdu.txt')
    tree = make_tree.import_graph(path)
    input = make_emd_input.tree_to_EMDU_input(tree)
    Tint, lint, nodes_in_order, EMDU_index_2_node = input.Tint, input.lint, input.basis, input.idx_to_node
    G = tree.main_tree
    nodes_samples = {
        'C': {'sample1': 1, 'sample2': 0},
        'B': {'sample1': 1, 'sample2': 1},
        'A': {'sample1': 0, 'sample2': 0},
        'temp0': {'sample1': 0, 'sample2': 1}}
    P = np.zeros(4)
    Q = np.zeros(4)
    for i in EMDU_index_2_node:
        P[i] = nodes_samples[EMDU_index_2_node[i]]['sample1']
        Q[i] = nodes_samples[EMDU_index_2_node[i]]['sample2']
    (Z, F, diffab) = solver.solve_with_flow(input, P, Q, weighted=False)
    print(F)
    print(diffab)
    assert Z == 0.5
    #assert F[(0, 3)] == 1
    #assert F[(1, 1)] == 1
    assert F[(1, 3)] == 1.
    assert F[(0, 0)] == 1.
    assert sum(F.values()) == 2.
    #assert diffab == {(2, 3): 0.29999999999999999, (0, 2): 0.20000000000000001}
    assert np.isclose(diffab[(2, 3)], 0.3)
    #assert np.isclose(diffab[(0, 2)], 0.2)
    assert np.isclose(diffab[(1, 2)], 0.2)

