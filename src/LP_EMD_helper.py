import itertools
import os.path

from pyemd import emd, emd_with_flow
import numpy as np
import pandas as pd
import time
import networkx as nx
import multiprocessing
from .CONSTANTS import BRITES


# FIXME: this says to delete, but it is used below
###delete
def parse_edge_list(file):
    df = pd.read_table(file)
    return df


def shortest_path_map_func_star(a_b_c):
    """Convert `f([1,2,3])` to `f(1,2,3)` call."""
    return shortest_path_map_func(*a_b_c)


def shortest_path_map_func(node, Gundir, edge_len_property):
    # compute the distance from the node to all leaf nodes
    len_dict = dict(nx.single_source_dijkstra_path_length(Gundir, node, weight=edge_len_property))
    return {node: len_dict}


def import_graph(edge_list_file, directed=True):
    """
    Import a graph from an edge list. The edge list can have 2 or 3 columns corresponding to parent, child,
    and (optionally) edge length.

    :param edge_list_file: text file with edge list, tab separated
    :param directed: boolean, whether the graph is directed or not
    :return: DiGraph
    """
    if not os.path.exists(edge_list_file):
        raise ValueError(f"Could not find edge list file {edge_list_file}")
    try:
        if directed:
            G = nx.read_edgelist(edge_list_file, delimiter='\t', nodetype=str, create_using=nx.DiGraph)
        else:
            G = nx.read_edgelist(edge_list_file, delimiter='\t', nodetype=str, create_using=nx.Graph)
    except:
        try:
            with open(edge_list_file, 'rb') as fid:
                # read the first line
                line = fid.readline()
                try:
                    line = line.decode('utf-8')
                    col1, col2, col3 = line.strip().split('\t')
                except:
                    raise ValueError('The first line of the file must have at most three columns: parent child '
                                     'edge_length')
                # import the graph. First two columns are the nodes, last column is the weight
                if directed:
                    G = nx.read_edgelist(fid, delimiter='\t', data=((col3, float),), create_using=nx.DiGraph)
                else:
                    G = nx.read_edgelist(fid, delimiter='\t', data=((col3, float),), create_using=nx.Graph)
        except:
            raise Exception(
                'Could not read edge list file. Make sure it is a tab-delimited file with two columns: parent & child'
                'or three columns: parent, child, edge_length')
    return G


def infer_edge_len_property(G):
    """
    From a networkx graph, find out the name of the edge length property

    :param G: graph
    :return: name of the edge length property
    """
    edge_properties = list(list(G.edges(data=True))[0][-1].keys())
    if len(edge_properties) > 1:
        raise ValueError(f'I found multiple edge properties: {edge_properties}. I don\'t know which one to use for '
                         f'edge lengths.')
    else:
        edge_len_property = edge_properties[0]
    return edge_len_property


def weighted_tree_to_EMDU_input(G: nx.DiGraph, edge_len_property=None):
    """
    Convert a weighted graph to the format required by EMDUniFrac: Tint, lint, and nodes_in_order.
    Since EMDUniFrac wants everything to be integers, also return a mapping from integers to node names.

    :param G: weighted networkx graph
    :return: (Tint dict, lint dict, nodes_in_order array, index_2_node dict)
    """
    if type(G) is not nx.DiGraph:
        raise ValueError('The graph must be a directed graph')
    # first determine the edge_length attribute
    if edge_len_property is None:
        edge_len_property = infer_edge_len_property(G)
    Tint = {}
    nodes = G.nodes()
    for node in nodes:
        pred = list(G.predecessors(node))
        if len(pred) > 1:
            raise Exception(f'Node {node} has more than one parent: {pred}')
        elif len(pred) == 1:
            Tint[node] = pred[0]
        else:
            pass
            #Tint[node] = "-1"  # TODO: check if this is how we do it in EMDUniFrac
    lint = {}
    for edge_data in G.edges(data=True):
        i, j, data = edge_data
        weight = data[edge_len_property]
        lint[(i, j)] = weight
        lint[(j, i)] = weight
    root = get_root_of_tree(G)
    nodes_in_order = list(nx.dfs_postorder_nodes(G, source=root))
    # Tint, lint, and nodes_in_order require these to be integers, so let's map everything to that
    node_2_EMDU_index = {node: i for i, node in enumerate(nodes_in_order)}
    # convert the Tint dict
    Tint = {node_2_EMDU_index[node]: node_2_EMDU_index[Tint[node]] for node in Tint}
    # convert the lint dict
    lint = {(node_2_EMDU_index[i], node_2_EMDU_index[j]): lint[(i, j)] for i, j in lint}
    # convert the nodes in order
    nodes_in_order = [node_2_EMDU_index[node] for node in nodes_in_order]
    # convert the integers to nodes
    EMDU_index_2_node = {i: node for node, i in node_2_EMDU_index.items()}
    return Tint, lint, nodes_in_order, EMDU_index_2_node


def get_distance_matrix_from_edge_list(edge_list_file, edge_len_property=None):
    '''
    Given an edge list file, with three columns: head, tail, length, return an ALL-PAIRS distance matrix and a list of
    nodes

    :param edge_list_file: file name of the edge list
    :param edege_length_property: name of the edge length property
    :return: distance matrix and a list of nodes indexing the distance matrix
    '''
    G = import_graph(edge_list_file, directed=False)
    if edge_len_property is None:
        edge_len_property = infer_edge_len_property(G)
    D_dict = dict(nx.all_pairs_dijkstra_path_length(G, weight=edge_len_property))
    node_list = list(G.nodes())
    distance_matrix = np.zeros((len(node_list), len(node_list)))
    for i, node1 in enumerate(node_list):
        for j, node2 in enumerate(node_list):
            distance_matrix[i, j] = D_dict[node1][node2]
    return distance_matrix, node_list


def get_distance_matrix_on_leaves_from_edge_list(edge_list_file, edge_len_property=None):
    '''
    Given an edge list file, with three columns: head, tail, length, return a distance matrix on all leaves and a
    list of nodes

    :param edge_list_file: file name of the edge list
    :param edege_length_property: name of the edge length property
    :return: distance matrix and a list of nodes indexing the distance matrix
    '''
    print('Reading edge list file...')
    Gdir = import_graph(edge_list_file, directed=True)
    Gundir = Gdir.to_undirected()
    # check if the graph is connected
    print('Checking if the graph is connected...')
    if not nx.is_connected(Gundir):
        raise ValueError('The graph is not connected. I cannot compute the distance matrix on leaves.')
    if edge_len_property is None:
        edge_properties = list(list(Gdir.edges(data=True))[0][-1].keys())
        if len(edge_properties) > 1:
            raise ValueError(f'I found multiple edge properties: {edge_properties}. I don\'t know which one to use for '
                             f'edge lengths.')
        else:
            edge_len_property = edge_properties[0]
    # get the leaf nodes
    print('Getting the leaf nodes...')
    leaf_nodes = get_leaf_nodes(Gdir)
    leaf_nodes_to_index = {n: i for i, n in enumerate(leaf_nodes)}
    # for each pair of leaf nodes, get the distance
    print('Initializing the distance matrix...')
    distance_matrix = np.zeros((len(leaf_nodes), len(leaf_nodes)))
    itr = 0
    len_leaves = len(leaf_nodes)
    print('Computing the distance matrix...')
    #TODO: specify the number of processes to use
    pool = multiprocessing.Pool(processes=int(multiprocessing.cpu_count() / 2))
    len_dicts = pool.map(shortest_path_map_func_star, zip(leaf_nodes, itertools.repeat(Gundir), itertools.repeat(edge_len_property)))
    pool.close()
    pool.join()
    # join up the results
    len_dict = {k:v for x in len_dicts for k,v in x.items()}
    for leaf1 in leaf_nodes:
        for leaf2 in leaf_nodes:
            i = leaf_nodes_to_index[leaf1]
            j = leaf_nodes_to_index[leaf2]
            distance_matrix[i, j] = len_dict[leaf1][leaf2]
            distance_matrix[j, i] = len_dict[leaf1][leaf2]
    return distance_matrix, leaf_nodes


def get_EMD_pyemd(P, Q, distance_matrix, with_flow=False):
    """
    Given two vectors P and Q, and a distance matrix, return the EMD between the two vectors

    :param P: 1D numpy array (float64)
    :param Q: 1D numpy array (float64)
    :param distance_matrix: 2D numpy array giving the distance between each pair of nodes (float64)
    :return: scalar EMD value
    """
    # check if P and Q are np arrays
    if type(P) is not np.ndarray:
        P = np.array(P)
    if type(Q) is not np.ndarray:
        Q = np.array(Q)
    if type(distance_matrix) is not np.ndarray:
        distance_matrix = np.array(distance_matrix)
    # check if P and Q are normalized to sum to 1, or close enough
    if not np.isclose(P.sum(), 1.0):
        raise ValueError('P is not normalized to sum to 1')
    if not np.isclose(Q.sum(), 1.0):
        raise ValueError('Q is not normalized to sum to 1')
    # cast everything to float64_t
    P = P.astype(np.float64)
    Q = Q.astype(np.float64)
    distance_matrix = distance_matrix.astype(np.float64)
    if with_flow:
        return emd_with_flow(P, Q, distance_matrix)
    return emd(P, Q, distance_matrix)


def get_root_of_tree(G:nx.DiGraph):
    """
    Returns the root node of a directed tree

    :param G: directed tree
    :return: root node name
    """
    roots = [n for n, d in G.in_degree() if d == 0]
    if len(roots) > 1:
        raise Exception(f"The graph has multiple roots: {roots}")
    return roots[0]


# get all leaf descendants of a certain node
def get_leaf_descendants_from_node(G, node):
    """
    Return all leaf descendants of a node, excluding the node itself.

    :param G: graph
    :param node: node
    :return: set of leaf nodes descending from the node
    """
    descendants = set()
    for n in nx.descendants(G, node):
        if G.out_degree(n) == 0:
            descendants.add(n)
    return descendants


def get_leaf_nodes(G):
    """
    Return all leaf nodes of a graph

    :param G: graph
    :return: set of leaf nodes
    """
    leaf_nodes = set()
    for n in G.nodes():
        if G.out_degree(n) == 0:
            leaf_nodes.add(n)
    leaf_nodes = sorted(list(leaf_nodes))
    return leaf_nodes


def get_descendants(G, node):
    """
    Return all descendants of a node, including the node itself.

    :param G: networkx graph
    :param node: name of a node
    :return: set of nodes
    """
    descendants = set()
    descendants.add(node)
    for n in nx.descendants(G, node):
        descendants.add(n)
    return descendants


def get_descendant(graph, v1, v2):
    """
    of the two nodes v1 and v2, ASSUMED TO BE ADJACENT, find out which one is the descendant of the other

    :param graph: networkx graph, directed
    :param v1: node name
    :param v2: node name
    :return: descendant node name
    """
    if v1 in graph.predecessors(v2):
        return v2
    elif v2 in graph.predecessors(v1):
        return v1
    else:
        raise ValueError(f"Nodes {v1} and {v2} are not adjacent")


def get_KO_labels_and_index(distances_labels_file, basis=None):
    """
    Given a file containing the basis for the rows of the pairwise KO distance matrix, return a list of labels and a
    dictionary mapping labels to indices. Only include the labels in the basis if specified.

    :param distances_labels_file: text file containing the labels from the output of sourmash compare
    :param basis: list of labels to include in the output
    :return: (list of labels, dictionary mapping labels to indices)
    """
    # import label names
    pairwise_dist_KOs = []
    with open(distances_labels_file, 'r') as f:
        for line in f.readlines():
            ko = line.strip().split('|')[-1]  # KO number is the last in the pip-delim list
            ko = ko.split(':')[-1]  # remove the ko: prefix
            pairwise_dist_KOs.append(ko)
    pairwise_dist_KO_index = {node: i for i, node in enumerate(pairwise_dist_KOs)}
    if basis is not None:
        # remove those that are not in the basis
        pairwise_dist_KOs = [x for x in pairwise_dist_KOs if x in basis]
        # also remove them from the pairwise dist KO index
        pairwise_dist_KO_index = {x: pairwise_dist_KO_index[x] for x in pairwise_dist_KOs}
    return pairwise_dist_KOs, pairwise_dist_KO_index


def get_graphs_and_index(edge_list_file, brite):
    """
    Given an edge list file, return a networkx graph and a dictionary mapping node names to indices

    :param edge_list_file: text file containing the edge list
    :param brite: which root to pick in the BRITE hierarchy
    :return: (directed networkx graph, undirected networkx graph, basis of node names, dictionary mapping node names to
    indices in that basis)
    """
    if brite not in BRITES:
        raise ValueError(f"brite must be one of {BRITES}. I received {brite}.")
    G = import_graph(edge_list_file, directed=True)
    # get the descendants of the brite
    descendants = get_descendants(G, brite)
    # add the root node back in (this will only have affect if you are using multiple brites. TODO: not yet implemented)
    descendants.add('root')  #TODO: see if I actually want to add the root back in or not
    # select the subgraph from the brite to the leaves
    G = G.subgraph(descendants)
    G_undirected = G.to_undirected()
    # set the basis for the tree, which is an ordering of the edges. I'll identify each edge by its terminal node
    basis = [x for x in G.nodes()]
    basis_index = {node: i for i, node in enumerate(basis)}
    return G, G_undirected, basis, basis_index


# Create a class that will generate simulated distributions
class LeafDistributionSimulator:
    def __init__(self, edge_list_file, brite):
        Gdir, Gundir, basis, basis_index = get_graphs_and_index(edge_list_file, brite)
        self.Gdir = Gdir
        self.Gundir = Gundir
        self.basis = basis
        self.basis_index = basis_index
        self.brite = brite
        self.leaf_nodes = get_leaf_descendants_from_node(Gdir, brite)
        self.leaf_nodes_index = {node: basis_index[node] for node in self.leaf_nodes}

    def get_random_dist_on_leaves(self):
        """
        Return a random distribution on the leaves of the tree

        :return: 1D numpy array, indexed by nodes in the graph, but only supported on the leaves
        """
        # Use a dirichlet distribution to generate a random distribution on the leaves, since this will sum to 1
        dist = np.random.dirichlet(np.ones(len(self.leaf_nodes)))
        # map the distribution to the basis
        P = np.zeros(len(self.basis))
        for i, node in enumerate(self.leaf_nodes):
            P[self.leaf_nodes_index[node]] = dist[i]
        return P


########################################################################################################################
# TODO: the following are works in progress
def get_EMD_from_edge_file(edge_file, branch_len_function):  #FIXME: branch_len_function is not used
    df = parse_edge_list(edge_file)
    df['length'] = [1.] * len(df)
    leaf_nodes = get_leaf_nodes(edge_file)
    distance_matrix, node_list = get_distance_matrix_from_edge_list(df)
    index_dict = get_ID_index_dict(node_list)
    P = simulate_leaf_supported_vector(leaf_nodes, len(node_list), index_dict)
    Q = simulate_leaf_supported_vector(leaf_nodes, len(node_list), index_dict)
    # print(norm_P)
    start_time = time.time()
    emd_value = get_EMD_pyemd(P, Q, distance_matrix)
    print(emd_value)
    print("Total process time:", time.time() - start_time)
    return


def get_leaf_nodes2(edge_list_file):
    df = parse_edge_list(edge_list_file)
    child_nodes = set(df['child'])
    parent_nodes = set(df['parent'])
    internal_nodes = parent_nodes.intersection(child_nodes)
    leaf_nodes = child_nodes.difference(internal_nodes)
    return leaf_nodes


def get_ID_index_dict(node_list):
    if type(node_list) is not list:
        node_list = list(node_list)
    index_dict = dict(zip(node_list, range(len(node_list))))
    return index_dict


def simulate_leaf_supported_vector(leaf_nodes, length, index_dict):
    vector = [0.]*length
    leaf_index = [index_dict[i] for i in leaf_nodes]
    for i in leaf_index:
        vector[i] = np.random.rand()
    norm_v = vector / np.linalg.norm(vector)
    return norm_v


#temp
def get_leaf_nodes_only_graph():
    leaf_nodes = get_leaf_nodes('data/kegg_ko_edge_df.txt')
    df = pd.DataFrame(columns=['parent', 'child'])
    root = 'root'
    df['child'] = list(leaf_nodes)
    df['parent'] = root
    # FIXME: this is hard coded
    df.to_csv('data/kegg_ko_leaf_only_df.txt', sep='\t', index=None)
    #print(df)

def get_profile_from_sourmash(sourmash_file, save_as, normalize = True):
    '''
    Given a sourmash file, extract the profile vector and write it to a file with the name given by param save_as
    :param sourmash_file: path to sourmash file
    :param name for the vector file to be saveved. Should be .txt format
    :return:
    '''
    df = pd.read_csv(sourmash_file)
    id = df['name']
    #id = list(map(lambda x: x.split(':')[1], id)) #if want to remove the "ko" in front
    abund = list(df['unique_intersect_bp'])
    df = pd.DataFrame(df, columns=['KO','rel_abund'])
    if normalize is True and not np.isclose(np.sum(abund), 1.0):
        print('sum of abundance before normalization ', np.sum(abund))
        abund = abund / np.sum(abund)
        print('sum of abundance after normalization ', np.sum(abund))
    df['KO'] = id
    df['rel_abund'] = abund
    df.to_csv(save_as, index=None, sep='\t')


def get_EMDUniFrac_from_functional_profiles(profile1, profile2, distance_matrix, node_list):
    start = time.time()
    print(profile1)
    df1 = pd.read_csv(profile1)
    id1 = df1['name']
    id1 = list(map(lambda x: x.split(':')[1], id1))
    df2 = pd.read_csv(profile2)
    id2 = df2['name']
    id2 = list(map(lambda x: x.split(':')[1], id2))
    abund1 = list(df1['unique_intersect_bp'])
    abund2 = list(df2['unique_intersect_bp'])
    print(abund1)
    sample_vector1 = [0.]*len(node_list)
    sample_vector2 = [0.]*len(node_list)

    for i,id in enumerate(node_list):
        if id in id1:
            sample_vector1[i] = abund1[id1.index(id)]
        elif id in id2:
            sample_vector2[2] = abund2[id2.index(id)]
    normed_sample1 = sample_vector1 / np.linalg.norm(sample_vector1)
    normed_sample2 = sample_vector2 / np.linalg.norm(sample_vector2)
    unifrac = emd(normed_sample1, normed_sample2, distance_matrix)
    print(unifrac)
    print(time.time() - start)
    return


# FIXME: these test don't use asserts and have a lot of stuff hard coded

def test_parse_edge_list():
    df = parse_edge_list('data/kegg_ko_edge_df.txt')
    print(df)


def test_get_matrix_from_edge_list():
    df = parse_edge_list('data/kegg_ko_edge_df.txt')
    df['length'] = [1.] * len(df)
    print(df)
    distance_matrix, node_list = get_distance_matrix_from_edge_list(df)
    print(distance_matrix)
    print(list(node_list)[:10])


def test_get_leaf_nodes():
    leaf_nodes = get_leaf_nodes('data/kegg_ko_edge_df.txt')
    print(leaf_nodes)


def test_get_EMD():
    df = parse_edge_list('data/kegg_ko_leaf_only_df.txt')
    df['length'] = [1.] * len(df)
    leaf_nodes = get_leaf_nodes('data/kegg_ko_leaf_only_df.txt')
    distance_matrix, node_list = get_distance_matrix_from_edge_list(df)
    index_dict = get_ID_index_dict(node_list)
    P = simulate_leaf_supported_vector(leaf_nodes, len(node_list), index_dict)
    Q = simulate_leaf_supported_vector(leaf_nodes, len(node_list), index_dict)
    #print(norm_P)
    start_time = time.time()
    emd_value = get_EMD_pyemd(P, Q, distance_matrix)
    print(emd_value)
    print("Total process time:", time.time() - start_time)


def test_simulate_leaf_supported_vector():
    leaf_nodes = get_leaf_nodes('data/kegg_ko_edge_df.txt')
    df = parse_edge_list('data/kegg_ko_edge_df.txt')
    df['length'] = [1.] * len(df)
    distance_matrix, node_list = get_distance_matrix_from_edge_list(df)
    index_dict = get_ID_index_dict(node_list)
    sparse_vector = simulate_leaf_supported_vector(leaf_nodes, len(distance_matrix), index_dict)
    print(sparse_vector)


def test_get_EMDUniFrac_from_profiles():
    profile1 = 'data/SRS1041031.denovo_duplicates_marked.trimmed_KOs_sketched_scaled_10.sig.zip_gather_k_5.csv'
    profile2 = 'data/SRS893174.denovo_duplicates_marked.trimmed_KOs_sketched_scaled_10.sig.zip_gather_k_5.csv'
    leaf_nodes = get_leaf_nodes('data/kegg_ko_edge_df.txt')
    df = parse_edge_list('data/kegg_ko_edge_df.txt')
    df['length'] = [1.] * len(df)
    distance_matrix, node_list = get_distance_matrix_from_edge_list(df)
    get_EMDUniFrac_from_functional_profiles(profile1, profile2, distance_matrix, node_list)


def make_edge_list_file_len_1_tmp():
    df = parse_edge_list('data/kegg_ko_edge_df.txt')
    df['length'] = [1.] * len(df)
    df.to_csv('data/kegg_ko_edge_df_len1.txt', sep='\t', index=None)


