#from phylodm import PhyloDM
from pyemd import emd
import numpy as np
import pandas as pd
import networkx
import time
import networkx as nx
import argparse

#not really necessary. May or may not be used. If not needed, remove phylodm as well.
def get_dm_from_tree_file(tree_file):
    pdm = PhyloDM.load_from_newick_path(tree_file)
    dm = pdm.dm(norm=False)
    return dm

###delete
def parse_edge_list(file):
    df = pd.read_table(file)
    return df
###

def get_distance_matrix_from_edge_list(edge_list_file):
    '''
    Given an edge list file, with three columns: head, tail, length, return a distance matrix and a list of nodes
    :param edge_list_file: file name of the edge list
    :return: distance matrix and a list of nodes indexing the distance matrix
    '''
    with open(edge_list_file, 'rb') as fid:
        # read the first line
        line = fid.readline()
        try:
            line = line.decode('utf-8')
            col1, col2, col3 = line.strip().split('\t')
        except:
            raise ValueError('The first line of the file must have three columns: parent child edge_length')
        # import the graph. First two columns are the nodes, last column is the weight
        G = nx.read_edgelist(fid, delimiter='\t', data=((col3, float),))
    D = nx.all_pairs_dijkstra_path_length(G, weight=col3)
    return D, G.nodes()

def get_EMD(P, Q, distance_matrix):
    return emd(P, Q, distance_matrix)

def get_EMD_from_edge_file(edge_file, branch_len_function):
    df = parse_edge_list(edge_file)
    df['length'] = [1.] * len(df)
    leaf_nodes = get_leaf_nodes(edge_file)
    distance_matrix, node_list = get_distance_matrix_from_edge_list(df)
    index_dict = get_ID_index_dict(node_list)
    P = simulate_leaf_supported_vector(leaf_nodes, len(node_list), index_dict)
    Q = simulate_leaf_supported_vector(leaf_nodes, len(node_list), index_dict)
    # print(norm_P)
    start_time = time.time()
    emd_value = get_EMD(P, Q, distance_matrix)
    print(emd_value)
    print("Total process time:", time.time() - start_time)
    return

def get_leaf_nodes(edge_list_file):
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
    df.to_csv('data/kegg_ko_leaf_only_df.txt', sep='\t', index=None)
    #print(df)

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


#tests
def test_get_dm_from_tree_file():
    dm = get_dm_from_tree_file('data/test_newick.tree')
    print(dm)

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
    emd_value = get_EMD(P, Q, distance_matrix)
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

