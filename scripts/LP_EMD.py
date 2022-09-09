from phylodm import PhyloDM
from pyemd import emd
import numpy as np
import pandas as pd
import networkx
import time
import argparse

def get_dm_from_tree_file(tree_file):
    pdm = PhyloDM.load_from_newick_path(tree_file)
    dm = pdm.dm(norm=False)
    return dm

def parse_edge_list(file):
    df = pd.read_table(file)
    return df

def get_matrix_from_edge_list(data):
    '''
    Assume has the form head, tail, length
    :param df:
    :return:
    '''
    edgeList = data.values.tolist()
    G = networkx.DiGraph()
    for i in range(len(edgeList)):
        G.add_edge(edgeList[i][0], edgeList[i][1], weight=edgeList[i][2])
    A = networkx.adjacency_matrix(G).A
    return A, G.nodes()

def get_EMD(P, Q, distance_matrix):
    return emd(P, Q, distance_matrix)

def get_EMD_from_edge_file(edge_file):
    df = parse_edge_list(edge_file)
    df['length'] = [1.] * len(df)
    leaf_nodes = get_leaf_nodes(edge_file)
    distance_matrix, node_list = get_matrix_from_edge_list(df)
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
    distance_matrix, node_list = get_matrix_from_edge_list(df)
    print(distance_matrix)
    print(list(node_list)[:10])

def test_get_leaf_nodes():
    leaf_nodes = get_leaf_nodes('data/kegg_ko_edge_df.txt')
    print(leaf_nodes)

def test_get_EMD():
    df = parse_edge_list('data/kegg_ko_leaf_only_df.txt')
    df['length'] = [1.] * len(df)
    leaf_nodes = get_leaf_nodes('data/kegg_ko_leaf_only_df.txt')
    distance_matrix, node_list = get_matrix_from_edge_list(df)
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
    distance_matrix, node_list = get_matrix_from_edge_list(df)
    index_dict = get_ID_index_dict(node_list)
    sparse_vector = simulate_leaf_supported_vector(leaf_nodes, len(distance_matrix), index_dict)
    print(sparse_vector)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get testing statistics of classification test.')
    parser.add_argument('-e', '--edge_file', type=str, help='An edge list file.')

    args = parser.parse_args()
    edge_file = args.edge_file
    get_EMD_from_edge_file(edge_file)
    #test_get_dm_from_tree_file()
    #test_parse_edge_list()
    #test_get_matrix_from_edge_list()
    #test_get_EMD()
    #test_get_leaf_nodes()
    #test_simulate_leaf_supported_vector()
    #get_leaf_nodes_only_graph()
