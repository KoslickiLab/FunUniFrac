from phylodm import PhyloDM
from pyemd import emd
import numpy as np
import pandas as pd
import networkx

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
    return A

def get_EMD(P, Q, distance_matrix):
    return emd(P, Q, distance_matrix)

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
    distance_matrix = get_matrix_from_edge_list(df)
    print(len(distance_matrix))


def test_get_EMD():
    df = parse_edge_list('data/kegg_ko_edge_df.txt')
    df['length'] = [1.] * len(df)
    distance_matrix = get_matrix_from_edge_list(df)
    P = np.random.rand(len(distance_matrix))
    norm_P = P/np.linalg.norm(P)
    Q = np.random.rand(len(distance_matrix))
    norm_Q = Q / np.linalg.norm(Q)
    #print(norm_P)
    emd_value = get_EMD(norm_P, norm_Q, distance_matrix)
    print(emd_value)


if __name__ == '__main__':
    #test_get_dm_from_tree_file()
    #test_parse_edge_list()
    #test_get_matrix_from_edge_list()
    test_get_EMD()