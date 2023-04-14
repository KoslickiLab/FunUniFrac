import os
import networkx as nx
from src.objects.func_tree import FuncTree, FuncTreeEmduInput


def import_graph(edge_list_file, directed=True) -> FuncTree:
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
    tree = FuncTree(G)
    return tree
