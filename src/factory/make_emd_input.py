from src.objects.func_tree import FuncTree, FuncTreeEmduInput
from src.objects.profile_vector import FuncProfileVector
import networkx as nx
import numpy as np
import pandas as pd
import multiprocessing
from itertools import repeat


def tree_to_EMDU_input(tree: FuncTree, brite=None, edge_len_property=None) -> FuncTreeEmduInput:
    """
    Convert a weighted graph to the format required by EMDUniFrac: Tint, lint, and nodes_in_order.
    Since EMDUniFrac wants everything to be integers, also return a mapping from integers to node names.

    :param G: weighted networkx graph
    :return: (Tint dict, lint dict, nodes_in_order array, index_2_node dict)
    """
    if brite is not None:
        G = tree.set_subtree(brite)
    else:
        G = tree.main_tree
    # first determine the edge_length attribute
    if edge_len_property is None:
        edge_len_property = FuncTree.infer_edge_len_property(G)
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
            # Tint[node] = "-1"  # TODO: check if this is how we do it in EMDUniFrac
    lint = {}
    for edge_data in G.edges(data=True):
        i, j, data = edge_data
        weight = data[edge_len_property]
        lint[(i, j)] = weight
        lint[(j, i)] = weight
    root = FuncTree.get_root_of_tree(G)
    nodes_in_order = list(nx.dfs_postorder_nodes(G, source=root))
    # Tint, lint, and nodes_in_order require these to be integers, so let's map everything to that
    node_2_EMDU_index = {node: i for i, node in enumerate(nodes_in_order)}
    # convert the Tint dict
    Tint = {node_2_EMDU_index[node]            : node_2_EMDU_index[Tint[node]] for node in Tint}
    # convert the lint dict
    lint = {(node_2_EMDU_index[i], node_2_EMDU_index[j])            : lint[(i, j)] for i, j in lint}
    # convert the nodes in order
    nodes_in_order = [node_2_EMDU_index[node] for node in nodes_in_order]
    # convert the integers to nodes
    EMDU_index_2_node = {i: node for node, i in node_2_EMDU_index.items()}
    return FuncTreeEmduInput(Tint, lint, nodes_in_order, EMDU_index_2_node)



def get_func_profiles_parallel(args):
    def map_func(file, input: FuncTreeEmduInput, abundance_key, unweighted):
        df = pd.read_csv(file)
        P: FuncProfileVector = functional_profile_to_vector(df, input, abundance_key=abundance_key, normalize=True)
        if unweighted:
            # if entries are positive, set to 1
            P[P > 0] = 1
        return file, P
    return map_func(*args)


def parallel_functional_profile_to_vector(num_threads, fun_files, input, abundance_key, unweighted)-> list[tuple[pd.DataFrame, FuncProfileVector]]:
    pool = multiprocessing.Pool(num_threads)
    results = pool.imap(get_func_profiles_parallel, 
                        zip(fun_files, 
                            repeat(input), 
                            repeat(abundance_key), 
                            repeat(unweighted)), 
                            chunksize=max(2, len(fun_files) // num_threads))
    pool.close()
    pool.join()
    return results


def functional_profile_to_vector(functional_profile: pd.DataFrame, input: FuncTreeEmduInput, abundance_key='median_abund', normalize=True):
    """This function will take a sourmash functional profile and convert it to a form that can be used by EMDUniFrac
    :param functional_profile_file: csv file output from `sourmash gather`
    :param EMDU_index_2_node: dictionary that translates between the indices used by EMDUniFrac and the actual graph
    node names
    :param abundance_key: key in the functional profile that contains the abundance information
    :param normalize: whether to normalize the abundance information
    :return: vector P
    """
    # reverse the dictionary for convenience
    EMDU_index_2_node = input.idx_to_node
    node_2_EMDU_index = {v: k for k, v in EMDU_index_2_node.items()}
    # import the functional profile
    df = functional_profile
    # get the functional profile as a vector
    P = np.zeros(len(EMDU_index_2_node))
    for i, row in df.iterrows():
        try:
            abundance = float(row[abundance_key])
        except ValueError:
            raise ValueError(
                f"The abundance key {abundance_key} gives a value that is not a number: {row[abundance_key]}")
        try:
            if row['name'].startswith('ko:'):
                ko = row['name'].split(':')[-1]  # ko:K12345 case
            else:
                # abc|xyz|lmnop|ko:K12345 case
                ko = row['name'].split('|')[-1].split(':')[-1]
            if ko not in node_2_EMDU_index:
                print(f"Warning: {ko} not found in EMDU index, skipping.")
            else:
                P[node_2_EMDU_index[ko]] = abundance
        except:
            raise Exception(
                f"Could not parse the name {row['name']} in the functional profile")
    if P.sum() == 0:
        raise Exception(f"Functional profile is empty! Perhaps try a different abundance key? "
                        f"You used {abundance_key}.")
    if normalize:
        P = P / P.sum()
    print(P.sum())
    return P
