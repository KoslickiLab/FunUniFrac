import src.objects.func_tree as func_tree
import networkx as nx


def weighted_tree_to_EMDU_input(tree: func_tree.FuncTree, brite=None, edge_len_property=None):
    """
    Convert a weighted graph to the format required by EMDUniFrac: Tint, lint, and nodes_in_order.
    Since EMDUniFrac wants everything to be integers, also return a mapping from integers to node names.

    :param G: weighted networkx graph
    :return: (Tint dict, lint dict, nodes_in_order array, index_2_node dict)
    """
    if brite is not None:
        G = tree.make_subtree(brite)
    else:
        G = tree.main_tree
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
