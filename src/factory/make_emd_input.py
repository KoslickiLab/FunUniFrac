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

