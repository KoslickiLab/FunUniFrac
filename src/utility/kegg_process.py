import networkx as nx
import src.utility.kegg_db as kegg_db


def make_nodes_readable(G):
    """
    This function takes a graph with nodes that are KEGG IDs and makes them readable.

    :param G: networkx graph
    :return: networkx graph with renamed nodes
    """
    # get the node names
    node_names = list(G.nodes())
    rename_map = dict()
    nodes_added = set()
    dup_iter = 0
    for old_node in node_names:
        # default to the old name
        rename_map[old_node] = old_node
        # try to see if it's a KO term
        new_node = kegg_db.instance.koids_human_names.get(old_node)
        if not new_node:
            # try if it's a BRITE term
            new_node = kegg_db.instance.brite_names.get(old_node)
        if not new_node:
            # see if it's a pathway and starts with a number
            first_part = old_node.split()[0]
            if first_part.isdigit():
                new_node = ' '.join(old_node.split()[1:])
        # if successful, use the new node name
        if new_node:
            # first check if we've seen this name before (IDs to names in KEGG are many to one)
            if new_node not in nodes_added:
                rename_map[old_node] = new_node.strip()
                nodes_added.add(new_node)
            else:
                new_node = new_node + f"dup_{dup_iter}"
                rename_map[old_node] = new_node.strip()
                nodes_added.add(new_node)
                dup_iter += 1
    Gret = nx.relabel_nodes(G, rename_map)
    return Gret



