import networkx as nx
import src.utils.kegg_db as kegg_db
import src.base.pairwise_dist as pairwise_dist
import numpy as np


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


def get_KO_pairwise_dist(distance_file, distances_label_file) -> pairwise_dist.PairwiseDistance:
    """ 
    Given KO distance files, return pairwise distance object.
    :param distances_file: npz file containing the distances from the output of sourmash compare
    :param distances_labels_file: text file containing the labels from the output of sourmash compare
    """
    pw_dist = np.load(distance_file)
    KO_dist_labels = []
    with open(distances_label_file, 'r') as f:
        for line in f.readlines():
            ko = line.strip().split('|')[-1]  # KO number is the last in the pip-delim list
            ko = ko.split(':')[-1]  # remove the ko: prefix
            KO_dist_labels.append(ko)
    KO_dist_indices = {node: i for i, node in enumerate(KO_dist_labels)}
    
    return pairwise_dist.PairwiseDistance(pw_dist, KO_dist_labels, KO_dist_indices)
