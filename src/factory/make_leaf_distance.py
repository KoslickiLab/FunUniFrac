import numpy as np
import src.objects.func_leaf_distance as func_leaf_distance


def get_KO_pairwise_dist(distance_file, distances_label_file) -> func_leaf_distance.LeafPairwiseDistances:
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
    
    return func_leaf_distance.LeafPairwiseDistances(pw_dist, KO_dist_labels, KO_dist_indices)