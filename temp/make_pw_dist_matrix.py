from Kegg_tree import get_KeggTree_from_edgelist
import argparse
import numpy as np
import time

"""
This script takes an edge list file with branch lengths and converts the distances between leaves into a pairwise
distance matrix and returns a pairwise distance matrix file and a label file.
"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="edge list file with branch lengths")
    parser.add_argument("-p", "--prefix", help="prefix for writing out files")
    args = parser.parse_args()

    kegg_tree = get_KeggTree_from_edgelist(args.file)
    kegg_tree.get_pw_dist() #in dict form
    pw_dist = np.ndarray((len(kegg_tree.leaf_nodes), len(kegg_tree.leaf_nodes)))
    for i, node in enumerate(kegg_tree.leaf_nodes):
        for j, another_node in enumerate(kegg_tree.leaf_nodes):
            if i == j:
                pw_dist[i][j] = 0
            else:
                pw_dist[i][j] = pw_dist[j][i] = kegg_tree.pw_dist[(node, another_node)] if node < \
                another_node else kegg_tree.pw_dist[(another_node, node)]
    print(pw_dist)
    np.save(args.prefix, pw_dist)
    label_file = open(args.prefix+"_labels.txt", 'w+')
    label_file.write('\n'.join(kegg_tree.leaf_nodes))


if __name__ == '__main__':
    main()
