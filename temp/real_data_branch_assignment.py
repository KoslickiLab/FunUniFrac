import networkx as nx
import pandas as pd
from itertools import combinations
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../data')
sys.path.append('..')
import argparse
from Kegg_tree import *


def main():
    parser = argparse.ArgumentParser(description="At this point the script takes in an edge-list file WITH branch lengths"
                                                 " and tries to recompute the edge lengths using the algorithm and compare "
                                                 "the accuracy. With that being said, any sub tree would need to be created"
                                                 " prior to passing into the script. ")
    parser.add_argument('-e', '--edge_list',
                        help='Input edge list file of the KEGG hierarchy. Must have lengths in the '
                             'third column.', required=True)
    parser.add_argument('-s', '--save', help='Path to save the output file.', default='edge_length_solution.png')
    parser.add_argument('-dm', '--dist_matrix', help='Pairwise KO distance matrix file.', default='')

    args = parser.parse_args()
    #edge_list = get_data_abspath(args.edge_list)
    edge_list = args.edge_list
    kegg_tree = get_KeggTree_from_edgelist(edge_list)
    edge_lengths_solution = {}
    pw_dist = kegg_tree.pw_dist
    print("preparation done.")
    assign_branch_lengths(kegg_tree, kegg_tree.leaf_nodes, pw_dist, edge_lengths_solution)

    visualize_diff(edge_lengths_solution, kegg_tree, args.save)


if __name__ == "__main__":
    main()