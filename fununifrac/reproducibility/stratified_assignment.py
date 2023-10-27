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
    parser.add_argument('-a', '--alpha', help='factor for stratified assignment.', type=float, default=0.1)

    args = parser.parse_args()
    #edge_list = get_data_abspath(args.edge_list)
    edge_list = args.edge_list
    kegg_tree = get_KeggTree_from_edgelist(edge_list, edge_length=False)
    kegg_tree.group_nodes_by_depth()
    kegg_tree.make_full_tree()
    kegg_tree.get_needed_pairs()

    edge_lengths_solution = stratified_assignment(kegg_tree, args.alpha)
    write_edge_list_preserve_order(edge_lengths_solution, edge_list, args.save)



if __name__ == "__main__":
    main()