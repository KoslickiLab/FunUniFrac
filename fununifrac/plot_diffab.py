import os, sys, argparse
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(ROOT_DIR)
import multiprocessing
import logging
from src.utility.differential_abundance import plot_diffab
from src.objects.func_tree import FuncTreeEmduInput
import src.factory.make_tree as make_tree
import src.factory.make_emd_input as make_emd_input
from src.algorithms.emd_unifrac import EarthMoverDistanceUniFracSolver
import numpy as np
import pandas as pd
import itertools as it


def make_fununifrac_inputs(raw_P, input, normalize=True):
    #convert an array into one that's suitable for use
    EMDU_index_2_node = input.idx_to_node
    node_2_EMDU_index = {v: k for k, v in EMDU_index_2_node.items()}

    if normalize:
        raw_P = raw_P/raw_P.sum()
    P = np.zeros(len(EMDU_index_2_node))
    for ko in raw_P.index:
        if ko not in node_2_EMDU_index:
            print(f"Warning: {ko} not found in EMDU index, skipping.")
        else:
            P_index = node_2_EMDU_index[ko]
            P[P_index] = raw_P[ko]
    return P

def main():
    parser = argparse.ArgumentParser(description="This script will take a directory of sourmash gather results and "
                                                 "a weighted edge list representing the KEGG hierarchy and compute "
                                                 "all pairwise functional unifrac distances.")
    parser.add_argument('-df', '--dataframe', help='Dataframe of combined samples')
    parser.add_argument('-e', '--edge_list',
                        help='Input edge list file of the KEGG hierarchy. Must have lengths in the '
                             'third column.', required=True)
    parser.add_argument('-fd', '--file_dir', help="Directory of sourmash files.", required=False)
    parser.add_argument('-fp', '--file_pattern', help="Pattern to match files in the directory. Default is "
                                                      "'*_gather.csv'", default='*_gather.csv')
    parser.add_argument('-o', '--out_dir', help='Output directory name.', required=True)
    parser.add_argument('-i', '--out_id', help='Test purpose: give an identifier to the output file so that tester can recognize it')
    parser.add_argument('-f', '--force', help='Overwrite the output file if it exists', action='store_true')
    parser.add_argument('-a', '--abundance_key',
                        help='Key in the gather results to use for abundance. Default is `f_unique_weighted`',
                        default='f_unique_weighted')
    parser.add_argument('-t', '--threads', help='Number of threads to use. Default is half the cores available.',
                        default=int(multiprocessing.cpu_count() / 2), type=int)
    parser.add_argument('-b', '--brite', help='Use the subtree of the KEGG hierarchy rooted at the given BRITE ID. '
                                              'eg. brite:ko00001', default="ko00001", type=str, required=False)
    parser.add_argument('--diffab', action='store_true', help='Also return the difference abundance vectors.')
    parser.add_argument('-v', help="Be verbose", action="store_const", dest="loglevel", const=logging.INFO,
                        default=logging.WARNING)
    parser.add_argument('--unweighted', help="Compute unweighted unifrac instead of the default weighted version",
                        action="store_true")
    parser.add_argument('--L2', help="Use L2 UniFrac instead of L1", action="store_true")
    parser.add_argument('--Ppushed', help="Flag indicating you want the pushed vectors to be saved.", action="store_true")


    # call main function
    args = parser.parse_args()
    solver = EarthMoverDistanceUniFracSolver()
    # Parse the graph and get the FunUniFrac objects Tint, length, and edge_list
    tree = make_tree.import_graph(args.edge_list, directed=True)
    input: FuncTreeEmduInput = make_emd_input.tree_to_EMDU_input(tree, args.brite)
    if args.dataframe.endswith('.csv'):
        sample_df = pd.read_csv(args.dataframe, index_col=0)
    else:
        sample_df = pd.read_table(args.dataframe, index_col=0)
    sample_dict = dict()
    for c in sample_df.columns:
        sample_dict[c] = sample_df[c]
    print(sample_dict)
    for p, q in it.combinations(sample_dict.keys(), 2):
        print(f"{p} vs {q}")
    (Z, F, diffab) = solver.solve(input, sample_dict[p], sample_dict[q], weighted=True)
    print(f"Z: {Z}")
    print(diffab)
    #plot_diffab(nodes_in_order, diffab, P_label, Q_label, plot_zeros=True)

if __name__ == '__main__':
    main()


