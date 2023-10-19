import argparse
import logging
import pandas as pd
from src.objects.func_tree import FuncTreeEmduInput
import src.factory.make_tree as make_tree
import src.factory.make_emd_input as make_emd_input
from src.algorithms.emd_unifrac import EarthMoverDistanceUniFracSolver
import numpy as np


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
    parser.add_argument('-e', '--edge_list',
                        help='Input edge list file of the KEGG hierarchy. Must have lengths in the '
                             'third column.')
    parser.add_argument('-f', '--file', help="Input file that contains a combined dataframe of all samples with column"
                                             "names being samples and row names being KOs.")
    parser.add_argument('-o', '--out_dir', help='Output directory name.')
    parser.add_argument('-i', '--out_id',
                        help='Test purpose: give an identifier to the output file so that tester can recognize it')
    parser.add_argument('-a', '--abundance_key',
                        help='Key in the gather results to use for abundance. Default is `f_unique_weighted`',
                        default='f_unique_weighted')
    parser.add_argument('-b', '--brite', help='Use the subtree of the KEGG hierarchy rooted at the given BRITE ID. '
                                              'eg. brite:ko00001', default="ko00001", type=str, required=False)
    parser.add_argument('--diffab', action='store_true', help='Also return the difference abundance vectors.')
    parser.add_argument('-v', help="Be verbose", action="store_const", dest="loglevel", const=logging.INFO,
                        default=logging.WARNING)
    parser.add_argument('--unweighted', help="Compute unweighted unifrac instead of the default weighted version",
                        action="store_true")
    parser.add_argument('--L2', help="Use L2 UniFrac instead of L1", action="store_true")
    parser.add_argument('--Ppushed', help="Flag indicating you want the pushed vectors to be saved.",
                        action="store_true")

    args = parser.parse_args()
    solver = EarthMoverDistanceUniFracSolver()
    # Parse the graph and get the FunUniFrac objects Tint, length, and edge_list
    tree = make_tree.import_graph(args.edge_list, directed=True)
    input: FuncTreeEmduInput = make_emd_input.tree_to_EMDU_input(tree, args.brite)
    sample_df = pd.read_csv(args.file, index_col='name')
    Ps_pushed = {}
    for col in sample_df.columns:
        P = make_fununifrac_inputs(sample_df[col], input)
        P_pushed = solver.push_up_L1(P, input)
        Ps_pushed[col] = P_pushed
    dists, diffabs_sparse = solver.pairwise_computation(Ps_pushed, sample_df.columns, input, args.diffab, args.L2)
    out_id = args.out_id if args.out_id else "fununifrac_out"
    out_file = f"{args.out_dir}/{out_id}.npy"
    np.save(out_file, dists)
    basis_out = f"{args.out_dir}/{out_id}.basis.npy"
    with open(basis_out, 'w') as f:
        for file in sample_df.columns:
            f.write(f"{file}\n")

if __name__ == '__main__':
    main()

