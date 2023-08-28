#!/usr/bin/env python
import os, sys
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print(ROOT_DIR)
sys.path.append(ROOT_DIR)
sys.path.append(ROOT_DIR + '../src')
sys.path.append('../../')
from src.commands import compute_all_pairwise_fununifrac
import src.utility.kegg_db as kegg_db
from src.algorithms.emd_unifrac import EarthMoverDistanceUniFracSolver
import argparse
import pandas as pd
from sklearn.cluster import KMeans
import src.factory.make_tree as make_tree
import src.factory.make_emd_input as make_emd_input


TREE_FILE = ROOT_DIR + '/../data/kegg_ko_edge_df_br_ko00001.txt_AAI_lengths_n_50_f_10_r_100.txt'

def parsearg():
    parser = argparse.ArgumentParser(description="This script clusters the samples according to fununifrac distance"
                                                 "for a given .csv file")
    parser.add_argument("-f", "--file", type=str, help="Path to the .csv file with row indexed by KO and "
                                                       "columns indexed by sample names.")
    parser.add_argument("--fast", action=argparse.BooleanOptionalAction, help="If specified, clustering will be done"
                                                                              "by pushing up all vectors and then "
                                                                              "cluster using Kmeans")
    parser.add_argument('-b', '--brite', help='Use the subtree of the KEGG hierarchy rooted at the given BRITE ID. '
                                              'eg. brite:ko00001', default="ko00001", type=str, required=False)
    parser.add_argument('-m', '--metadata_file', help='Path to the metadata file.')
    return parser.parse_args()

def extract_vectors_from_csv(file):
    df = pd.read_csv(file, index_col=0)
    df.fillna(0, inplace=True)
    KO_index = df.index
    sample_list = df.columns
    vectors = [df[sample_id].to_list() for sample_id in sample_list]
    print(len(vectors[0]))
    return KO_index, vectors

def cluster():
    '''
    push up all vectors
    :return:
    '''
    pass



if  __name__ == "__main__":
    args = parsearg()
    if args.fast:
        print("fast: using KMeans")
    else:
        print("slow: compute pairwise FunUniFrac distances first.")
    sample_file = args.file
    brite = args.brite

    solver = EarthMoverDistanceUniFracSolver()
    tree = make_tree.import_graph(TREE_FILE, directed=True)
    FuncTreeEmduInput = make_emd_input.tree_to_EMDU_input(tree, brite)
    KO_index, vectors = extract_vectors_from_csv(sample_file)
    for P in vectors:
        P_pushed = solver.push_up_L1(P, FuncTreeEmduInput)
        print(P_pushed)
        break