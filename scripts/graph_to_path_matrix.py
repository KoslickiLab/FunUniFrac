import os, sys, argparse
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(ROOT_DIR)
from src import graph_to_path_matrix


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Given the KEGG hierarchy, first this will select the subtree consisting of'
                    ' all the ancestors of the given brite_id. Then, it will create a matrix '
                    'where the (i,j) entry will be 1 iff for the ith pair of KOs in the '
                    '--distances file: (KO1, KO2), edge j is on the shortest path from '
                    'KO1 to KO2')
    parser.add_argument('-e', '--edge_list', help='Input edge list file of the KEGG hierarchy', required=True)
    parser.add_argument('-d', '--distances', help='File containing all pairwise distances between KOs. Use sourmash '
                                                  'compare', required=True)
    parser.add_argument('-o', '--out_dir', help='Output directory', required=True)
    parser.add_argument('-b', '--brite_id', help='Brite ID of the KEGG hierarchy you want to focus on. Eg. ko00001',
                        required=True)
    # call main function
    graph_to_path_matrix.main(parser.parse_args())
