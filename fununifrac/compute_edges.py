import os, sys, argparse
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(ROOT_DIR)
from src import compute_func_edge_lengths


def main():
    """
    This function takes the matrix output of graph_to_path_matrix.py as well as the pairwise distance matrix, and solve
    the least squares problem of inferring the edge lengths of the graph.

    :return: an edge list file with branch lengths, ouput into user-specified directory.
    """
    parser = argparse.ArgumentParser(
        description='This will take the matrix made by graph_to_path_matrix.py and the all '
                    'pairwise distance matrix and solve the least squares problem of '
                    'inferring the edge lengths of the graph.')
    parser.add_argument('-e', '--edge_file', help='Input edge list file of the KEGG hierarchy', required=True)
    parser.add_argument('-d', '--distance_file', help='File containing all pairwise distances between KOs. Use sourmash '
                                                  'compare', required=True)
    parser.add_argument('-o', '--out_dir', help='Output directory: the location to place the output file with edge list with lengths in the last column',
                        required=True)
    parser.add_argument('-i', '--out_id', help='Test purpose: give an identifier to the output file so that tester can recognize it')
    parser.add_argument('-b', '--brite_id', help='Brite ID of the KEGG hierarchy you want to focus on. Eg. ko00001',
                        required=True)
    parser.add_argument('-n', '--num_iter', help='Number of random selections on which to perform the NNLS',
                        default=100)
    parser.add_argument('-f', '--factor', help='Selects <--factor>*(A.shape[1]) rows for which to do the NNLS',
                        default=5)
    parser.add_argument('-r', '--reg_factor', help='Regularization factor for the NNLS', default=1)
    parser.add_argument('--distance',
                        help='Flag indicating that the input matrix is a distance (0=identical). If not, it is assumed to be a similarity (1=identical).',
                        action='store_true')
    # call main function
    compute_func_edge_lengths.main(parser.parse_args())

    
if __name__ == '__main__':
    main()
