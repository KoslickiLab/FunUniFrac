import os, sys, argparse
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(ROOT_DIR)
import multiprocessing
import logging
from src import make_all_pw_fununifrac


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This script will take a directory of sourmash gather results and "
                                                 "a weighted edge list representing the KEGG hierarchy and compute "
                                                 "all pairwise functional unifrac distances.")
    parser.add_argument('-e', '--edge_list',
                        help='Input edge list file of the KEGG hierarchy. Must have lengths in the '
                             'third column.', required=True)
    parser.add_argument('-fp', '--file_pattern', help="Pattern to match files in the directory. Default is "
                                                      "'*_gather.csv'", default='*_gather.csv')
    parser.add_argument('-o', '--out_file', help='Output file name: will be a numpy array.', required=True)
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
    make_all_pw_fununifrac.main(parser.parse_args())


