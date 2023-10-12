import os, sys, argparse
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(ROOT_DIR)
import logging
from src.factory.make_emd_input import functional_profile_to_vector

def main():
    parser = argparse.ArgumentParser(description="Given a group of samples, compute the average.")
    parser.add_argument('-e', '--edge_list',
                        help='Input edge list file of the KEGG hierarchy. Must have lengths in the '
                             'third column.', required=True)
    parser.add_argument('-fd', '--file_dir', help="Directory of sourmash files.", required=True)
    parser.add_argument('-fp', '--file_pattern', help="Pattern to match files in the directory. Default is "
                                                      "'*_gather.csv'", default='*_gather.csv')
    parser.add_argument('-o', '--out_dir', help='Output directory name.', required=True)
    parser.add_argument('-i', '--out_id', help='Test purpose: give an identifier to the output file so that tester can recognize it')
    parser.add_argument('-f', '--force', help='Overwrite the output file if it exists', action='store_true')
    parser.add_argument('-a', '--abundance_key',
                        help='Key in the gather results to use for abundance. Default is `f_unique_weighted`',
                        default='f_unique_weighted')
    parser.add_argument('--diffab', action='store_true', help='Also return the difference abundance vectors.')
    parser.add_argument('-v', help="Be verbose", action="store_const", dest="loglevel", const=logging.INFO,
                        default=logging.WARNING)
    parser.add_argument('--unweighted', help="Compute unweighted unifrac instead of the default weighted version",
                        action="store_true")
    parser.add_argument('--L2', help="Use L2 UniFrac instead of L1", action="store_true")
    parser.add_argument('--Ppushed', help="Flag indicating you want the pushed vectors to be saved.",
                        action="store_true")

    #take the average directly
