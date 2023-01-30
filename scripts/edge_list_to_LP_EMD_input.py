#!/usr/bin/env python
import argparse
import os
import numpy as np
from LP_EMD_helper import get_distance_matrix_from_edge_list, parse_edge_list


def argument_parser():
    parser = argparse.ArgumentParser(description='Given an edge list file, convert it into a sparse matrix to be used'
                                                 'as an input for LP_EMD.py')
    parser.add_argument('-e', '--edge_list', help='Input edge list file of the KEGG hierarchy', required=True)
    parser.add_argument('-o', '--out_dir', help='Output directory', required=True)
    parser.add_argument('--force', help='Overwrite the output file if it exists', action='store_true')
    return parser


def main():
    parser = argument_parser()
    args = parser.parse_args()

    edge_list_file = args.edge_list
    out_dir = args.out_dir
    force = args.force
    if not os.path.exists(edge_list_file):
        raise FileNotFoundError(f"Could not find {edge_list_file}")

    basename = os.path.splitext(os.path.basename(edge_list_file))[0]
    file_name = f"{basename}_D.npy"
    full_file_name = os.path.join(out_dir, file_name)
    if os.path.exists(full_file_name) and not force:
        raise FileExistsError(f"{full_file_name} already exists. Please delete it or choose another name, or use --force.")

    A, node_list = get_distance_matrix_from_edge_list(edge_list_file)
    np.save(full_file_name, A)
    print(f"{file_name} saved")


if __name__ == '__main__':
    main()