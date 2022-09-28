#!/usr/bin/env python
import argparse
import os
import sys
sys.path.append('./src')
import numpy as np
from LP_EMD_helper import get_matrix_from_edge_list, parse_edge_list



parser = argparse.ArgumentParser(description='Given an edge list file, convert it into a sparse matrix to be used'
                                             'as an input for LP_EMD.py')
parser.add_argument('-e', '--edge_list', help='Input edge list file of the KEGG hierarchy', required=True)
parser.add_argument('-o', '--out_dir', help='Output directory', required=True)
parser.add_argument('--force', help='Overwrite the output file if it exists', action='store_true')
args = parser.parse_args()

edge_list = args.edge_list
out_dir = args.out_dir
force = args.force
if not os.path.exists(edge_list):
    raise FileNotFoundError(f"Could not find {edge_list}")

basename = os.path.splitext(os.path.basename(edge_list))[0]
file_name = f"{basename}_A.npy"
full_file_name = os.path.join(out_dir, file_name)
if os.path.exists(full_file_name) and not force:
    raise FileExistsError(f"{full_file_name} already exists. Please delete it or choose another name, or use --force.")

df = parse_edge_list(edge_list)
A, node_list = get_matrix_from_edge_list(df)
np.save(full_file_name, A)
print(f"{file_name} saved")
