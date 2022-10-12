#!/usr/bin/env python
import argparse
import os
import sys
import numpy as np
try:
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    sys.path.append(os.path.dirname(SCRIPT_DIR))
except:
    pass
from src.LP_EMD_helper import get_profile_from_sourmash

def argument_parser():
    parser = argparse.ArgumentParser(description="Given a sourmash gather file, extract into a relative abundance vector"
                                                 "and write into a .txt file.")
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input file, a sourmash gather file.')
    parser.add_argument('-od', '--out_dir', type=str, help='Path to the output directory')
    parser.add_argument('-p','--prefix', type=str, help='Prefix of the output file. If not given, will be the same as the input file.')
    return parser

def main():
    parser = argument_parser()
    args = parser.parse_args()
    sourmash_file = args.input
    if not args.prefix:
        #get basename of input
        prefix = os.path.splitext(os.path.basename(sourmash_file))[0]
    output_file = os.path.join(args.out_dir, prefix + ".txt")
    get_profile_from_sourmash(sourmash_file, output_file, normalize=True)


if __name__ == '__main__':
    main()