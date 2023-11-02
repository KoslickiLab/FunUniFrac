import argparse
import pandas as pd
import random


def main():
    parser = argparse.ArgumentParser(description="Assign branch lengths randomly to a given tree.")
    parser.add_argument('-e', '--edge_list', help='Input edge list file of the KEGG hierarchy.',
                        default='data/kegg_trees/kegg_ko00001_no_edge_lengths_filtered.txt')
    parser.add_argument('-o', '--out_file', help='Path to save the output file.', default='data/kegg_trees/kegg_tree_random.txt')

    args = parser.parse_args()
    df = pd.read_table(args.edge_list)
    df['edge_length'] = [random.random() for _ in range(len(df['child']))]
    df.to_csv(args.out_file)


main()