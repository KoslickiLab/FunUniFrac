import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import itertools as it
from Kegg_tree import L1_norm


def main():
    parser = argparse.ArgumentParser(description="Compare pw distances between two pw matrix files.")
    parser.add_argument('-r', '--recovered', help='edge list file of computed branch lengths.')
    parser.add_argument('-i', '--original', help='Original pw matrix file.')
    parser.add_argument('-o', '--out_file', help='Output file name. Should be .png file.')
    parser.add_argument('-lo', '--label_original', help='Labels for file1.')
    parser.add_argument('-lr', '--label_recovered', help='Labels for file2.')

    args = parser.parse_args()
    pw_dist_original = np.load(args.original)
    pw_dist_recovered = np.load(args.recovered)
    labels_original = [line.strip() for line in open(args.label_original, 'r')]
    labels_recovered = [line.strip() for line in open(args.label_recovered, 'r')]
    labels_original_index_dict = {l:i for i, l in enumerate(labels_original)}
    labels_recovered_index_dict = {l:i for i, l in enumerate(labels_recovered)}

    branches = [(i, j) for (i, j) in it.combinations(labels_recovered, 2)]

    df = pd.DataFrame(columns=['branch', 'Original_distance', 'Reconstructed_distance'])
    original_dist = []
    recovered_dist = []
    for (node1, node2) in branches:
        original_dist.append(pw_dist_original[labels_original_index_dict[node1]][labels_original_index_dict[node2]])
        recovered_dist.append(pw_dist_recovered[labels_recovered_index_dict[node1]][labels_recovered_index_dict[node2]])
    df['branch'] = branches
    df['Reconstructed_distance'] = recovered_dist
    df['Original_distance'] = original_dist
    print(L1_norm(df['Reconstructed_distance'], df['Original_distance']))
    sns.scatterplot(data=df, x='Original_distance', y='Reconstructed_distance', s=1)
    plt.savefig(args.out_file)


if __name__ == '__main__':
    main()
