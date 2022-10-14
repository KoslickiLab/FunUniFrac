#!/usr/bin/env python
# explore the pairwise distances. This script is intended to be run interactively
import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
import argparse
matplotlib.use('Agg')

def argument_parser():
    parser = argparse.ArgumentParser(
        description='Given a matrix of pairwise FunUnifrac distances along with metadata associated '
                    'with the study, this script will plot the distances in a dendrogram. It is assumed that the'
                    'metadata file has a column called "sample_name" which contains the sample names which correspond'
                    'to the samples in the distance matrix. The distance matrix is assumed to be a numpy array')
    parser.add_argument('-d', '--distances', help='File of the pairwise unifrac distances', required=True)
    parser.add_argument('-o', '--out_dir', help='Output directory', required=True)
    parser.add_argument('-m', '--metadata', help='Metadata file', required=True)
    parser.add_argument('-k', '--key', help='Metadata key to use to label dendrogram leaves', required=False, default='disease_type')
    return parser


def main():
    args = argument_parser().parse_args()
    pairwise_dists_file = args.distances
    out_dir = args.out_dir
    meta_data_file = args.metadata
    key = args.key
    # check if files exist
    if not os.path.exists(pairwise_dists_file):
        raise FileNotFoundError(f"Could not find {pairwise_dists_file}")
    if not os.path.exists(meta_data_file):
        raise FileNotFoundError(f"Could not find {meta_data_file}")
    # check if basis file exists
    basis_file = f"{pairwise_dists_file}.basis.txt"
    if not os.path.exists(basis_file):
        raise FileNotFoundError(f"Could not find {basis_file}")
    # load the metadata
    meta_data = pd.read_csv(meta_data_file, sep='\t')
    # load the distances
    pairwise_dists = np.load(pairwise_dists_file)
    # scipy expects a condensed distance matrix (i.e. not the full matrix)
    pairwise_dists_cond = squareform(pairwise_dists)
    # load the basis
    with open(basis_file, 'r') as f:
        basis = [os.path.basename(line.strip()) for line in f]
    # get the sample name from the file names
    sample_names = [os.path.basename(name).split('.')[0] for name in basis]
    prefix = meta_data['sample_name'][0].split('.')[0]  # get the prefix of the sample names. eg "13984."
    prefix = prefix + "."
    sample_names_with_prefix = [prefix + x for x in sample_names]
    print(f"sample_names_with_prefix: {sample_names_with_prefix}")
    # get the metadata for the samples
    #meta_data_of_dist = meta_data[meta_data['sample_name'].isin(sample_names_with_prefix)]
    # cluster the pairwise distances
    Z = linkage(pairwise_dists_cond, method='ward')
    # plot the dendrogram
    plt.figure(figsize=(10, 40), dpi=600)
    plt.title(f"Hierarchical Clustering of\n {pairwise_dists_file}")
    plt.xlabel('sample index')
    plt.ylabel('distance')
    #dendrogram(Z, labels=sample_names_with_prefix)
    labels = [meta_data.loc[meta_data['sample_name'] == sample_names_with_prefix[x]][key].values[0] for x in
              range(len(sample_names_with_prefix))]
    dendrogram(Z, labels=labels, orientation='left', leaf_font_size=7)
    # make the labels bigger
    plt.rcParams.update({'font.size': 22})
    plt.savefig(os.path.join(out_dir,f"dendro_{pairwise_dists_file}.png"), dpi=600)


if __name__ == '__main__':
    main()
