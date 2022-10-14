#!/usr/bin/env python
# explore the pairwise distances. This script is intended to be run interactively
import sys

import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, fcluster
from scipy.spatial.distance import squareform
import argparse
import json
matplotlib.use('Agg')

def argument_parser():
    parser = argparse.ArgumentParser(
        description='Given a matrix of pairwise FunUnifrac distances along with metadata associated '
                    'with the study, this script will plot the distances in a dendrogram. It is assumed that the'
                    'metadata file has a column called "sample_name" which contains the sample names which correspond'
                    'to the samples in the distance matrix. The distance matrix is assumed to be a numpy array.'
                    ''
                    'This script will also return a JSON object with the following information: '
                    'As we move a cutoff threshold from left to right (top to bottom) of the dendrogram,'
                    'keep track of how the underlying samples cluster. Eg. just below the root, the dendrogram splits'
                    'into 2 subtrees. Thus the values of this JSON object[2] will be two clusters of sample names'
                    'corresponding to these two groups. ')
    parser.add_argument('-d', '--distances', help='File of the pairwise unifrac distances', required=True)
    parser.add_argument('-o', '--out_dir', help='Output directory', required=True)
    parser.add_argument('-m', '--metadata', help='Metadata file', required=True)
    parser.add_argument('-k', '--key', help='Metadata key to use to label dendrogram leaves', required=False,
                        default='disease_type')
    parser.add_argument('-s', '--split', help='The number of splits to retain (i.e. how far right to move the threshold'
                                              'when saving the JSON recording of clusters).', required=False,
                        default=5, type=int)
    return parser


def main():
    args = argument_parser().parse_args()
    pairwise_dists_file = args.distances
    out_dir = args.out_dir
    meta_data_file = args.metadata
    key = args.key
    split = args.split
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
    # get the metadata for the samples
    #meta_data_of_dist = meta_data[meta_data['sample_name'].isin(sample_names_with_prefix)]
    # cluster the pairwise distances
    Z = linkage(pairwise_dists_cond, method='ward')
    # plot the dendrogram
    plt.figure(figsize=(20, 40), dpi=600)
    plt.title(f"Hierarchical Clustering of\n {pairwise_dists_file}")
    plt.xlabel('sample index')
    plt.ylabel('distance')
    #dendrogram(Z, labels=sample_names_with_prefix)
    labels = [meta_data.loc[meta_data['sample_name'] == sample_names_with_prefix[x]][key].values[0] for x in
              range(len(sample_names_with_prefix))]
    dendrogram(Z, labels=labels, orientation='left', leaf_font_size=7)
    # make the labels bigger
    plt.rcParams.update({'font.size': 22})
    plt.savefig(os.path.join(out_dir, f"dendro_{pairwise_dists_file}.png"), dpi=600)

    # prepare the json dump of the clusters as we move the cutoff to the right
    # These are the values on the x-axis of the dendrogram where clusters are merged
    cuts = Z[:, 2]
    # The midpoints between the cuts
    cut_mids = (cuts[1:] + cuts[:-1]) / 2
    # This is a dictionary going from the top/left of the dendrogram to the bottom/right
    # i.e. most clustered to least clustered
    cut_ind_to_mids = dict()
    for i, cut in enumerate(cut_mids[::-1]):
        cut_ind_to_mids[i] = cut

    cut_to_clusters = dict()
    for cut_ind in range(min(split, len(cut_mids))):
        mid = cut_ind_to_mids[cut_ind]
        cut_to_clusters[int(cut_ind)] = dict()
        clust_inds = fcluster(Z, mid, criterion='distance')
        clusts = dict()
        for i, sample in zip(clust_inds, basis):
            if i not in clusts:
                clusts[int(i)] = []
            clusts[int(i)].append(sample)
        cut_to_clusters[int(cut_ind)] = clusts

    # write the clusters to a json file
    with open(os.path.join(out_dir, f"dendro_{pairwise_dists_file}_clusters.json"), 'w') as f:
        json.dump(cut_to_clusters, f, indent=4)


if __name__ == '__main__':
    main()


