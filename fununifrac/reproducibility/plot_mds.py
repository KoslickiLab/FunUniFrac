#!/usr/bin/env python
from sklearn.manifold import MDS
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path
from matplotlib.colors import ListedColormap

def parsearg():
    parser = argparse.ArgumentParser(description="This script clusters the samples according to fununifrac distance"
                                                 "for a given .csv file")
    parser.add_argument("-pd", "--pairwise_distance", type=str, help="Path to pairwise distance file.")
    parser.add_argument("-l", "--label_file", type=str, help="Path to the label file")
    parser.add_argument('-m', '--metadata_file', help='Path to the metadata file.')
    parser.add_argument('-t', '--title', help='Title of the plot.')
    parser.add_argument('-o', '--output', help='Output file name.')
    parser.add_argument('-c', '--condition', help='Column key for the condition based on which MDS is performed.',
                        default='study_full_name')
    parser.add_argument('-id', '--sample_id', help='Column name of the sample id.', default='f_uid')
    parser.add_argument('-factor', '--factor', help='Multiply everything in the input matrix by a constant.', type=int)
    return parser.parse_args()


if __name__ == "__main__":
    args = parsearg()
    if args.metadata_file.endswith('.csv'):
        metadata = pd.read_csv(args.metadata_file)
    else:
        metadata = pd.read_table(args.metadata_file)
    with open(args.label_file, 'r') as f:
        labels = f.readlines()
        labels = [l.strip() for l in labels]
    labels = [Path(l).stem for l in labels]
    labels = [l.replace("sourmash_gather_out_scale_1000_k_11_", "") for l in labels]
    labels = [l.split('.')[0] for l in labels]
    print(labels)
    metadata_dict = {x:y for (x, y) in zip(metadata[args.sample_id], metadata[args.condition])}
    pw_dist = np.load(args.pairwise_distance)

    #filter according to metadata
    if len(metadata_dict) < len(labels): #not enough metadata
        remove_labels = []
        remove_index = []
        for i, l in enumerate(labels):
            if l not in metadata_dict:
                remove_labels.append(l)
                remove_index.append(i)
        labels = [l for l in labels if l not in remove_labels]
        print(f"Remove {len(remove_labels)} samples")
        pw_dist = np.delete(pw_dist, remove_index, axis=0)
        pw_dist = np.delete(pw_dist, remove_index, axis=1)

    if args.factor:
        pw_dist = pw_dist*args.factor
    fig = plt.figure()
    #ax = fig.add_subplot(projection='3d')
    mds = MDS(dissimilarity='precomputed', n_components=2)
    coordinates = mds.fit_transform(pw_dist)
    colors = [metadata_dict[l] for l in labels]
    sns.scatterplot(x=coordinates[:, 0], y=coordinates[:, 1], hue=colors)
    #cmap = ListedColormap(sns.color_palette("husl", 256).as_hex())
    #ax.scatter(coordinates[:, 0], coordinates[:, 1], coordinates[:, 2], cmap=cmap, alpha=1)
    if args.title:
        plt.title(args.title)
    plt.savefig(args.output)
    plt.show()