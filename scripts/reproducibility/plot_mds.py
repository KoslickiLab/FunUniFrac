#!/usr/bin/env python
import os, sys
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print(ROOT_DIR)
sys.path.append(ROOT_DIR)
sys.path.append(ROOT_DIR + '../src')
sys.path.append('../../')
from sklearn.manifold import MDS
import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path

def parsearg():
    parser = argparse.ArgumentParser(description="This script clusters the samples according to fununifrac distance"
                                                 "for a given .csv file")
    parser.add_argument("-pd", "--pairwise_distance", type=str, help="Path to pairwise distance file.")
    parser.add_argument("-l", "--label_file", type=str, help="Path to the label file")
    parser.add_argument('-m', '--metadata_file', help='Path to the metadata file.')
    return parser.parse_args()


if  __name__ == "__main__":
    args = parsearg()

    metadata = pd.read_table(args.metadata_file)
    with open(args.label_file, 'r') as f:
        labels = f.readlines()
        labels = [l.strip() for l in labels]
    labels = [Path(l).stem for l in labels]
    labels = [l.replace("sourmash_gather_out_scale1000_k_11_", "") for l in labels]
    metadata_dict = {x:y for (x, y) in zip(metadata["f_uid"], metadata["study_full_name"])}
    pw_dist = np.load(args.pairwise_distance)
    mds = MDS(dissimilarity='precomputed')
    coordinates = mds.fit_transform(pw_dist)
    colors = [metadata_dict[l] for l in labels]
    sns.scatterplot(x=coordinates[:, 0], y=coordinates[:, 1], hue=colors)
    plt.show()