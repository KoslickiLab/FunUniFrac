import scipy.spatial as sp, scipy.cluster.hierarchy as hc
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
    parser.add_argument('-t', '--title', help='Title of the plot.')
    parser.add_argument('-o', '--output', help='Output file name.')
    return parser.parse_args()


if __name__ == "__main__":
    args = parsearg()
    pw_distance = np.load(args.pairwise_distance)
    print(pw_distance)
    linkage = hc.linkage(pw_distance, method='average')
    metadata = pd.read_table(args.metadata_file)

    lut = dict(zip(metadata['study_full_name'].unique(), "rbg"))
    row_colors = metadata['study_full_name'].map(lut)
    sns.clustermap(pw_distance, row_linkage=linkage, col_linkage=linkage)
    plt.show()


    # from sklearn.datasets import load_iris
    # iris = load_iris()
    # X, y = iris.data, iris.target
    # DF = pd.DataFrame(X, index=["iris_%d" % (i) for i in range(X.shape[0])], columns=iris.feature_names)
    # DF_corr = DF.T.corr()
    # DF_dism = 1 - DF_corr  # distance matrix
    # linkage = hc.linkage(sp.distance.squareform(DF_dism), method='average')
    # sns.clustermap(DF_dism, row_linkage=linkage, col_linkage=linkage)
