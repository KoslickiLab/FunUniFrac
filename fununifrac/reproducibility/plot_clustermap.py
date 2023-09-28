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
    print(pw_distance.shape)
    linkage = hc.linkage(pw_distance, method='average')
    metadata = pd.read_table(args.metadata_file)

    with open(args.label_file, 'r') as f:
        labels = f.readlines()
        labels = [l.strip() for l in labels]
    labels = [Path(l).stem for l in labels]
    labels = [l.replace("sourmash_gather_out_scale1000_k_11_", "") for l in labels]
    metadata_dict = {x: y for (x, y) in zip(metadata["f_uid"], metadata["study_full_name"])}

    #phenotypes = metadata.pop('study_full_name')
    phenotypes = [metadata_dict[i] for i in labels]
    lut = dict(zip(['IBD', 'T2D', 'HHS'], "rbg"))
    #row_colors = phenotypes.map(lut)
    row_colors = [lut[i] for i in phenotypes]


    #sns.clustermap(pw_distance, row_colors=row_colors, row_linkage=linkage, col_linkage=linkage)
    sns.clustermap(pd.DataFrame(pw_distance, columns=phenotypes), row_linkage=linkage, col_linkage=linkage, row_colors=row_colors)
    # iris = sns.load_dataset("iris")
    # species = iris.pop("species")
    # sns.clustermap(iris)
    # lut = dict(zip(species.unique(), "rbg"))
    # row_colors = species.map(lut)
    # sns.clustermap(iris, row_colors=row_colors)
    plt.savefig(args.output)
    plt.show()

