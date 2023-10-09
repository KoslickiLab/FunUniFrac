import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path
from skbio import DistanceMatrix #to install: pip install scikit-bio
from skbio.stats.ordination import pcoa


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
    pw_dist = np.load(args.pairwise_distance)
    dm = DistanceMatrix(pw_dist, labels)
    filtered_meta_df = metadata[metadata[args.sample_id].isin(labels)]
    filtered_meta_df.set_index(args.sample_id, inplace=True)
    dist_pc = pcoa(dm)
    dist_pc.plot(df=filtered_meta_df, column=args.condition, cmap='Set1', axis_labels=('PC1', 'PC2', 'PC3'))
    plt.show()
