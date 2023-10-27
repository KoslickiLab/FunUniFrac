import argparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path
from skbio import DistanceMatrix #to install: pip install scikit-bio
from skbio.stats.ordination import pcoa
from scipy.spatial.distance import pdist


def parsearg():
    parser = argparse.ArgumentParser(description="This script clusters the samples according to fununifrac distance"
                                                 "for a given .csv file")
    parser.add_argument("-f", "--file", type=str, help="File containing vectors.", required=True)
    parser.add_argument('-m', '--metadata_file', help='Path to the metadata file.')
    parser.add_argument('-t', '--title', help='Title of the plot.')
    parser.add_argument('-o', '--output', help='Output file name.')
    parser.add_argument('-c', '--condition', help='Column key for the condition based on which MDS is performed.',
                        default='study_full_name')
    parser.add_argument('-id', '--sample_id', help='Column name of the sample id.', default='f_uid')
    parser.add_argument('-2D', '--two_D', action='store_true')
    return parser.parse_args()


if __name__ == "__main__":
    args = parsearg()

    df = pd.read_csv(args.file, index_col=0)
    labels = df.columns
    print(labels)
    if args.metadata_file.endswith('.csv'):
        metadata = pd.read_csv(args.metadata_file, index_col=0)
    else:
        metadata = pd.read_table(args.metadata_file, index_col=0)
    metadata_dict = {x: y for (x, y) in zip(metadata[args.sample_id], metadata[args.condition])}
    # filter according to metadata
    if len(metadata_dict) < len(labels):  # not enough metadata
        remove_labels = []
        for l in labels:
            if l not in metadata_dict:
                remove_labels.append(l)
        print(f"df shape before remove: {df.shape}")
        for l in remove_labels:
            df.drop(l, axis=1, inplace=True)
        print(f"df.shape after remove: {df.shape}")
        labels = [l for l in labels if l not in remove_labels]
        print(f"Remove {len(remove_labels)} samples")

    filtered_meta_df = metadata[metadata[args.sample_id].isin(labels)]
    filtered_meta_df.set_index(args.sample_id, inplace=True)
    pw_dist = pdist(df.transpose()) #rows:samples, columns:KO
    dm = DistanceMatrix(pw_dist, labels)
    dist_pc = pcoa(dm)
    if args.two_D:
        ordination = dist_pc.samples[['PC1', 'PC2']]
        merged_df = ordination.merge(metadata, left_index=True, right_on=args.sample_id, how="left")
        cleaned_merged_df = merged_df[merged_df['PC2'] < 0.01]
        print(cleaned_merged_df.shape)
        fig = sns.scatterplot(data=cleaned_merged_df, x="PC1", y="PC2", hue=args.condition)
        fig = fig.get_figure()
        if args.output:
            plt.savefig(args.output)
    else:
        dist_pc.plot(df=filtered_meta_df, column=args.condition, cmap='Set1', axis_labels=('PC1', 'PC2', 'PC3'))
        if args.output:
            plt.savefig(args.output)
    plt.show()
