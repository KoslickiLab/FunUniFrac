import matplotlib.pyplot as plt
import numpy as np
import argparse
import pandas as pd


def parsearg():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-pd", "--pairwise_distance", type=str, help="Path to pairwise distance file.")
    parser.add_argument('-o', '--output', help='Output file name.')
    parser.add_argument('-c', '--column', help='Column of the dataframe of which the histogram is computed.')
    parser.add_argument('-f', '--file', help='A dataframe. Most likely ending in .txt')
    return parser.parse_args()


if __name__ == "__main__":
    args = parsearg()
    if args.pairwise_distance:
        pw_distance = np.load(args.pairwise_distance, allow_pickle=True)
        data = pw_distance.flatten()
    elif args.file:
        df = pd.read_table(args.file)
        data = df[args.column]
    else:
        data=None
    plt.hist(data, bins=50)
    plt.savefig(args.output)
    plt.show()


