import matplotlib.pyplot as plt
import numpy as np
import argparse


def parsearg():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-pd", "--pairwise_distance", type=str, help="Path to pairwise distance file.")
    parser.add_argument('-o', '--output', help='Output file name.')

    return parser.parse_args()


if __name__ == "__main__":
    args = parsearg()
    pw_distance = np.load(args.pairwise_distance, allow_pickle=True)
    data = pw_distance.flatten()
    print(data.max())
    print(data.min())
    print(data.mean())
    plt.hist(data, bins=50)
    plt.savefig(args.output)
    plt.show()