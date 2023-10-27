import seaborn as sns
import argparse
import pandas as pd
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description='Plot box plot based on given  dataframe')
    parser.add_argument("-df", "--dataframe", type=str, help="Path to pairwise distance file.")
    parser.add_argument("-x", "--x", type=str, help="X coordinate of the plot")
    parser.add_argument("-y", "--y", type=str, help="Y coordinate of the plot")
    parser.add_argument("-t", "--title", type=str, help="Title of the plot")
    parser.add_argument("-hue", "--hue", type=str, help="Hue of the plot")
    parser.add_argument("-o", "--output", type=str, help="Output file name")
    parser.add_argument("-b", "--box", type=str, action="store_true", help="If given boxplot instead of line")

    args = parser.parse_args()
    df = pd.read_table(args.dataframe)
    print(df)
    if args.box:
        sns.boxplot(df, x=args.x, y=args.y, hue=args.hue)
    else:
        sns.lineplot(df, x=args.x, y=args.y, hue=args.hue)
    if args.title:
        plt.title(args.title)
    plt.savefig(args.output)
    plt.show()

main()