import argparse
import pandas as pd
import sys

EPSILON = sys.float_info.epsilon


def parsearg():
    parser = argparse.ArgumentParser(description="Temporary script. Fix an edge list by giving a very small value"
                                                 "to any branch that has a negative value.")
    parser.add_argument('-e', '--edge_list',
                        help='Input edge list file of the KEGG hierarchy. Must have lengths in the '
                             'third column.', required=True)
    parser.add_argument('-s', '--save', help='Path to save the output file.', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parsearg()
    df = pd.read_table(args.edge_list)
    for (index, row) in df.iterrows():
        if row['edge_length'] <= 0:
            df.at[index, 'edge_length'] = EPSILON
    print(df)
