import argparse
import pandas as pd


def parsearg():
    parser = argparse.ArgumentParser(description="This script is produced because some edge lists have KOs that lack "
                                                 "necessary information for its edge length to be computed. As such, "
                                                 "we will have to remove those nodes from the original edge list. It's "
                                                 "sad, but I don't have a better way yet.")
    parser.add_argument('-e', '--edge_list',
                        help='Input edge list file of the KEGG hierarchy. Must have lengths in the '
                             'third column.', required=True)
    parser.add_argument('-l', '--label_file', help='Label file for the distance matrix. KOs in edge list but not in '
                                                   'the label file will be removed.', required=True)
    parser.add_argument('-s', '--save', help='Path to save the output file.', required=True)
    return parser.parse_args()


if __name__ == "__main__":
    args = parsearg()
    df = pd.read_table(args.edge_list)
    with open(args.label_file, 'r') as f:
        KOs = f.readlines()
        KOs = [k.strip() for k in KOs]
    print(f"{len(KOs)} in label file")
    for (index, row) in df.iterrows():
        if row['child'].startswith('K') and row['child'] not in KOs:
            print(f"remove {row['child']}")
            df.drop(index, axis=0, inplace=True)
    df.to_csv(args.save, sep='\t', index=False)




