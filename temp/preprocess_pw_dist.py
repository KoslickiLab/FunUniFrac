from Kegg_tree import get_KeggTree_from_edgelist
import time
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", help="edge list file with no branch lengths")
    parser.add_argument("-d", "--pw_dist", help="pairwise distance matrix file")
    parser.add_argument("-l", "--label_file", help="label for the pairwise distance matrix")
    parser.add_argument("-s", "--save", help="Save preprocessed distance file as")
    args = parser.parse_args()

    kegg_tree = get_KeggTree_from_edgelist(args.file, edge_length=False)
    kegg_tree.preprocess_pw_dist(args.pw_dist, args.label_file)
    kegg_tree.write_pw_dist(args.save)


if __name__ == "__main__":
    main()


