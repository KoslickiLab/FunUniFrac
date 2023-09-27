from Kegg_tree import visualize_diff
import argparse


def main():
    parser = argparse.ArgumentParser(description="Compare the branch lengths between two edge lists"
                                                 "the two files should have the same branches in the same"
                                                 "order.")
    parser.add_argument('-i', '--inferred', help='edge list file of computed branch lengths.')
    parser.add_argument('-r', '--reference', help='edge list file of reference. Branches should be identical as that '
                                                  'of the inferred file.')
    parser.add_argument('-o', '--out_file', help='Output file name. Should be .png file.')
    args = parser.parse_args()

    visualize_diff(args.inferred, args.reference, args.out_file)


if __name__ == '__main__':
    main()