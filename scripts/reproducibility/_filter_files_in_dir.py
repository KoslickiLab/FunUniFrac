import pandas as pd
import argparse
import os


def parsearg():
    parser = argparse.ArgumentParser(description="Creates a new directory and move files in source directory "
                                                 "to the target directory based on filtering conditions. For the "
                                                 "running of compute_fununifrac.py.")
    parser.add_argument('-m', '--metadata_file', help='Path to the metadata file.')
    parser.add_argument('-sd', '--source_dir', help='Path to the source directory where original files are.')
    parser.add_argument('-td', '--target_dir', help='Path to the target directory where original files are.')
    return parser.parse_args()


if __name__ == "__main__":
    args = parsearg()
    all_files = os.listdir(args.source_dir)
    meta_data = pd.read_table(args.metadata_file)
    print(all_files[:5])
    #for (index, row) in meta_data.iterrows():
    #    file
