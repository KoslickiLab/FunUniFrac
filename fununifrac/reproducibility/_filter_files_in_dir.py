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
    ids = meta_data['f_uid']
    selected_files = [f'sourmash_gather_out_scale1000_k_11_{i}.csv' for i in ids]
    for i, file in enumerate(all_files):
        if os.path.basename(file) in selected_files:
            os.rename(f"{args.source_dir}/{file}", f"{args.target_dir}/{file}")

