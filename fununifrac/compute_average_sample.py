import os, sys, argparse
import pandas as pd
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(ROOT_DIR)
import numpy as np
from src.utility.merge_fun_files import merge_files
import re


def normalize(vector):
    norm = vector/np.sum(vector)
    return norm

def L2_average(list_of_vectors):
    return np.mean(list_of_vectors, axis=0)

def main():
    parser = argparse.ArgumentParser(description="Given a group of samples, compute the average.")
    parser.add_argument('-fd', '--file_dir', help="Directory of sourmash files.", required=True)
    parser.add_argument('-f', '--file', type=str,help='File containing the combined profiles.')
    parser.add_argument('-fp', '--file_pattern', help="Pattern to match files in the directory. Default is "
                                                      "'*_gather.csv'", default='*_gather.csv')
    parser.add_argument('-o', '--out_file', help='Output file name.', required=True)
    parser.add_argument('-a', '--abundance_key',
                        help='Key in the gather results to use for abundance. Default is `f_unique_weighted`',
                        default='f_unique_weighted')
    parser.add_argument('-m', '--metadata', help='Metadata file')
    parser.add_argument('-c', '--condition', help='Column key for the condition based on which averaging is performed.',
                        default='study_full_name')
    parser.add_argument('-id', '--sample_id', help='Column name of the sample id.', default='f_uid')

    #take the average directly
    args = parser.parse_args()
    if args.file:
        df = pd.read_csv(args.file)
    else:
        df = merge_files(args.file_dir, args.file_pattern, args.abundance_key, 'name')

    df.fillna(0, inplace=True)
    if args.metadata.endswith('.csv'):
        meta_df = pd.read_csv(args.metadata)
    else:
        meta_df = pd.read_table(args.metadata)
    print(meta_df)
    meta_dict = dict(zip(meta_df[args.sample_id], meta_df[args.condition]))
    conditions = set(meta_df[args.condition])
    condition_dict = dict()
    for c in conditions:
        pattern = r'(SRR|DRR|ERR)\d+'
        df.columns = [re.search(pattern, b).group() for b in df.columns]
        samples = [df[sample_id].to_list() for sample_id in df.columns if meta_dict[sample_id] == c]
        print(f"condition {c} has {len(samples)} samples")
        condition_dict[c] = samples

    average_vector_df = pd.DataFrame(index=df.index, columns=list(condition_dict.keys()))
    for c in condition_dict:
        normalized_vectors = [normalize(v) for v in condition_dict[c]]
        average_vector = L2_average(normalized_vectors)
        average_vector_df[c] = average_vector
    average_vector_df.to_csv(args.out_file, sep='\t')


if __name__ == '__main__':
    main()