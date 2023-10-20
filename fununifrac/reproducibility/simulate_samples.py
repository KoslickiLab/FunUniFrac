import argparse
import pandas as pd
import random
import numpy as np

COL_NUM = 100


#get hold of a template
def parsearg():
    parser = argparse.ArgumentParser(description="This script generates simulated data")
    parser.add_argument("-f", "--template_file", type=str, help="Template file")
    parser.add_argument("-o", "--out_dir", type=str, help="Output directory")
    return parser.parse_args()


def main():
    args = parsearg()

    df = pd.read_csv(args.template_file, index_col='name')
    print(f"No. of samples: {len(df.columns)}")
    print(f"No. of KO: {len(df.index)}")
    df = df[df.columns[:COL_NUM]]
    meta_dict = dict()
    for col in df.columns[:50]:
        meta_dict[col] = "environment 1"
    for col in df.columns[50:]:
        meta_dict[col] = "environment 2"
    meta_df = pd.DataFrame(meta_dict.items(), columns=['sample', 'env'])
    meta_df.to_csv(f"{args.out_dir}/simulated_metadata.csv")

    sim_dict = {
        0.1: 'low',
        0.5: 'medium',
        0.9: 'high',
    }
    for i in range(100):
        for proportion in sim_dict:
            partition = int(len(df.index)*proportion)
            env1_distribution_vector = np.ones(len(df.index))
            env1_distribution_vector[:partition] = np.random.exponential(scale=50, size=partition)
            env2_distribution_vector = np.ones(len(df.index))
            env2_distribution_vector[partition:] = np.random.exponential(scale=50, size=partition)
            file_name = f"{args.out_dir}/sim_sample_{sim_dict[proportion]}_{i}.csv"
            for col in df.columns[:50]:
                df.columns[col] = np.random.dirichlet(env1_distribution_vector, 1)[0]
            for col in df.columns[50:]:
                df.columns[col] = np.random.dirichlet(env2_distribution_vector, 1)[0]
            df.to_csv(file_name)


if __name__ == '__main__':
    main()
