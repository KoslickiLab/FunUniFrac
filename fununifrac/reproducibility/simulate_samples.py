import argparse
import pandas as pd
import numpy as np

COL_NUM = 100
sim_dict = {
    0.1: '10',
    0.5: '50',
    0.9: '90',
    0.95: '95',
    0.99: '99',
}
mean_pos = [6000, 7000, 8000, 8500]
peak_heights = [100, 250, 500, 1000]

similarity_levels = {
    6000: 'low',
    7000: 'medium',
    8000: 'high',
    8500: 'very_high',
}


#get hold of a template
def parsearg():
    parser = argparse.ArgumentParser(description="This script generates simulated data")
    parser.add_argument("-f", "--template_file", type=str, help="Template file")
    parser.add_argument("-o", "--out_dir", type=str, help="Output directory")
    parser.add_argument("-d", "--distribution", type=str, choices=['normal', 'exponential'], default='exponential')
    return parser.parse_args()

def univariate_normal(x, mean, variance):
    """pdf of the univariate normal distribution."""
    return ((1. / np.sqrt(2 * np.pi * variance)) *
            np.exp(-(x - mean)**2 / (2 * variance)))


def main():
    args = parsearg()

    df = pd.read_csv(args.template_file, index_col='name')
    df = df[df.columns[:COL_NUM]]
    meta_dict = dict()
    for col in df.columns[:50]:
        meta_dict[col] = "environment 1"
    for col in df.columns[50:]:
        meta_dict[col] = "environment 2"
    meta_df = pd.DataFrame(meta_dict.items(), columns=['sample', 'env'])
    meta_df.to_csv(f"{args.out_dir}/simulated_metadata.csv")

    for i in range(50):
        if args.distribution == 'exponential':
            for proportion in sim_dict:
                partition = int(len(df.index)*proportion)
                env1_distribution_vector = np.ones(len(df.index))
                env1_distribution_vector[:partition] = list(range(partition, 0, -1))
                env2_distribution_vector = np.ones(len(df.index))
                env2_distribution_vector[len(df.index)-partition:] = list(range(partition, 0, -1))
                file_name = f"{args.out_dir}/sim_sample_{sim_dict[proportion]}_{i}.csv"
                for col in df.columns[:50]:
                    df[col] = np.random.dirichlet(env1_distribution_vector, 1)[0]
                for col in df.columns[50:]:
                    df[col] = np.random.dirichlet(env2_distribution_vector, 1)[0]
                df.to_csv(file_name)
        else:
            for pos in mean_pos:
                for height in peak_heights:
                    env1_distribution_vector = np.ones(len(df.index))
                    env1_start_pos = pos - height
                    env1_end_pos = pos + height
                    distribution_vector = list(range(1, height+1, 1)) + list(range(height, 0, -1))
                    env1_distribution_vector[env1_start_pos: env1_end_pos] = distribution_vector
                    env2_distribution_vector = np.ones(len(df.index))
                    env2_start_pos = len(df.index) - pos - height
                    env2_end_pos = env2_start_pos + height * 2
                    env2_distribution_vector[env2_start_pos:env2_end_pos] = distribution_vector
                    file_name = f"{args.out_dir}/sim_sample_sim_{similarity_levels[pos]}_spread_{height*2}_{i}.csv"
                    print(file_name)
                for col in df.columns[:50]:
                    df[col] = np.random.dirichlet(env1_distribution_vector, 1)[0]
                for col in df.columns[50:]:
                    df[col] = np.random.dirichlet(env2_distribution_vector, 1)[0]
                df.to_csv(file_name)



if __name__ == '__main__':
    main()
