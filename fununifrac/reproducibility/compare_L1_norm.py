import pandas as pd
from Kegg_tree import L1_norm

DIR = './data'
DATA_SIZE = [50, 100, 500, 1000, 2000, 3000]
R = [2, 3, 4, 10, 12, 15]

for r in R:
    data_dict = {'Size': [], 'Method': [], 'L1': []}
    for d_size in DATA_SIZE:
        edge_list = f"data/size_{d_size}_{r}ary_tree.txt"
        ko_edge_list = f"data/ko00001_size_{d_size}_{r}ary_tree.txt"
        determ_output_file = f"data/size_{d_size}_{r}ary_tree_assigned.txt"
        rand_output_file = f"data/ko00001_size_{d_size}_{r}ary_tree_assigned.txt"

        inferred_df_deterministic = pd.read_table(determ_output_file, index_col=0)
        inferred_df_random = pd.read_table(rand_output_file, index_col=0)
        reference_df = pd.read_table(edge_list, index_col=0)

        reference_lengths = list(reference_df['edge_length'])
        inferred_df_random.drop('ko00001', inplace=True)
        inferred_lengths_deterministic = list(inferred_df_deterministic['edge_length'])
        data_dict["Size"].append(d_size)
        data_dict["Method"].append('Deterministic')
        deterministic_L1 = L1_norm(inferred_lengths_deterministic, reference_lengths)
        data_dict['L1'].append(deterministic_L1)
        data_dict["Size"].append(d_size)
        data_dict["Method"].append('Randomized')
        inferred_lengths_random = list(inferred_df_random['edge_length'])
        random_L1 = L1_norm(inferred_lengths_random, reference_lengths)
        data_dict["L1"].append(random_L1)

    df = pd.DataFrame(data_dict)
    output_file_name = f"data/compare_L1_{r}ary_trees.csv"
    print(output_file_name)
    print(df)
    df.to_csv(output_file_name)