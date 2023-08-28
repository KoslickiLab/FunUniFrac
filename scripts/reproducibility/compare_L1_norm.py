import pandas as pd
from Kegg_tree import L1_norm

DIR = 'data'
DATA_SIZE = [19, 32, 83, 207, 302, 581]

df = pd.DataFrame(columns=['Size', 'Method', 'L1'])
L1_col = []
size_col = []
method_col = []


for d_size in DATA_SIZE:
    inferred_file_deterministic = f'{DIR}/subtree_size{d_size}_assigned.txt'
    reference_file = f'{DIR}/subtree_size{d_size}.txt'
    inferred_file_random = f'{DIR}/subtree_size{d_size}_ko00001_assigned.txt'

    inferred_df_deterministic = pd.read_table(inferred_file_deterministic, index_col=0)
    inferred_df_random = pd.read_table(inferred_file_random, index_col=0)
    reference_df = pd.read_table(reference_file, index_col=0)

    reference_lengths = list(reference_df['edge_length'])
    inferred_df_random.drop('ko00001', inplace=True)
    inferred_lengths_deterministic = list(inferred_df_deterministic['edge_length'])
    size_col.append(d_size)
    method_col.append('Deterministic')
    deterministic_L1 = L1_norm(inferred_lengths_deterministic, reference_lengths)
    L1_col.append(deterministic_L1)
    size_col.append(d_size)
    method_col.append('Random')
    inferred_lengths_random = list(inferred_df_random['edge_length'])
    random_L1 = L1_norm(inferred_lengths_random, reference_lengths)
    L1_col.append(random_L1)

df['Size'] = size_col
df['Method'] = method_col
df['L1'] = L1_col

df.to_csv('data/L1_comparison_df.csv')
