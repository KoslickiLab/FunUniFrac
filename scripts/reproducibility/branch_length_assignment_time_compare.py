import os
import time
import pandas as pd


#cur_dir = 'reproducibility'
DIR = './data'
DATA_SIZE = [19, 32, 83, 207, 302, 581]
MY_PATH = '../scripts/reproducibility'

#my method
df = pd.DataFrame(columns=['Size', 'Method', 'Time'])
sizes = []
methods = []
times = []

for d_size in DATA_SIZE:
    print(f"size = {d_size}")
    edge_list = f"data/subtree_size{d_size}_no_edge_lengths.txt"
    dm_file = f"data/subtree_size{d_size}_pw_dist.npy"
    label_file = f"data/subtree_size{d_size}_pw_dist.npy.labels.txt"
    output_file = f"data/subtree_size{d_size}_assigned.txt"
    edge_list_with_brite = f"data/subtree_size{d_size}_ko00001_no_edge_lengths.txt"
    A_matrix_file = f"data/ko00001_subtree_size{d_size}_ko00001_pw_dist.npy_A.npz"
    output_file_with_brite = f"data/subtree_size{d_size}_ko00001_assigned.txt"
    dm_file_with_brite = f"data/subtree_size{d_size}_ko00001_pw_dist.npy"

    files = [edge_list, dm_file, label_file, edge_list_with_brite]
    for file in files:
        if not os.path.isfile(file):
            print(f"file {file} not found")

    sizes.append(d_size)
    methods.append('Deterministic')
    start_time = time.time()
    os.system(f"python real_data_branch_assignment.py -e {edge_list} -dm {dm_file} -l {label_file} -s {output_file}")
    end_time = time.time()
    times.append(end_time - start_time)

    #randomized method
    sizes.append(d_size)
    methods.append('Randomized')
    os.system(f'python ../graph_to_path_matrix.py -e {edge_list_with_brite} -d {dm_file_with_brite} -o data'
              f' -b ko00001')
    start_time = time.time()
    os.system(f'python ../create_edge_lengths.py -d {MY_PATH}/{dm_file_with_brite} -e {MY_PATH}/{edge_list_with_brite} '
              f'-A {MY_PATH}/{A_matrix_file} -b ko00001 -o {MY_PATH}/{output_file_with_brite}')
    end_time = time.time()
    times.append(end_time - start_time)

df['Size'] = sizes
df['Method'] = methods
df['Time'] = times

print(df)
df.to_csv('data/comparison_df.csv')








