import os
import time
import pandas as pd


#cur_dir = 'reproducibility'
DIR = './data'
DATA_SIZE = [19, 32, 83]
MY_PATH = '../scripts/reproducibility/data'

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
    edge_list_with_brite = f"{MY_PATH}/subtree_size{d_size}_ko00001_no_edge_lengths.txt"
    A_matrix_file = f"{MY_PATH}/ko00001_subtree_size{d_size}_ko00001_pw_dist.npy_A.npz"
    output_file_with_brite = f"{MY_PATH}/subtree_size{d_size}_ko00001_assigned.txt"
    dm_file_with_brite = f"{MY_PATH}/subtree_size{d_size}_pw_dist.npy"

    files = [edge_list, dm_file, label_file]
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
    start_time = time.time()
    os.system(f'python ../create_edge_lengths.py -d {dm_file_with_brite} -e {edge_list_with_brite} '
              f'-A {A_matrix_file} -b ko00001 -o {output_file_with_brite}')
    end_time = time.time()
    times.append(end_time - start_time)

df['Size'] = sizes
df['Method'] = methods
df['Time'] = times

print(df)
df.to_csv('data/comparison_df.csv')








