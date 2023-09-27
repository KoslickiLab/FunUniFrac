import os
import time
import pandas as pd


#cur_dir = 'reproducibility'
DIR = './data'
DATA_SIZE = [50, 100, 500, 1000, 2000, 3000]
R = [2, 3, 4, 10, 12, 15]
MY_PATH = '../scripts/reproducibility'

for r in R:
    data_dict = {'Size': [], 'Method': [], 'Time': []}
    for d_size in DATA_SIZE:
        edge_list = f"data/size_{d_size}_{r}ary_tree.txt"
        edge_list_no_edge_lengths = f"data/size_{d_size}_{r}ary_tree_no_edge_lengths.txt"
        dm_file = f"data/size_{d_size}_{r}ary_tree_pw_dist.npy"
        label_file = f"data/size_{d_size}_{r}ary_tree_pw_dist.npy.labels.txt"
        ko_edge_list = f"data/ko00001_size_{d_size}_{r}ary_tree.txt"
        ko_edge_list_no_edge_lengths = f"data/ko00001_size_{d_size}_{r}ary_tree_no_edge_lengths.txt"
        ko_dm_file = f"data/ko00001_size_{d_size}_{r}ary_tree_pw_dist.npy"
        determ_output_file = f"data/size_{d_size}_{r}ary_tree_assigned.txt"
        A_matrix_file = f"data/ko00001_ko00001_size_{d_size}_{r}ary_tree_pw_dist.npy_A.npz"
        rand_output_file = f"data/ko00001_size_{d_size}_{r}ary_tree_assigned.txt"

        print(edge_list)
        #create needed files
        #add brite, remove edge lengths, create dm files
        df = pd.read_table(edge_list)
        df_no_lengths = df.drop('edge_length', axis=1)
        df_no_lengths.to_csv(edge_list_no_edge_lengths, sep='\t', index=False)
        df.loc[len(df.index)] = ['ko00001', '0', 0] #this df is altered.
        df.to_csv(ko_edge_list, sep='\t', index=False)
        ko_df_no_lengths = df.drop('edge_length', axis=1)
        ko_df_no_lengths.to_csv(ko_edge_list_no_edge_lengths, sep='\t', index=False)
        #generate pairwise distance files
        if not os.path.isfile(dm_file):
            os.system(f"python make_pw_dist_matrix.py -f {edge_list} -p {dm_file.split('.')[0]}")
        if not os.path.isfile(ko_dm_file):
            os.system(f"python make_pw_dist_matrix.py -f {ko_edge_list} -p {ko_dm_file.split('.')[0]}")

        #my method
        data_dict["Size"].append(d_size)
        data_dict["Method"].append('Deterministic')
        start_time = time.time()
        os.system(f"python real_data_branch_assignment.py -e {edge_list_no_edge_lengths} "
                  f"-dm {dm_file} -l {label_file} -s {determ_output_file}")
        end_time = time.time()
        data_dict["Time"].append(end_time - start_time)

        #prof's method
        data_dict["Size"].append(d_size)
        data_dict["Method"].append('Randomized')
        if not os.path.isfile(f"{MY_PATH}/{A_matrix_file}"):
            os.system(f'python ../graph_to_path_matrix.py -e {ko_edge_list_no_edge_lengths} -d {ko_dm_file} -o data'
                  f' -b ko00001')
        start_time = time.time()
        os.system(f'python ../create_edge_lengths.py -d {MY_PATH}/{ko_dm_file} -e {MY_PATH}/{ko_edge_list_no_edge_lengths} '
                  f'-A {MY_PATH}/{A_matrix_file} -b ko00001 -o {MY_PATH}/{rand_output_file} --force --distance')
        end_time = time.time()
        data_dict["Time"].append(end_time - start_time)


    df = pd.DataFrame(data_dict)
    print(df)
    output_file_name = f"data/compare_time_{r}ary_trees.csv"
    df.to_csv(output_file_name)











