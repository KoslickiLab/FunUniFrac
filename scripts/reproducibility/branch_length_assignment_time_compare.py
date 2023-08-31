import os
import time
import pandas as pd


#cur_dir = 'reproducibility'
DIR = './data'
#DATA_SIZE = [50, 100, 500, 1000, 2000, 3000]
#DATA_SIZE = ['50', '100', '500']
DATA_SIZE = ['1000', '2000', '3000']
R = [2, 3, 4]
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
        os.system(f"python make_pw_dist_matrix.py -f {edge_list} -p {dm_file.split('.')[0]}")
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
        os.system(f'python ../graph_to_path_matrix.py -e {ko_edge_list_no_edge_lengths} -d {ko_dm_file} -o data'
                  f' -b ko00001')
        start_time = time.time()
        os.system(f'python ../create_edge_lengths.py -d {MY_PATH}/{ko_dm_file} -e {MY_PATH}/{ko_edge_list_no_edge_lengths} '
                  f'-A {MY_PATH}/{A_matrix_file} -b ko00001 -o {MY_PATH}/{rand_output_file} --force')
        end_time = time.time()
        data_dict["Time"].append(end_time - start_time)


    df = pd.DataFrame(data_dict)
    print(df)
    output_file_name = f"data/compare_time_{r}ary_trees_big.csv"
    df.to_csv(output_file_name)


# for d_size in DATA_SIZE:
#     print(f"size = {d_size}")
#     edge_list = f"data/subtree_size{d_size}_no_edge_lengths.txt"
#     dm_file = f"data/subtree_size{d_size}_pw_dist.npy"
#     label_file = f"data/subtree_size{d_size}_pw_dist.npy.labels.txt"
#     output_file = f"data/subtree_size{d_size}_assigned.txt"
#     edge_list_with_brite = f"data/subtree_size{d_size}_ko00001_no_edge_lengths.txt"
#     A_matrix_file = f"data/ko00001_subtree_size{d_size}_ko00001_pw_dist.npy_A.npz"
#     output_file_with_brite = f"data/subtree_size{d_size}_ko00001_assigned.txt"
#     dm_file_with_brite = f"data/subtree_size{d_size}_ko00001_pw_dist.npy"
#
#     files = [edge_list, dm_file, label_file, edge_list_with_brite]
#     for file in files:
#         if not os.path.isfile(file):
#             print(f"file {file} not found")
#
#     sizes.append(d_size)
#     methods.append('Deterministic')
#     start_time = time.time()
#     os.system(f"python real_data_branch_assignment.py -e {edge_list} -dm {dm_file} -l {label_file} -s {output_file}")
#     end_time = time.time()
#     times.append(end_time - start_time)
#
#     #randomized method
#     sizes.append(d_size)
#     methods.append('Randomized')
#     os.system(f'python ../graph_to_path_matrix.py -e {edge_list_with_brite} -d {dm_file_with_brite} -o data'
#               f' -b ko00001')
#     start_time = time.time()
#     os.system(f'python ../create_edge_lengths.py -d {MY_PATH}/{dm_file_with_brite} -e {MY_PATH}/{edge_list_with_brite} '
#               f'-A {MY_PATH}/{A_matrix_file} -b ko00001 -o {MY_PATH}/{output_file_with_brite}')
#     end_time = time.time()
#     times.append(end_time - start_time)
#
# df['Size'] = sizes
# df['Method'] = methods
# df['Time'] = times
#
# print(df)
# df.to_csv('data/comparison_df.csv')








