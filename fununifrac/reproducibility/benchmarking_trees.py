import sys, os
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(ROOT_DIR)
sys.path.append('../../')
sys.path.append('../../src')
import glob
import pandas as pd
from src.objects.func_tree import FuncTreeEmduInput
import src.factory.make_tree as make_tree
import src.factory.make_emd_input as make_emd_input
from src.algorithms.emd_unifrac import EarthMoverDistanceUniFracSolver
import numpy as np
from sklearn.metrics import silhouette_score


RAND_TREE = 'data/kegg_trees/fununifrac_edge_lengths_kegg_ko00001_randomized_method.csv'
UNIFORM_TREE = 'data/kegg_trees/kegg_ko00001_edge_length_1.txt'
DETERMINISTIC_TREE = 'data/kegg_trees/kegg_ko00001_scaled_10_k_5_assigned_positivity_enforced.txt'
RAND_BR_TREE = 'data/kegg_trees/kegg_ko00001_edge_length_random_1.txt'
ADJUSTED_TREE = 'data/kegg_trees/kegg_ko00001_scaled_10_k_5_assigned_adjusted.txt'

BRITE = 'ko00001'

metadata_file = 'data/simulated_data/simulated_metadata.csv'
trees = {
    RAND_TREE: 'randomized_tree',
    UNIFORM_TREE: 'uniform_tree',
    DETERMINISTIC_TREE: 'deterministic_tree',
    RAND_BR_TREE: 'randomly_assigned_tree',
    ADJUSTED_TREE: 'adjusted_tree'
}
input_dir = 'data/simulated_data'
similarity_levels = ['high', 'medium', 'low']
meta = pd.read_csv(metadata_file)
meta_dict = dict(zip(meta['sample'], meta['env']))
print(meta_dict)

def make_fununifrac_inputs(raw_P, input, normalize=True):
    #convert an array into one that's suitable for use
    EMDU_index_2_node = input.idx_to_node
    node_2_EMDU_index = {v: k for k, v in EMDU_index_2_node.items()}
    if normalize:
        raw_P = raw_P/raw_P.sum()
    P = np.zeros(len(EMDU_index_2_node))
    for ko in raw_P.index:
        if ko not in node_2_EMDU_index:
            print(f"Warning: {ko} not found in EMDU index, skipping.")
        else:
            P_index = node_2_EMDU_index[ko]
            P[P_index] = raw_P[ko]
    return P

def compute_pw_fununifrac(tree_path, dataframe_file):
    #compute pw_fununifrac of 1 file
    solver = EarthMoverDistanceUniFracSolver()
    tree = make_tree.import_graph(tree_path, directed=True)
    input: FuncTreeEmduInput = make_emd_input.tree_to_EMDU_input(tree, BRITE)
    sample_df = pd.read_csv(dataframe_file, index_col='name')
    Ps_pushed = {}
    for col in sample_df.columns:
        P = make_fununifrac_inputs(sample_df[col], input)
        P_pushed = solver.push_up_L1(P, input)
        Ps_pushed[col] = P_pushed
    dists, diffabs_sparse = solver.pairwise_computation(Ps_pushed, sample_df.columns, input, False, False)
    return dists, sample_df.columns

df_dict = {
    'tree': [],
    'score': [],
    'percentage': [],
}

percentage_overlap = {
    'high': 80,
    'medium': 50,
    'low': 0
}

for sim in similarity_levels:
    files = glob.glob(f"{input_dir}/sim_*{sim}*.csv")
    for tree in trees:
        for file in files:
            print(file)
            dist_matrix, sample_ids = compute_pw_fununifrac(tree, file)
            labels = [meta_dict[i] for i in sample_ids]
            sil_score = silhouette_score(dist_matrix, labels, metric='precomputed')
            df_dict['tree'].append(trees[tree])
            df_dict['score'].append(sil_score)
            df_dict['percentage'].append(percentage_overlap[sim])
df = pd.DataFrame.from_dict(df_dict)
print(df)
out_file_name = f"data/simulated_data/df_combined.tsv"
df.to_csv(out_file_name, sep='\t', index=False)




