import src.LP_EMD_helper as LH
import src.EMDU as EMDU
import multiprocessing
from itertools import combinations, repeat
import glob
from src.CONSTANTS import BRITES
from src.LP_EMD_helper import get_descendants
import os
edge_list_file = 'experiments/KOtree/kegg_ko_edge_df_br_ko00001.txt_lengths_n_50_f_10_r_100.txt'
directory = 'experiments/QIITA_study/small_gather_results'
file_pattern = "*.csv"
force = True
abundance_key = 'median_abund'
brite = 'ko00001'
Gdir = LH.import_graph(edge_list_file, directed=True)
# Select the subtree of the KEGG hierarchy rooted at the given BRITE ID
descendants = get_descendants(Gdir, brite)
# add the root node back in (this will only have affect if you are using multiple brites. TODO: not yet implemented)
descendants.add('root')
# select the subgraph from the brite to the leaves
Gdir = Gdir.subgraph(descendants)

# Then create the inputs for EMDUniFrac
Tint, lint, nodes_in_order, EMDU_index_2_node = LH.weighted_tree_to_EMDU_input(Gdir)

# Then import all the P vectors
Ps = dict()
fun_files = glob.glob(os.path.join(directory, file_pattern))
fun_files = sorted(fun_files)
for file in fun_files:
    P = EMDU.functional_profile_to_EMDU_vector(file, EMDU_index_2_node,
                                               abundance_key=abundance_key, normalize=True)
    Ps[file] = P
