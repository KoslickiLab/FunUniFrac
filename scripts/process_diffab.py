import src.LP_EMD_helper as LH
import src.EMDU as EMDU
import multiprocessing
from itertools import combinations, repeat
import glob
from src.CONSTANTS import BRITES
from src.LP_EMD_helper import get_descendants
import os
import pandas as pd
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

P = Ps[fun_files[0]]
Q = Ps[fun_files[1]]
Z, diffab = EMDU. EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)

EMDU.plot_diffab(nodes_in_order, diffab, "P_label", "Q_label", plot_zeros=False, thresh=0.0001)

# now make a data frame with the differential abundance values contained in it
# Note that the keys of the diffab are in the nodes_in_order order. So the first entry of a
# key tuple is the lower node. i.e. the one I'm using to identify the edge
diffab_node_names = dict()
for (x,y), val in diffab.items():
    diffab_node_names[EMDU_index_2_node[x]] = val
diffab_node_names_df = pd.DataFrame.from_dict(diffab_node_names, orient='index', columns=['diffab'])
# TODO: I will want to change this to be a table with two columns: the names of the P and Q and then the abs(val) and
#  0 (or reversed if negative)
pd.DataFrame({"a":{"b":1,"c":2},"d":{"e":1,"f":3}}).fillna(0)