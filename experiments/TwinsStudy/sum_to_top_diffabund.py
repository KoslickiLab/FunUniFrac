# This script will take the differential abundance vectors between DZ and MZ, take the top few from the means
# or vars, and then sum the gather results up to this internal node. Then we can perform a t-test on these internal
# nodes to see if they are significantly different between DZ and MZ.
import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sys import platform
import sparse
import os
import src.EMDU as EMDU
import src.LP_EMD_helper as LH
from src.KEGG_helpers import make_nodes_readable
from networkx.drawing.nx_pydot import graphviz_layout
import matplotlib.pyplot as plt
import networkx as nx
from scipy.stats import ttest_ind
# import the metadata
if platform == "win32":
    os.chdir("C:/Users/dmk333/PycharmProjects/FunUniFrac/experiments/TwinsStudy")
else:
    import matplotlib
    matplotlib.use('TkAgg')
    os.chdir("/Users/dmk333/Dropbox/Repositories/FunUniFrac/experiments/TwinsStudy")
metadata = pd.read_csv("metadata/metadata_with_file_prefix_and_sample_name.csv", sep=None, engine='python')



# Promising
metric = "L1"
scale = 1000
method = "motifs"
kSize = 7
approach = "weighted"

gather_dir = f"data/merged/sketches_{scale}/gather_{kSize}"

# another promising one
#metric = "L2"
#scale = 10000
#method = "uniform"
#kSize = 11
#approach = "weighted"
# import the pw unifrac distances
#pairwise_dists_file = "output/merged_pw_fu_motifs_scale_10000_k_11_f_unique_weighted.np.npy"

if approach == "unweighted":
    to_insert = "unweighted_"
else:
    to_insert = ""
if metric == "L1":
    pairwise_dists_file = f"output/merged_pw_fu_{to_insert}{method}_scale_{scale}_k_" \
                          f"{kSize}_f_unique_weighted.np.npy"
else:
    pairwise_dists_file = f"output/merged_pw_fu_{to_insert}{method}_scale_{scale}_k_" \
                          f"{kSize}_f_unique_weighted_L2.np.npy"
pw_unifrac = np.load(pairwise_dists_file)
if metric == "L2":
    pw_unifrac = np.sqrt(pw_unifrac)
# import the basis, handling the case that numpy may auto-add the .npy extension
basis_file = f"{pairwise_dists_file}.basis.txt"
if not os.path.exists(basis_file):  # handle the case when the file extension npy wasn't included
    pw_dists_no_ext = os.path.splitext(pairwise_dists_file)[0]
    basis_file = f"{pw_dists_no_ext}.basis.txt"
with open(basis_file, 'r') as f:
    basis = [os.path.basename(line.strip()) for line in f]
# next, match the sample names with the basis
sample_names = set(metadata['sample_name'])
basis_revised = []
for basis_elem in basis:
    found_flag = False
    for i in range(len(basis_elem)):
        basis_prefix = basis_elem[:-i]
        if basis_prefix in sample_names:
            found_flag = True
            basis_revised.append(basis_prefix)
            break
    if not found_flag:
        raise Exception(f"Could not find a match for {basis_elem} in the metadata")

# remove from the metadata any samples that are not in the basis
metadata = metadata[metadata['sample_name'].isin(basis_revised)]

# convert the pw unifrac distances to a dataframe
pw_unifrac_df = pd.DataFrame(pw_unifrac, index=basis_revised, columns=basis_revised)

sample_to_index = dict(zip(basis_revised, range(len(basis_revised))))

# get the DZ samples and MZ samples
DZ_sample_names = set(metadata[metadata['Zygosity'] == 'DZ']['sample_name'])
MZ_sample_names = set(metadata[metadata['Zygosity'] == 'MZ']['sample_name'])

# now do the same thing, but only for the pairs of samples that are within the same family
# get the DZ twin pairs (horribly inefficient, but gets the job done)
name_to_twin = dict(zip(metadata['sample_name'], metadata['Twin pair number']))

DZ_twin_pairs = []
for DZ_sample_name1 in DZ_sample_names:
    for DZ_sample_name2 in DZ_sample_names:
        if DZ_sample_name1 != DZ_sample_name2 and name_to_twin[DZ_sample_name1] == name_to_twin[DZ_sample_name2]:
            DZ_twin_pairs.append((DZ_sample_name1, DZ_sample_name2))

# And MZ
MZ_twin_pairs = []
for MZ_sample_name1 in MZ_sample_names:
    for MZ_sample_name2 in MZ_sample_names:
        if MZ_sample_name1 != MZ_sample_name2 and name_to_twin[MZ_sample_name1] == name_to_twin[MZ_sample_name2]:
            MZ_twin_pairs.append((MZ_sample_name1, MZ_sample_name2))

unifrac_of_DZ_twin_pairs = []
for twin_pair in DZ_twin_pairs:
    unifrac_of_DZ_twin_pairs.append(pw_unifrac_df.loc[twin_pair[0], twin_pair[1]])

unifrac_of_MZ_twin_pairs = []
for twin_pair in MZ_twin_pairs:
    unifrac_of_MZ_twin_pairs.append(pw_unifrac_df.loc[twin_pair[0], twin_pair[1]])

unifrac_of_unrelated_pairs = []
for sample_name in basis_revised:
    for other_sample_name in basis_revised:
            if sample_name != other_sample_name and name_to_twin[sample_name] != name_to_twin[other_sample_name]:
                unifrac_of_unrelated_pairs.append(pw_unifrac_df.loc[sample_name, other_sample_name])

pairwise_dists_file_no_npy = os.path.splitext(pairwise_dists_file)[0]
diffabs_file = pairwise_dists_file_no_npy + ".diffab.npz"
diffabs_basis_file = diffabs_file + ".nodes.txt"
diffabs = sparse.load_npz(diffabs_file)

basis_to_index = {file[:-23]: i for i, file in enumerate(basis)}  # remove the .sig.zip_gather_k_7.csv extension

diffab_3rd_dim_basis = []
with open(diffabs_basis_file) as f:
    for line in f.readlines():
        diffab_3rd_dim_basis.append(line.strip())

twin_pair_to_dab = dict()
for twin_pair in DZ_twin_pairs:
    twin_pair_to_dab[twin_pair] = diffabs[basis_to_index[twin_pair[0]], basis_to_index[twin_pair[1]]].todense()

DZ_ave_diffabs = np.mean(list(twin_pair_to_dab.values()), axis=0)
DZ_var_diffabs = np.var(list(twin_pair_to_dab.values()), axis=0)

twin_pair_to_dab = dict()
for twin_pair in MZ_twin_pairs:
    twin_pair_to_dab[twin_pair] = diffabs[basis_to_index[twin_pair[0]], basis_to_index[twin_pair[1]]].todense()

MZ_ave_diffabs = np.mean(list(twin_pair_to_dab.values()), axis=0)
MZ_var_diffabs = np.var(list(twin_pair_to_dab.values()), axis=0)


large_ave_DZ = np.argmax(DZ_ave_diffabs)
large_ave_MZ = np.argmax(MZ_ave_diffabs)
print(f"Large mean DZ: {diffab_3rd_dim_basis[large_ave_DZ]}")
print(f"Large mean MZ: {diffab_3rd_dim_basis[large_ave_MZ]}")
#important_node = "ko01001"
important_node = "K03654"

# get the tree so we can keep track of internal node descendants
edge_file = os.path.join('data', "kegg_ko_edge_df_br_ko00001.txt_motifs_lengths_n_50_f_10_r_100.txt")
Gdir = LH.import_graph(edge_file, directed=True)
leaves = LH.get_leaf_descendants_from_node(Gdir, important_node)
leaves.add(important_node)
leaves_with_prefix = [f"ko:{leaf}" for leaf in leaves]

# Now sum up the relative abundances of these leaf nodes for each sample
sample_to_important_node_abundances = {}
for sample_name in basis_revised:
    sample_to_important_node_abundances[sample_name] = 0
    gather_results = pd.read_csv(os.path.join(gather_dir, f"{sample_name}.sig.zip_gather_k_{kSize}.csv"), sep=',',
                                 engine='pyarrow')
    sample_to_important_node_abundances[sample_name] = np.sum(gather_results[gather_results['name'].isin(
        leaves_with_prefix)]['f_unique_weighted'])

# plot a histogram of the relative abundances of the important node
plt.figure()
plt.hist(sample_to_important_node_abundances.values(), bins=100)
plt.xlabel("Relative abundance of important node")
plt.ylabel("Number of samples")
plt.show()

# now seperate out the DZ and MZ samples
important_node_abundances_DZ = []
important_node_abundances_MZ = []
for sample_name in MZ_sample_names:
    important_node_abundances_MZ.append(sample_to_important_node_abundances[sample_name])
for sample_name in DZ_sample_names:
    important_node_abundances_DZ.append(sample_to_important_node_abundances[sample_name])

# make boxplots of the relative abundances of the important node
plt.figure()
sns.set_style("whitegrid")
sns.set_context("paper")
fig, ax = plt.subplots()
ax = sns.boxplot(data=[important_node_abundances_MZ, important_node_abundances_DZ])
ax.set_xticklabels(['MZ', 'DZ'])
ax.set_ylabel(f"Total relative abundance of important node {important_node}")
ax.set_xlabel("Zygosity")
plt.show()

# run a t-test to see if the relative abundances of the important node are significantly different between MZ and DZ
t, p = ttest_ind(important_node_abundances_MZ, important_node_abundances_DZ)

# Ok, not much is showing up with the differential abundance vectors. Let's go ahead and do it on the pushed up vectors,
# each internal node
