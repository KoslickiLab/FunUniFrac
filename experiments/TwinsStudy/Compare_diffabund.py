# This will compare the diffabund vectors of the most promising FUF plots
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


print(f"DZ_mean_paired: {np.mean(unifrac_of_DZ_twin_pairs)}")
print(f"MZ_mean_paired: {np.mean(unifrac_of_MZ_twin_pairs)}")
print(f"unrelated_mean_paired: {np.mean(unifrac_of_unrelated_pairs)}")

from scipy.stats import mannwhitneyu
p01 = mannwhitneyu(unifrac_of_DZ_twin_pairs, unifrac_of_MZ_twin_pairs)
p02 = mannwhitneyu(unifrac_of_DZ_twin_pairs, unifrac_of_unrelated_pairs)
p12 = mannwhitneyu(unifrac_of_MZ_twin_pairs, unifrac_of_unrelated_pairs)
print(p01.pvalue)
print(p02.pvalue)
print(p12.pvalue)

from scipy.stats import ranksums
p01 = ranksums(unifrac_of_DZ_twin_pairs, unifrac_of_MZ_twin_pairs)
p02 = ranksums(unifrac_of_DZ_twin_pairs, unifrac_of_unrelated_pairs)
p12 = ranksums(unifrac_of_MZ_twin_pairs, unifrac_of_unrelated_pairs)
print(p01.pvalue)
print(p02.pvalue)
print(p12.pvalue)

# add these p-values to the plot
sns.set_style("whitegrid")
sns.set_context("paper")
fig, ax = plt.subplots()
ax = sns.boxplot(data=[unifrac_of_DZ_twin_pairs, unifrac_of_MZ_twin_pairs, unifrac_of_unrelated_pairs])
ax.set_xticklabels(['DZ', 'MZ', 'UN'])
ax.set_ylabel("Weighted FunUniFrac Distance")
ax.set_xlabel("Zygosity")
# add a title
if approach == "unweighted":
    ax.set_title(f"Unweighted UniFrac, Tree: {method}, kSize: {kSize}, scale: {scale}, {metric}")
else:
    ax.set_title(f"Weighted UniFrac, Tree: {method}, kSize: {kSize}, scale: {scale}, {metric}")
max_val = max([max(unifrac_of_DZ_twin_pairs), max(unifrac_of_MZ_twin_pairs), max(unifrac_of_unrelated_pairs)])
x1, x2 = 0, 1
y, h, col = max_val*(1 + 0.03), 0.01*max_val*(1 + 0.03), 'k'
ax.plot([x1, x1, x2 - 0.05, x2 - 0.05], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p01.pvalue), ha='center', va='bottom', color=col)
x1, x2 = 0, 2
y, h, col = max_val*(1 + 0.04), 0.01*max_val*(1 + 0.04), 'k'
ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p02.pvalue), ha='center', va='bottom', color=col)
x1, x2 = 1, 2
y, h, col = max_val*(1 + 0.03), 0.01*max_val*(1 + 0.03), 'k'
ax.plot([x1+.05, x1+.05, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p12.pvalue), ha='center', va='bottom', color=col)
plt.show()


# Get the differential abundance vectors for each twin pair
# remove the npy file extension
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

twin_pair_to_dab = dict()
for twin_pair in MZ_twin_pairs:
    twin_pair_to_dab[twin_pair] = diffabs[basis_to_index[twin_pair[0]], basis_to_index[twin_pair[1]]].todense()

MZ_ave_diffabs = np.mean(list(twin_pair_to_dab.values()), axis=0)

# FIXME: Something is off, the diffabs aren't summing to the unifrac value
print(f"Sum of abs diffab: {np.sum(np.abs(twin_pair_to_dab[twin_pair]))}")
print(f"FUF: {pw_unifrac_df.loc[twin_pair[0], twin_pair[1]]}")


edge_file = os.path.join('data', "kegg_ko_edge_df_br_ko00001.txt_motifs_lengths_n_50_f_10_r_100.txt")
Gdir = LH.import_graph(edge_file, directed=True)
Tint, lint, nodes_in_order, EMDU_index_2_node = LH.weighted_tree_to_EMDU_input(Gdir)
diffab_indexer = EMDU.DiffabArrayIndexer(diffabs, nodes_in_order, basis, EMDU_index_2_node)

val = 0
for i in range(len(basis)):
    new_val = np.sum(np.abs(diffabs[0, i, :]))
    if new_val > val:
        print(f"{i}: {new_val}")
        val = new_val

