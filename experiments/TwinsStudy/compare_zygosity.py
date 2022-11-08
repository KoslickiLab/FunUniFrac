import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sys import platform

# import the metadata
if platform == "win32":
    os.chdir("C:/Users/dmk333/PycharmProjects/FunUniFrac/experiments/TwinsStudy")
else:
    import matplotlib
    matplotlib.use('TkAgg')
    os.chdir("/Users/dmk333/Dropbox/Repositories/FunUniFrac/experiments/TwinsStudy")
metadata = pd.read_csv("metadata/metadata_with_file_prefix_and_sample_name.csv", sep=None, engine='python')
# import the pw unifrac distances
#pairwise_dists_file = "output/merged_pw_fu_motifs_scale_10000_k_11_f_unique_weighted.np.npy"
pairwise_dists_file = "output/merged_pw_fu_AAI_scale_10000_k_15_f_unique_weighted.np.npy"
pw_unifrac = np.load(pairwise_dists_file)
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

# get the pw unifrac for all (DZ,MZ) pairs
DZ_MZ_pw_unifrac = pw_unifrac_df.loc[:,list(MZ_sample_names)].loc[list(DZ_sample_names),:].values.flatten()
DZ_MZ_mean = np.mean(DZ_MZ_pw_unifrac)

# get the pw unifrac for all (DZ,DZ) pairs
# Do it properly
DZ_DZ_pw_unifrac = pw_unifrac_df.loc[:,list(DZ_sample_names)].loc[list(DZ_sample_names),:].values[np.triu_indices(len(
    DZ_sample_names), k=1)]
DZ_DZ_mean = np.mean(DZ_DZ_pw_unifrac)

MZ_MZ_pw_unifrac = pw_unifrac_df.loc[:,list(MZ_sample_names)].loc[list(MZ_sample_names),:].values[np.triu_indices(len(
    MZ_sample_names), k=1)]
MZ_MZ_mean = np.mean(MZ_MZ_pw_unifrac)


print(f"DZ_MZ_mean: {DZ_MZ_mean}")
print(f"DZ_DZ_mean: {DZ_DZ_mean}")
print(f"MZ_MZ_mean: {MZ_MZ_mean}")

# Make a boxplot of each of the sets of samples
#sns.set_style("whitegrid")
#sns.set_context("paper")
#fig, ax = plt.subplots()
# use violin plot
#ax = sns.violinplot(data=[DZ_MZ_pw_unifrac, DZ_DZ_pw_unifrac, MZ_MZ_pw_unifrac], inner="quartile")
#ax.set_xticklabels(['DZ-MZ', 'DZ-DZ', 'MZ-MZ'])
#ax.set_ylabel("Weighted FunUniFrac Distance")
#ax.set_xlabel("Zygosity")
# set a title
#ax.set_title("Weighted FunUniFrac Distance Between all pairs of individuals")
#plt.savefig("output/zygosity_boxplot.png", dpi=300, bbox_inches='tight')
#plt.show()

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

# run a wilcoxon signed rank test (doesn't work since I don't have matched data)
#from scipy.stats import wilcoxon
#print(wilcoxon(unifrac_of_DZ_twin_pairs, unifrac_of_MZ_twin_pairs))
#print(wilcoxon(unifrac_of_DZ_twin_pairs, unifrac_of_unrelated_pairs))
#print(wilcoxon(unifrac_of_MZ_twin_pairs, unifrac_of_unrelated_pairs))

# run a t-test
from scipy.stats import ttest_ind
p01 = ttest_ind(unifrac_of_DZ_twin_pairs, unifrac_of_MZ_twin_pairs)
p02 = ttest_ind(unifrac_of_DZ_twin_pairs, unifrac_of_unrelated_pairs)
p12 = ttest_ind(unifrac_of_MZ_twin_pairs, unifrac_of_unrelated_pairs)
print(p01.pvalue)
print(p02.pvalue)
print(p12.pvalue)

from scipy.stats import mannwhitneyu
p01 = mannwhitneyu(unifrac_of_DZ_twin_pairs, unifrac_of_MZ_twin_pairs)
p02 = mannwhitneyu(unifrac_of_DZ_twin_pairs, unifrac_of_unrelated_pairs)
p12 = mannwhitneyu(unifrac_of_MZ_twin_pairs, unifrac_of_unrelated_pairs)
print(p01.pvalue)
print(p02.pvalue)
print(p12.pvalue)

from scipy.stats import kruskal
p01 = kruskal(unifrac_of_DZ_twin_pairs, unifrac_of_MZ_twin_pairs)
p02 = kruskal(unifrac_of_DZ_twin_pairs, unifrac_of_unrelated_pairs)
p12 = kruskal(unifrac_of_MZ_twin_pairs, unifrac_of_unrelated_pairs)
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
ax = sns.violinplot(data=[unifrac_of_DZ_twin_pairs, unifrac_of_MZ_twin_pairs, unifrac_of_unrelated_pairs], inner="quartile")
ax.set_xticklabels(['DZ', 'MZ', 'UN'])
ax.set_ylabel("Weighted FunUniFrac Distance")
ax.set_xlabel("Zygosity")
x1, x2 = 0, 1
y, h, col = pw_unifrac_df.max().max() + 0.02, 0.01, 'k'
ax.plot([x1, x1, x2 - 0.05, x2 - 0.05], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p01.pvalue), ha='center', va='bottom', color=col)
x1, x2 = 0, 2
y, h, col = pw_unifrac_df.max().max() + 0.1, 0.01, 'k'
ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p02.pvalue), ha='center', va='bottom', color=col)
x1, x2 = 1, 2
y, h, col = pw_unifrac_df.max().max() + 0.02, 0.01, 'k'
ax.plot([x1+.05, x1+.05, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p12.pvalue), ha='center', va='bottom', color=col)
plt.savefig("output/zygosity_boxplot_paired.png", dpi=300, bbox_inches='tight')
plt.show()
