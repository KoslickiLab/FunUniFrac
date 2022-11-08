# This script will mirror the analysis in compare_zygosity.py, but use Bray Curtis and Spearman instead of
# Funspearman
import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sys import platform
from scipy.spatial.distance import braycurtis
import glob
from scipy.stats import spearmanr
# Get all the metadata
# import the metadata
if platform == "win32":
    os.chdir("C:/Users/dmk333/PycharmProjects/FunUniFrac/experiments/TwinsStudy")
else:
    import matplotlib
    matplotlib.use('TkAgg')
    os.chdir("/Users/dmk333/Dropbox/Repositories/FunUniFrac/experiments/TwinsStudy")
metadata = pd.read_csv("metadata/metadata_with_file_prefix_and_sample_name.csv", sep=None, engine='python')
sample_names = [x for x in metadata['sample_name'] if x]
gather_dir = "data/merged/sketches_10000/gather_7"
# pull off the relative abundances for each of the samples
# first, get all the file names, using glob
gather_files = glob.glob(f"{gather_dir}/*gather_k_7.csv")
name_to_abund = dict()
name_ind = 9
rel_abund_ind = 4
for gather_file in gather_files:
    file_prefix = ".".join(os.path.basename(gather_file).split(".")[0:5])
    name_to_abund[file_prefix] = dict()
    with open(gather_file) as f:
        first_line = f.readline()
        for line in f.readlines():
            line = line.strip().split(',')
            KO_name = line[name_ind]
            rel_abund = float(line[rel_abund_ind])
            name_to_abund[file_prefix][KO_name] = rel_abund

# get a list of all the KO names
KO_names = set()
for file_prefix in name_to_abund.keys():
    KO_names.update(name_to_abund[file_prefix].keys())
KO_names = list(KO_names)
name_to_abund_vector = dict()
for file_prefix in name_to_abund.keys():
    name_to_abund_vector[file_prefix] = np.array([name_to_abund[file_prefix].get(KO_name,0) for KO_name in KO_names])

# now normalize all the vectors to sum to 1
#NOTE: don't do this since I don't for the FunUnifrac approach either
#for file_prefix in name_to_abund_vector.keys():
#    name_to_abund_vector[file_prefix] = name_to_abund_vector[file_prefix]/sum(name_to_abund_vector[file_prefix])

# now partition the data into the three zygosity groups
DZ_sample_names = set(metadata[metadata['Zygosity'] == 'DZ']['sample_name'])
MZ_sample_names = set(metadata[metadata['Zygosity'] == 'MZ']['sample_name'])
name_to_twin = dict(zip(metadata['sample_name'], metadata['Twin pair number']))

DZ_twin_pairs = []
for DZ_sample_name1 in DZ_sample_names:
    if DZ_sample_name1 in name_to_abund_vector:
        for DZ_sample_name2 in DZ_sample_names:
            if DZ_sample_name2 in name_to_abund_vector:
                if DZ_sample_name1 != DZ_sample_name2 and name_to_twin[DZ_sample_name1] == name_to_twin[DZ_sample_name2]:
                    DZ_twin_pairs.append((DZ_sample_name1, DZ_sample_name2))

# And MZ
MZ_twin_pairs = []
for MZ_sample_name1 in MZ_sample_names:
    if MZ_sample_name1 in name_to_abund_vector:
        for MZ_sample_name2 in MZ_sample_names:
            if MZ_sample_name2 in name_to_abund_vector:
                if MZ_sample_name1 != MZ_sample_name2 and name_to_twin[MZ_sample_name1] == name_to_twin[MZ_sample_name2]:
                    MZ_twin_pairs.append((MZ_sample_name1, MZ_sample_name2))

spearman_of_DZ_twin_pairs = []
for twin_pair in DZ_twin_pairs:
    r = spearmanr(name_to_abund_vector[twin_pair[0]], name_to_abund_vector[twin_pair[1]])[0]
    spearman_of_DZ_twin_pairs.append(r)

spearman_of_MZ_twin_pairs = []
for twin_pair in MZ_twin_pairs:
    r = spearmanr(name_to_abund_vector[twin_pair[0]], name_to_abund_vector[twin_pair[1]])[0]
    spearman_of_MZ_twin_pairs.append(r)

spearman_of_unrelated_pairs = []
for sample_name in name_to_abund_vector.keys():
    if sample_name in name_to_abund_vector:
        for other_sample_name in name_to_abund_vector.keys():
            if other_sample_name in name_to_abund_vector:
                if sample_name in name_to_twin and other_sample_name in name_to_twin:
                    if sample_name != other_sample_name and name_to_twin[sample_name] != name_to_twin[other_sample_name]:
                        r = spearmanr(name_to_abund_vector[sample_name], name_to_abund_vector[other_sample_name])[0]
                        spearman_of_unrelated_pairs.append(r)

# then plot the results
print(f"DZ_mean_paired: {np.mean(spearman_of_DZ_twin_pairs)}")
print(f"MZ_mean_paired: {np.mean(spearman_of_MZ_twin_pairs)}")
print(f"unrelated_mean_paired: {np.mean(spearman_of_unrelated_pairs)}")

# add the wilcoxon test
from scipy.stats import ranksums
p01 = ranksums(spearman_of_DZ_twin_pairs, spearman_of_MZ_twin_pairs)
p02 = ranksums(spearman_of_DZ_twin_pairs, spearman_of_unrelated_pairs)
p12 = ranksums(spearman_of_MZ_twin_pairs, spearman_of_unrelated_pairs)
print(p01.pvalue)
print(p02.pvalue)
print(p12.pvalue)

max_val = max(max(spearman_of_DZ_twin_pairs), max(spearman_of_MZ_twin_pairs), max(spearman_of_unrelated_pairs))
# add these p-values to the plot
sns.set_style("whitegrid")
sns.set_context("paper")
fig, ax = plt.subplots()
ax = sns.violinplot(data=[spearman_of_DZ_twin_pairs, spearman_of_MZ_twin_pairs, spearman_of_unrelated_pairs], inner="quartile")
ax.set_xticklabels(['DZ', 'MZ', 'UN'])
ax.set_ylabel("Spearman R")
ax.set_xlabel("Zygosity")
# add a title
ax.set_title("Spearman R of twin pairs")
x1, x2 = 0, 1
y, h, col = max_val + 0.02, 0.01, 'k'
ax.plot([x1, x1, x2 - 0.05, x2 - 0.05], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p01.pvalue), ha='center', va='bottom', color=col)
x1, x2 = 0, 2
y, h, col = max_val + 0.1, 0.01, 'k'
ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p02.pvalue), ha='center', va='bottom', color=col)
x1, x2 = 1, 2
y, h, col = max_val + 0.02, 0.01, 'k'
ax.plot([x1+.05, x1+.05, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p12.pvalue), ha='center', va='bottom', color=col)
plt.savefig("output/zygosity_spearman_boxplot_paired.png", dpi=300, bbox_inches='tight')
plt.show()

# now do the same thing, but for the Bray Curtis distance
def BC_distance(x, y):
    total = np.sum(np.minimum(x, y))
    return 1 - 2*total/(sum(x) + sum(y))

BC_of_DZ_twin_pairs = []
for twin_pair in DZ_twin_pairs:
    r = BC_distance(name_to_abund_vector[twin_pair[0]], name_to_abund_vector[twin_pair[1]])
    BC_of_DZ_twin_pairs.append(r)

BC_of_MZ_twin_pairs = []
for twin_pair in MZ_twin_pairs:
    r = BC_distance(name_to_abund_vector[twin_pair[0]], name_to_abund_vector[twin_pair[1]])
    BC_of_MZ_twin_pairs.append(r)

BC_of_unrelated_pairs = []
for sample_name in name_to_abund_vector.keys():
    if sample_name in name_to_abund_vector:
        for other_sample_name in name_to_abund_vector.keys():
            if other_sample_name in name_to_abund_vector:
                if sample_name in name_to_twin and other_sample_name in name_to_twin:
                    if sample_name != other_sample_name and name_to_twin[sample_name] != name_to_twin[other_sample_name]:
                        r = BC_distance(name_to_abund_vector[sample_name], name_to_abund_vector[other_sample_name])
                        BC_of_unrelated_pairs.append(r)

# then plot the results
print(f"DZ_mean_paired: {np.mean(BC_of_DZ_twin_pairs)}")
print(f"MZ_mean_paired: {np.mean(BC_of_MZ_twin_pairs)}")
print(f"unrelated_mean_paired: {np.mean(BC_of_unrelated_pairs)}")

# add the wilcoxon test
from scipy.stats import ranksums
p01 = ranksums(BC_of_DZ_twin_pairs, BC_of_MZ_twin_pairs)
p02 = ranksums(BC_of_DZ_twin_pairs, BC_of_unrelated_pairs)
p12 = ranksums(BC_of_MZ_twin_pairs, BC_of_unrelated_pairs)
print(p01.pvalue)
print(p02.pvalue)
print(p12.pvalue)



max_val = max(max(BC_of_DZ_twin_pairs), max(BC_of_MZ_twin_pairs), max(BC_of_unrelated_pairs))
# add these p-values to the plot
sns.set_style("whitegrid")
sns.set_context("paper")
fig, ax = plt.subplots()
ax = sns.violinplot(data=[BC_of_DZ_twin_pairs, BC_of_MZ_twin_pairs, BC_of_unrelated_pairs], inner="quartile")
ax.set_xticklabels(['DZ', 'MZ', 'UN'])
ax.set_ylabel("Bray Curtis")
ax.set_xlabel("Zygosity")
# add a title
ax.set_title("Bray Curtis distance between twin pairs")
x1, x2 = 0, 1
y, h, col = max_val + 0.2, 0.01, 'k'
ax.plot([x1, x1, x2 - 0.05, x2 - 0.05], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p01.pvalue), ha='center', va='bottom', color=col)
x1, x2 = 0, 2
y, h, col = max_val + 0.1, 0.01, 'k'
ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p02.pvalue), ha='center', va='bottom', color=col)
x1, x2 = 1, 2
y, h, col = max_val + 0.2, 0.01, 'k'
ax.plot([x1+.05, x1+.05, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p12.pvalue), ha='center', va='bottom', color=col)
plt.savefig("output/zygosity_BC_boxplot_paired.png", dpi=300, bbox_inches='tight')
plt.show()
