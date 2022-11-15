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
from src.KEGG_helpers import make_nodes_readable
from networkx.drawing.nx_pydot import graphviz_layout
import matplotlib.pyplot as plt
import networkx as nx
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
y, h, col = max_val*(1 + 0.05), 0.01*max_val*(1 + 0.05), 'k'
ax.plot([x1, x1, x2 - 0.05, x2 - 0.05], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p01.pvalue), ha='center', va='bottom', color=col)
x1, x2 = 0, 2
y, h, col = max_val*(1 + 0.09), 0.01*max_val*(1 + 0.09), 'k'
ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
ax.text((x1+x2)*.5, y+h, "p=" + "{0:.5g}".format(p02.pvalue), ha='center', va='bottom', color=col)
x1, x2 = 1, 2
y, h, col = max_val*(1 + 0.05), 0.01*max_val*(1 + 0.05), 'k'
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
DZ_var_diffabs = np.var(list(twin_pair_to_dab.values()), axis=0)

twin_pair_to_dab = dict()
for twin_pair in MZ_twin_pairs:
    twin_pair_to_dab[twin_pair] = diffabs[basis_to_index[twin_pair[0]], basis_to_index[twin_pair[1]]].todense()

MZ_ave_diffabs = np.mean(list(twin_pair_to_dab.values()), axis=0)
MZ_var_diffabs = np.var(list(twin_pair_to_dab.values()), axis=0)

# There's an interesting outlier in the variances
sns.scatterplot(x=DZ_var_diffabs, y=MZ_var_diffabs)
plt.show()

large_var_DZ = np.argmax(DZ_var_diffabs)
large_var_MZ = np.argmax(MZ_var_diffabs)
print(f"Large var DZ: {diffab_3rd_dim_basis[large_var_DZ]}")
print(f"Large var MZ: {diffab_3rd_dim_basis[large_var_MZ]}")



# Let's see if we can visualize these on a tree
edge_file = os.path.join('data', "kegg_ko_edge_df_br_ko00001.txt_motifs_lengths_n_50_f_10_r_100.txt")
Gdir = LH.import_graph(edge_file, directed=True)
Tint, lint, nodes_in_order, EMDU_index_2_node = LH.weighted_tree_to_EMDU_input(Gdir)
# set all edge lengths to zero
for u, v in Gdir.edges():
    Gdir[u][v]['edge_length'] = 0
    Gdir[u][v]['weight'] = .00000000000000000000001
    Gdir[u][v]['color'] = 'k'


# So it looks like the largest variance ones in DZ are: protein kinases, iron complex outermembrane receptors,
# bacterial motility proteins, and quorum sensing.
# Same with the MZ ones, interestingly

indexer = EMDU.DiffabArrayIndexer(diffabs, nodes_in_order, basis, EMDU_index_2_node)
#for twin_pair in MZ_twin_pairs:
#    print(indexer.get_diffab_for_node(twin_pair[0]+".sig.zip_gather_k_7.csv", twin_pair[1]+".sig.zip_gather_k_7.csv",
#                                      "K00463").todense())


topN = 20
diffab_of_interest = MZ_ave_diffabs
# Get the top 10 diffabs
top_diffabs = np.argsort(np.abs(diffab_of_interest))[-topN:]
# add all the weights to Gdir
for i, diff_val in enumerate(diffab_of_interest):
    u = diffab_3rd_dim_basis[i]
    if u != 'root':
        v = list(Gdir.predecessors(u))[0]
        Gdir[v][u]['weight'] = np.abs(diff_val)
        if diff_val > 0:
            Gdir[v][u]['color'] = 'r'
        else:
            Gdir[v][u]['color'] = 'b'
# then select the nodes corresponding to the top_diffabs
important_vertices = set()
for i in top_diffabs:
    u = diffab_3rd_dim_basis[i]
    if u != 'root':
        v = list(Gdir.predecessors(u))[0]
        important_vertices.add(u)
        important_vertices.add(v)
    else:
        important_vertices.add(u)
# also add all the ancestors of the top diffabs
ancestors_of_importants = set()
for u in important_vertices:
    ancestors_of_importants.update(list(nx.ancestors(Gdir, u)))
important_vertices.update(ancestors_of_importants)
T = Gdir.subgraph(important_vertices)
T = make_nodes_readable(T)
# rename nodes to escape : in the names
T = nx.relabel_nodes(T, {node: node.replace(':', '_') for node in T.nodes()})

node_size_by_degree = [ T.degree(node)*15 for node in T.nodes() ]
new_labels = {}
for u in T.nodes():
    if T.degree(u) <= 5:
        new_labels[u] = u
        # new_labels[u] = ''
    else:
        new_labels[u] = u
widths_orig = np.array([T[u][v]['weight'] for u, v in T.edges()])
widths = widths_orig / np.max(widths_orig) * 30
#widths = np.array([1e7*T[u][v]['weight'] for u, v in T.edges()])
widths += 1
colors = [T[u][v]['color'] for u, v in T.edges()]
pos = graphviz_layout(T, prog="twopi")
plt.figure(figsize=(50, 50))
nx.draw(T, pos, node_size=node_size_by_degree, alpha=0.7, with_labels=True, arrows=False, arrowsize=0, width=widths,
        edge_color=colors, labels=new_labels, font_size=20)
plt.savefig('MZ_ave_diffabs_top_20.png')


