#!/usr/bin/env python
# explore the pairwise distances. This script is intended to be run interactively
import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram

# get the pairwise distances
#pairwise_dists_file = 'real_data/uniform_pw_fu_ko00001.npy'
#pairwise_dists_file = 'real_data/uniform_pw_fu_ko00001_all.npy'
#pairwise_dists_file = 'AAI_pw_fu_ko00001_all.npy'
#pairwise_dists_file = 'uniform_pw_fu_ko00001_all.npy'
pairwise_dists_file = 'motifs_pw_fu_ko00001_all.npy'

#matplotlib.use('MacOSX')
#matplotlib.use('Qt5Agg')
matplotlib.use('Agg')
#meta_data_file = 'real_data/13984_20211224-094712.txt'
meta_data_file = '13984_20211224-094712.txt'
meta_data = pd.read_csv(meta_data_file, sep='\t')

pairwise_dists = np.load(pairwise_dists_file)
# get the basis
basis_file = f"{pairwise_dists_file}.basis.txt"
with open(basis_file, 'r') as f:
    basis = [os.path.basename(line.strip()) for line in f]
# get the sample name from the file names
sample_names = [os.path.basename(name).split('.')[0] for name in basis]
sample_names_with_prefix = ['13984.' + x for x in sample_names]
# get the metadata for the samples
meta_data_of_dist = meta_data[meta_data['sample_name'].isin(sample_names_with_prefix)]
# cluster the pairwise distances
Z = linkage(pairwise_dists, method='ward')
# plot the dendrogram
plt.figure(figsize=(10, 40))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance')
#dendrogram(Z, labels=sample_names_with_prefix)
labels = [meta_data.loc[meta_data['sample_name'] == sample_names_with_prefix[x]]['disease_type'].values[0] for x in
          range(len(sample_names_with_prefix))]
dendrogram(Z, labels=labels, orientation='left', leaf_font_size=7)
# make the labels bigger
plt.rcParams.update({'font.size': 22})
#plt.show()
plt.savefig(f"dendro_{pairwise_dists_file}.png")
