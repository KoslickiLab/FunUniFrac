import gzip
from line_profiler_pycharm import profile
import pandas as pd
from itertools import combinations
import glob
import numpy as np

data_dict = dict()
files = glob.glob('/scratch/shared_data/KEGG_FTP/kegg/genes/organisms/*/*ent.gz')
print(files)
#{KO:[set of motifs]

@profile
def parse_file(file):
    print(file)
    with gzip.open(file, 'rt') as f:
        data = f.read()
    entries = data.split('///')
    for entry in entries:
        KO = ""
        list_of_motifs = []
        lines = entry.split('\n')
        for line in lines:
            if line.startswith('ORTHOLOGY'):
                KO = line.split()[1]
            if line.startswith('MOTIF'):
                list_of_motifs = line.split()[2:]
        if KO and len(list_of_motifs) > 0:
            if KO in data_dict:
                data_dict[KO] = data_dict[KO].union(list_of_motifs)
            else:
                data_dict[KO] = set(list_of_motifs)


for file in files:
    parse_file(file)


#compute pairwise distances using jaccard index
def compute_jaccard(KO1, KO2):
    jaccard = len(data_dict[KO1].intersection(data_dict[KO2]))/len(data_dict[KO1].union(data_dict[KO2]))
    return jaccard

KO_list = list(data_dict.keys())
df = pd.DataFrame(index=KO_list, columns=KO_list)
for KO in KO_list:
    df[KO][KO] = 0.
for (KO1, KO2) in combinations(KO_list, 2):
    df[KO1][KO2] = df[KO2][KO1] = compute_jaccard(KO1, KO2)

df.to_csv('data/pw_dist_by_motifs_all.csv')
np_df = df.to_numpy()
print(np_df)
np.save('data/pw_dist_by_motifs_all.npy', np_df)
with open('data/pw_dist_by_motifs_all_labels.txt', 'w') as f:
    for KO in KO_list:
        f.write(f"{KO}\n")