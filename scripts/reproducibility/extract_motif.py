import gzip
from line_profiler_pycharm import profile
import json
import pandas as pd
from itertools import combinations


data_dict = dict()
files = '/scratch/shared_data/KEGG_FTP/kegg/genes/organisms/*/*ent.gz'
#{KO:[set of motifs]

@profile
def parse_file(file):
    print(file)
    with gzip.open(file, 'rt') as f:
        data = f.read()
    entries = data.split('///')
    print(len(entries))
    for entry in entries:
        lines = entry.split('\n')
        for line in lines:
            if line.startswith('ORTHOLOGY'):
                KO = line.split()[1]
            if line.startswith('MOTIF'):
                list_of_motifs = line.split()[2:]
        if KO in data_dict:
            data_dict[KO] += data_dict[KO].union(list_of_motifs)
        else:
            data_dict[KO] = set(list_of_motifs)


for file in files:
    parse_file()

with open('data/motifs_for_KO.json', 'w') as f:
    json.dump(data_dict, f, indent=4)

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
df.to_numpy('data/pw_dist_by_motifs_all.npy')