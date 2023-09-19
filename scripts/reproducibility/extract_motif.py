import gzip
from line_profiler_pycharm import profile

data_dict = dict()
#{KO:[list of motif]

@profile
def parse_file():
    file = 'data/T07275.ent.gz'
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
            data_dict[KO] += list_of_motifs
        else:
            data_dict[KO] = list_of_motifs


parse_file()
print(len(data_dict))