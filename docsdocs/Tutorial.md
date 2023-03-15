
# Tutorial
A tutorial on the usage of Functional UniFrac.

## Installation
To install the package, run the following command:
```bash
conda create -n fununifrac -c bioconda -c conda-forge -c conda --file requirements.txt
conda activate fununifrac
```

## Usage

### What you will need to produce pairwise functional UniFrac distances
This is the core function of the functional unifrac package. It is performed by running the `make_all_pw_fununifrac` 
command. It requires three inputs: 1. a file representing the underlying functional hierarchy in the form of an edge list.
For now, the accepted data is a `KEGG` hierarchy file. An example is provided in the `data` directory, with the name
`kegg_ko_edge_df_br_ko00001.txt`. 2. a directory containing sourmash files. An additional argument `-fp` can be added to
indicate the file pattern to match in this directory. An example is `*_gather.csv`, which is also the default option.
The `KEGG` hierarchy consists of many trees each rooted at a brite (for more information on brites, refer to
https://www.genome.jp/kegg/brite.html). 

### Example 
```bash
python make_all_pw_fununifrac.py -e kegg_ko_edge_df_br_ko00001.txt_AAI_lengths_n_50_f_10_r_100.txt -fp '*.csv' -b ko00001 -o results.npy
```

#### Producing functional profiles