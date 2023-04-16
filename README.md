# FunUniFrac
A repository to implement UniFrac, but on functional profiles of metagenomic data.

## Installation
To install the package, run the following command:
```bash
# conda create -n fununifrac -c bioconda -c conda-forge -c conda --file requirements.txt
conda env update -f environment.yml
conda activate fununifrac
```

## Running FunUniFrac
### Required inputs
This is the core function of the functional unifrac package. It is performed by running the `make_all_pw_fununifrac` 
command. It requires three inputs: 1. a file representing the underlying functional hierarchy in the form of an edge list.
For now, the accepted data is a `KEGG` hierarchy file. An example is provided in the `data` directory, with the name
`kegg_ko_edge_df_br_ko00001.txt`. 2. a directory containing sourmash files. An additional argument `-fp` can be added to
indicate the file pattern to match in this directory. An example is `*_gather.csv`, which is also the default option.
The `KEGG` hierarchy consists of many trees each rooted at a brite (for more information on brites, refer to
https://www.genome.jp/kegg/brite.html). 

### Example 
```bash
python fununifrac.py -d sourmash_results -e kegg_ko_edge_df_br_ko00001.txt_AAI_lengths_n_50_f_10_r_100.txt -fp '*.csv' -b ko00001 -o results.npy
```

## A brief explanation of data structures underlying FunUniFrac

### Setting edge lengths

#### Obtaining the underlying tree/DAG structure
The underlying graph structure on which the earth mover's distance
is based on is obtained from the KEGG database. The KEGG database 
[extraction repo](https://github.com/KoslickiLab/KEGG_data_extraction)
is used to create the edge-list of this graph. In particular:
[get_ko_hierarchy.py](https://github.com/KoslickiLab/KEGG_data_extraction/blob/master/python_scripts/get_ko_hierarchy.py)

Currently, to keep things simple, I am dealing with just sub**trees** of the KEGG database. 
As such, you need to pick a parent [BRITE node](https://rest.kegg.jp/list/brite) to start with. Later, we will extract the tree
of all the children of this parent node. For example, you can choose the BRITE: ko00001 which 
will correspond to "KEGG Orthology (KO)". In the following, call the resulting tree the "KO tree".
In the end, you should have a file like `kegg_ko_edge_df_br:ko00001.txt`.

#### Creating pairwise distances
To infer edge lengths on this tree, you will need a matrix of all pairwise distances between the KEGG leaves
of the KO tree. This can be accomplished via:
1. Sketch all the KOs (note, this will take forever, so just use the pre-computed ones on the server:
```bash
#!/usr/bin/env bash
set -e
set -u
CSVINPUT=/data/shared_data/KEGG_data/PartitionByKO/data/KOs_for_sketch_fromfile.txt
OUTDIR=/data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched
#for scaleFactor in 10 100 500 1000
for scaleFactor in 1
do
        /usr/bin/time sourmash sketch fromfile -p protein,k=5,k=7,k=11,abund,scaled=${scaleFactor} -o ${OUTDIR}/KOs_sketched_scaled_${scaleFactor}.sig.zip $CSVINPUT
done
```
2. Compute the pairwise distances between all the KOs
```bash
#!/usr/bin/env bash
set -e
set -u
DATADIR=/data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched
for kSize in 5 7 11 15
do
        sourmash compare --protein --no-dna -o ${DATADIR}/compare_k_${kSize} --ani -k ${kSize} -p 50 ${DATADIR}/KOs_sketched_scaled_10.sig.zip
done
```
At the end of this, you will have a pairwise distance matrix for all the KOs, along with a text file telling you which
KO corresponds to which row/column in the matrix. Eg. of such files are: `KOs_sketched_scaled_10_compare_5`
and `KOs_sketched_scaled_10_compare_5.labels.txt`.

#### Creating the matrix of all edges traversed between all pairs of KOs
In order to infer the edge lengths, we need to know which edges are traversed between all pairs of KOs.
For example, if we have a binary tree with 2 leaves:
```
    1
L1 / \ L2
  2   3
``` 
and edge lengths `L1` and `L2` respectively, then the distance between the two leaves is `L1 + L2`, resulting in a 
matrix equation of the form:
```
    [1, 1] [L1;L2] = [d(a,b)] 
``` 
where `d(a,b)` is the distance between the two leaves `a` and `b` derived in the previous section.
The coefficient matrix of all edges traversed between all pairs of KOs is obtained by running:
```bash
./scripts/./graph_to_path_matrix.py -e <edge list: kegg_ko_edge_df_br:ko00001.txt> -d <distance matrix: 
KOs_sketched_scaled_10_compare_5> -o <output directory> -b <BRITE: ko00001> 
```
The resulting matrix will have a name such as: `ko00001_KOs_sketched_scaled_10_compare_5_A.npz` and will be a 
massive (233M rows, ~20K columns) sparse matrix called `A` below.

#### Inferring edge lengths
We will employ a regularized, randomized nonnegative least squares approach to infer the edge lengths.
In practice, we've observed that the matrix `A` is rank deficient, yet over-determined (more rows than columns).
Thus there is no unique solution to `Ax=y` where `y` is the vector of pairwise distances between all KOs. We can, 
however, find a unique "shortest length" solution by adding a regularization term to the objective function:
```
min ||Ax - y||^2 + lambda * ||x||^2
s.t. x>=0
```
where `lambda` is a regularization parameter. This can be cast as a NNLS problem via:
```
min ||A'x - y'||^2
s.t. x>=0
where
A' = [A; sqrt(lambda)*ones(1, n)]
y' = [y; 0]
```
We can then solve this problem via a randomized approach, selecting enough rows of `A'` to ensure that the rank is 
large enough. After a bunch of iterates, the average is then taken. An example of running the code is:
```bash
scripts/./create_edge_lengths.py -e kegg_ko_edge_df_br\:ko00001.txt -d KOs_sketched_scaled_10_compare_5 -A 
ko00001_KOs_sketched_scaled_10_compare_5_A.npz -b ko00001 -n 100 -f 10 -r 1 --force -o kegg_ko_edge_df_br\:ko00001.txt
```
The output will be an edge list with an additional column specifying the length of that edge. Eg. `kegg_ko_edge_df_br:ko00001.txt`

# Testing
To run the tests, you will need to install the `pytest` package (via conda). Then, you can run the tests via:
```bash
cd test
pytest -v .
```
