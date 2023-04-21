# FunUnifrac
A repository to implement UniFrac, but on functional profiles of metagenomic data.

## Installation
To install the package, run the following command:
```bash
conda env update -f environment.yml
conda activate fununifrac
```

## Preprocessed data support
All required resrouces are shared by Onedrive currently. [Please check here.](https://pennstateoffice365-my.sharepoint.com/:f:/g/personal/akp6031_psu_edu/Ele7H9T9yuRLuAdBY_WDS90BuTCYwk1z6VYqkmdJj9ATnQ?e=ccM4gy)
Detailed processes of generating the files are described at the "Preprocessing" section.

# Script1: compute_edges.py
The original UniFrac uses a phylogenetic tree, allowing one to measure the difference in composition between two 
samples. For FunUniFrac, we use a KEGG Orthology (KO) tree and functional profiles instead of an OTU table as 
input samples. The only obstacle is that the KO tree does not naturally come with branch lengths. To overcome 
this, we use pairwise AAI (average amino acid identity) as a proxy of the pairwise distances between leaf nodes 
and assign the rest of the branch lengths by solving the linear system as demonstrated below.  
<img width="708" alt="image" src="https://user-images.githubusercontent.com/90921267/233392729-db2874b3-f68e-4481-ac62-eb4175c70ef7.png">

### Input
* edge list: The KEGG hierarchy where all the ancestors of the given brite_id. 
* pairwise distances: (KO1, KO2), edge j is on the shortest path from '
                    'KO1 to KO2'
                    
```bash
wget -O "edge_ko00001.txt" "https://pennstateoffice365-my.sharepoint.com/:t:/g/personal/akp6031_psu_edu/ETUrru_StwJFray9WJ8q4lAB9OwzSwF2Zhze_OvXv6HP4A?download=1"
# 1.59GB
wget -O "KOs_sketched_scaled_10_compare_5" "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/ESd-RVfvWbZKpoky2lE6K8EBNktiaw8elnSbxZGZ9mbzHQ?download=1"
wget -O "KOs_sketched_scaled_10_compare_5.labels.txt" "https://pennstateoffice365-my.sharepoint.com/:t:/g/personal/akp6031_psu_edu/EW39Mh7xqsdKsVAIxM5FWVEBnRMkR4bWs5oLvJlbR1y8fw?download=1"
```
### Run
```python
python ./scripts/compute_edges.py -e edge_ko00001.txt -d KOs_sketched_scaled_10_compare_5 -o ../fununifrac_out -b ko00001 -n 50 -f 10 -r 100 --distance
```

### Output
* fununifrac_edge_lengths_{datetime}.csv: The edge lengths added to the original edge list

# Script2(Optional): create_edge_matrix.py
This is a sub-script called inside of compute_edges.py.
The same input is used to generate the intermediate matrix.
### Run
```python
python ./scripts/create_edge_matrix.py -e edge_ko00001.txt -d KOs_sketched_scaled_10_compare_5 -o ../fununifrac_out/large -b ko00001
```
### Output
* fununifrac_A.npz: The matrix describing the least sequare problem for computing edge lengths.
* fununifrac_A_basis.txt: The column information.

# Script3: compute_fununifrac.py
This is the core function of the functional unifrac package. It takes a directory containing functional profiles 
of the sample, in the format of sourmash gather files, as well as 
an edge list file (in this case, for instance, the output of Script1 in the preprocessing steps above) describing the 
underlying the KEGG tree. compute_fununifrac.py computes the pairwise functional UniFrac values among the sourmash gather
files provided by the user.

### Input
* A directory containing functional profiles with the option to specify the file pattern to match to. The default 
file pattern is `*_gather.csv`.
* An edge list file describing the underlying KO tree
* A path specifying where the output directory

### Run
```python
python ./scripts/compute_fununifrac.py -e {edge} -fd . -o {output} --diffab --force -b ko00001 -a median_abund --L2
```
### Output
* Pairwise functional UniFrac distances matrix saved in .npy format
* Basis of the above matrix in .npy format

### Example
We provide a minimum example to demonstrate the usage of this function below.
#### Download KO tree and sample data
```bash
wget -O "edge_ko00001_lengths.txt" "https://pennstateoffice365-my.sharepoint.com/:t:/g/personal/akp6031_psu_edu/EVbB-ieWK7xDqux7Y7u_4lEBpMi-jCyA7oDkq3RhtrueXQ?download=1"

wget -O "f1077_ihmp_IBD_MSM5LLGF_P_k_5_gather.csv" "https://pennstateoffice365-my.sharepoint.com/:x:/g/personal/akp6031_psu_edu/EeSp1SOl-aJPsJOtcDylQc8B_cSWK_pxKiuPumW5CBJevQ?download=1"
wget -O "f2103_ihmp_IBD_HSM7J4Q3_P_k_5_gather.csv" "https://pennstateoffice365-my.sharepoint.com/:x:/g/personal/akp6031_psu_edu/EVlFAOKUEDhCslA5KBk8s2QBOcUdhrlkhDJXsFRap72cKA?download=1"
wget -O "f3158_ihmp_IBD_PSM6XBW1_P_k_5_gather.csv" "https://pennstateoffice365-my.sharepoint.com/:x:/g/personal/akp6031_psu_edu/EWnhFoHwindPp1q5JxmYbS4BQ9p2OArBGGxQwa0XCP9cSg?download=1"
```
#### Run
```bash
python ./scripts/compute_fununifrac.py -e ./edge_ko00001_lengths.txt -d . -o fununifrac_results --diffab
```
#### Output
* fununifrac_results/fununifrac_out_20230420-150954.main.basis.npy
* fununifrac_results/fununifrac_out_20230420-150954.main.npy
* fununifrac_results/ununifrac_out_20230420-150954.diffab.nodes.npy
* fununifrac_results/fununifrac_out_20230420-150954.diffab.npz


# Preprocessing
## 0. KO tree
The KEGG Orthology (KO) hierarchy tree was donwloaded using the [KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html#list). The [KEGG BRITE Database](https://www.genome.jp/kegg/brite.html) is a collection of hierarchical classification system in which [ko00001](https://www.genome.jp/brite/ko00001) is the KO hierarchy. Therefore, we utilized the [get](https://www.kegg.jp/kegg/rest/keggapi.html#get) operation of KEGG API to access the [KO hierarchy](https://rest.kegg.jp/get/br:ko00001). For each KO leaf node in this hierarchy tree, we get its corresponding lineage. For instance, the lineage of KO `K03382` is "09100 Metabolism; 09111 Xenobiotics biodegradation and metabolism; 00791 Atrazine degradation; K03382  atzB; hydroxydechloroatrazine ethylaminohydrolase". With this lineage information, we save the KO hierarhical tree into a file `kegg_ko00001_edges.txt` with an edge list format.  

### Download
```bash
wget -O "kegg_ko00001_edges.txt" "https://pennstateoffice365-my.sharepoint.com/:t:/g/personal/akp6031_psu_edu/EWWFfRZc4o5OltW7YZFH8yYBCQB8HW8IW_zLveO7eAeMPQ?download=1"
```
### Script
```bash
python ./scripts/reproducibility/generate_ko_hierarchy.py --outdir [Path to save the output file]
```

## 1. KO sketches
From the KEGG database, we downloaded all the KO's and used sourmash sketch to produce the sketch for each of the KO's.
These sketches can be used as a reference to efficiently generate the functional profiles for metagenomic samples, which 
will be explained in the Section 3. Sample sketches below.
To use the pre-generated KO sketches directly, simply download using the command below. 
Alternatively, the user can download the raw data and perform sourmash sketch using the parameters of their own choice.
For usage of sourmash, refer to the Appendix section.

### Download sketched KO
```bash
wget -O "KOs_sketched_scaled_10.sig.zip" "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/ESayLPk49QhIhmPExvgNyrIBHhSyk4hzkdg5dPNHMGL5Mg?download=1"
```
### Generate KO sketches from raw KO

Below is command by which we generated the KO sketches above, where $INPUT is the raw KO data downloaded from KEGG.
```bash
sourmash sketch fromfile -p protein,k=5,k=7,k=11,abund,scaled=${scaleFactor} -o ${OUTDIR}/KOs_sketched_scaled_${scaleFactor}.sig.zip $INPUT
```

## 2. KO leaf pairwise distances
[explain]
### Download
```bash
wget -O "KOs_sketched_scaled_10_compare_5" "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/ESsA10lqAXtOrv0Z_mVVedABTBdgHLAXzeE3N8bZ6q9f7Q?download=1"
wget -O "KOs_sketched_scaled_10_compare_5.labels.txt" "https://pennstateoffice365-my.sharepoint.com/:t:/g/personal/akp6031_psu_edu/Ea4iI9gCtl5Nue0bVABSqGQBIR1cv79l3fUk2bHSv6pxMA?download=1"
```
### Script
```bash
```

## 3. sequence sketches
[explain]
### Download
```bash
# sequences
wget -O "f1077_ihmp_IBD_MSM5LLGF_P.fastq" "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/EeFqDGl6lZZIgGFPme4GkBMBK20P_y0qrjqDWDuWnZ6ImA?download=1"
wget -O "f2103_ihmp_IBD_HSM7J4Q3_P.fastq" "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/EaRRHDpkIK5Nid9pFf9BRAsB8Y-KT8lA4Zp1cdYv1eGJ1Q?download=1"
wget -O "f3158_ihmp_IBD_PSM6XBW1_P.fastq" "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/ER7xvoZtBSRKlw7YoYJrfZ8BiLkxdqa-43suUolSUS86Ng?download=1"
# outputs
wget -O "f1077_ihmp_IBD_MSM5LLGF_P.fastq.gz.sig" "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/Ebzry6ZzMrpAocLia1gb8zgBzSThqyvc7YZN_lcJ8Fkofg?download=1"
wget -O "f2103_ihmp_IBD_HSM7J4Q3_P.fastq.sig" "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/EdSy4EMYhJVOmorMSHyOBCMBu5HshZMur5BKG49D7DgSAA?download=1"
wget -O "f3158_ihmp_IBD_PSM6XBW1_P.fastq.sig" "https://pennstateoffice365-my.sharepoint.com/:u:/g/personal/akp6031_psu_edu/EQsLhGlmu09KuY6QDm6ymLMBi9PgIRADjL0ZLQCHSFEmwg?download=1"
```
### Script
```bash
```

## 4. sequence and KOs comparisons
[explain]
### Download
```bash
wget -O "f1077_ihmp_IBD_MSM5LLGF_P_k_5_gather.csv" "https://pennstateoffice365-my.sharepoint.com/:x:/g/personal/akp6031_psu_edu/EYe7vV9JDBhIgY0v0otSz8IBzJ9sp-LXFZQkM2qV7n5d0Q?download=1"
wget -O "f2103_ihmp_IBD_HSM7J4Q3_P_k_5_gather.csv" "https://pennstateoffice365-my.sharepoint.com/:x:/g/personal/akp6031_psu_edu/Ee-QM-tEoQhPsM-2sKOv3y4BdpHZWIJW_PcItrx7JPTQyQ?download=1"
wget -O "f3158_ihmp_IBD_PSM6XBW1_P_k_5_gather.csv" "https://pennstateoffice365-my.sharepoint.com/:x:/g/personal/akp6031_psu_edu/EVN4JF7jkVhGmdHoCKoSTY4Byz3ybRFAsVtBDSv2VM2wEQ?download=1"
```

### Script
```bash
```

## Appendix
### sourmash
Installation directions and usage of sourmash can be found at https://sourmash.readthedocs.io/en/latest/.

# Outdated
## Producing sourmash files
To produce the sourmash gather files of protein sequences mentioned above from DNA reads, one can do the following.

### Installing sourmash
Installation directions of sourmash can be found at https://sourmash.readthedocs.io/en/latest/.

### An example usage
The detailed explanations of the usage can be found at the same site above. To give a brief example, one can run something like
```
nohup sourmash sketch translate -p k=5,k=7,abund,scaled=100 -o sketches/all_sketches_k_5_7_scaled_100.sig.zip --from-file fasta_files.txt &
```
This process translates the DNA sequences into protein sequences and creates sketches for these protein sequences.

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
./scripts/./edges_preprocess_ko.py -e <edge list: kegg_ko_edge_df_br:ko00001.txt> -d <distance matrix: 
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
scripts/./edges_computation.py -e kegg_ko_edge_df_br\:ko00001.txt -d KOs_sketched_scaled_10_compare_5 -A 
ko00001_KOs_sketched_scaled_10_compare_5_A.npz -b ko00001 -n 100 -f 10 -r 1 --force -o kegg_ko_edge_df_br\:ko00001.txt
```
The output will be an edge list with an additional column specifying the length of that edge. Eg. `kegg_ko_edge_df_br:ko00001.txt`

# Testing
To run the tests, you will need to install the `pytest` package (via conda). Then, you can run the tests via:
```bash
cd test
pytest -v .
```
