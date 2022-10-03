#!/usr/bin/env bash
set -e
set -u
thisDir=$(pwd)
scriptsDir="../scripts"
pwdist=KOs_sketched_scaled_10_compare_5
edgeList=kegg_ko_edge_df.txt
BRITE=ko00001
edgeLenOut='${edgeList}_lengths_n_50_f_10_r_100.txt'
# create the pairwise distance matrix
python ${scriptsDir}/graph_to_path_matrix.py -e ${edgeList} -d ${pwdist} -o . -b ${BRITE}

# infer the edge lengths
python ${scriptsDir}/create_edge_lengths.py -e ${edgeList} -d ${pwdist} -o ${edgeLenOut} -b ${BRITE} -A
${BRITE}_${pwdist}_A.npz -n 50 -f 10 -r 100 --force