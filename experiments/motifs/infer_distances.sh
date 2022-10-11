#!/usr/bin/env bash
set -e
set -u
thisDir=$(pwd)
scriptsDir="../../scripts"
pwdist=final_jindex_results_dists_pw_mat.npy
edgeList=kegg_ko_edge_df.txt
BRITE=ko00001
edgeLenOut=${edgeList}_lengths_n_50_f_10_r_100_motifs.txt
AName=${BRITE}_${pwdist}_A.npz
# create the pairwise distance matrix if it hasn't already been created
if [ ! -f ${AName} ] ; then
  python ${scriptsDir}/graph_to_path_matrix.py -e ${edgeList} -d ${pwdist} -o . -b ${BRITE}
fi

# infer the edge lengths
python ${scriptsDir}/create_edge_lengths.py -e ${edgeList} -d ${pwdist} -o ${edgeLenOut} -b ${BRITE} -A ${AName} -n 50 -f 10 -r 100 --force
