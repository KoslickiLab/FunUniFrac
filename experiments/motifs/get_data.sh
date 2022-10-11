#!/usr/bin/env bash
set -e
set -u
echo "This script won't work until I get sudo back"
exit
thisDir=$(pwd)
edgeListDir=/data/shared_data/KEGG_KO_hierarchy/results/
cp ${edgeListDir}/kegg_ko_edge_df.txt .
cp /data/chunyuma/pathon_detection/data/KEGG/motifs_investigation/final_jindex_results.pkl .
