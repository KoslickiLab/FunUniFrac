#!/usr/bin/env bash
set -e
set -u
thisDir=$(pwd)
koSketchDir=/data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/
edgeListDir=/data/shared_data/KEGG_KO_hierarchy/results/
cp ${koSketchDir}/KOs_sketched_scaled_10_compare_5 .
cp ${koSketchDir}/KOs_sketched_scaled_10_compare_5_labels.txt .
cp ${edgeListDir}/kegg_ko_edge_df.txt .