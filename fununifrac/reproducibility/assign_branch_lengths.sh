#!/bin/bash

#assigning branch lengths to using randomized method and deterministic method, using pw distance scaled 10 k 5
#current dir = reproducibility

kegg_tree_no_lengths="data/kegg_trees/kegg_ko00001_no_edge_lengths_filtered.txt"
ko_dm_file="data/pairwise_distances/KOs_sketched_scaled_10_k_5"
label_file="data/pairwise_distances/KOs_sketched_scaled_10_k_5.labels.txt"
out_dir="data/kegg_trees"
deterministic_output="data/kegg_trees/kegg_ko00001_scaled_10_k_5_assigned_positivity_enforced.txt"
A_matrix="data/kegg_trees/fununifrac_A.npz"
rand_output="data/kegg_trees/kegg_ko00001_randomized_method.txt"

#deterministic method
#python real_data_branch_assignment.py -e $kegg_tree_no_lengths -dm $ko_dm_file -l $label_file -s $deterministic_output

#randomized method
#python ../create_edge_matrix.py -e $kegg_tree_no_lengths -d $ko_dm_file -o $out_dir -b ko00001
python ../compute_edges.py -d "../../fununifrac/reproducibility/"$ko_dm_file -e "d../../fununifrac/reproducibility/$kegg_tree_no_lengths" -b ko00001 -o "../../fununifrac/reproducibility/data/$rand_output" --force --distance
