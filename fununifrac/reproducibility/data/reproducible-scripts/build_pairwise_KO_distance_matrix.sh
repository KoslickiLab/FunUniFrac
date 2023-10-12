#!/bin/bash

input_file=~/FunUniFrac/fununifrac/reproducibility/data/test_reproducibles/KO_AA_sketch_scaled_10_k_5.sig.zip
output_dir=~/FunUniFrac/fununifrac/reproducibility/data/test_reproducibles

/usr/bin/time -av -o ${output_dir}/runlog_scale_10_k_5_sourmash_compare_distance sourmash compare --ksize 5 --distance-matrix --protein --estimate-ani -o ${output_dir}/KO_scale_10_k_5_pw_dist.npy ${input_file}
