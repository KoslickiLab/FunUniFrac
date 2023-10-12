#!/bin/bash

input_file=/home/grads/wjw5274/FunUniFrac/fununifrac/reproducibility/data/test_reproducibles/KO_faa_for_sketch.csv
out_dir=~/FunUniFrac/fununifrac/reproducibility/data/test_reproducibles

/usr/bin/time -av -o runlog_scale_10_k_5 sourmash sketch fromfile -p protein,k=5,abund,scaled=10 -o ${out_dir}/KO_AA_sketch_scaled_10_k_5.sig.zip ${input_file}
