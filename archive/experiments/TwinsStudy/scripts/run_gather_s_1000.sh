#!/usr/bin/env bash
# Some of these will fail, but that is to be expected (things like samples that they didn't actually provide metadata for, blank lines being passed to parallel, etc.)
#set -e
#set -u
#set -o pipefail
#REF=/data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_1000.sig.zip
cd /data/shared_data/TwinsStudy/data/DZ
tail -n +2 ../../metadata/metadata_DZ.csv | rev | cut -d',' -f1 | rev | parallel -n 1 -P 5 'sourmash gather --protein -k 5 --estimate-ani-ci --threshold-bp 100 -o sketches_1000/gather_5/{}.sig.zip_gather_k_5.csv sketches_1000/{}.sig.zip /data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_1000.sig.zip';
tail -n +2 ../../metadata/metadata_DZ.csv | rev | cut -d',' -f1 | rev | parallel -n 1 -P 5 'sourmash gather --protein -k 7 --estimate-ani-ci --threshold-bp 100 -o sketches_1000/gather_7/{}.sig.zip_gather_k_7.csv sketches_1000/{}.sig.zip /data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_1000.sig.zip';
tail -n +2 ../../metadata/metadata_DZ.csv | rev | cut -d',' -f1 | rev | parallel -n 1 -P 5 'sourmash gather --protein -k 11 --estimate-ani-ci --threshold-bp 100 -o sketches_1000/gather_11/{}.sig.zip_gather_k_11.csv sketches_1000/{}.sig.zip /data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_1000.sig.zip';
tail -n +2 ../../metadata/metadata_DZ.csv | rev | cut -d',' -f1 | rev | parallel -n 1 -P 5 'sourmash gather --protein -k 15 --estimate-ani-ci --threshold-bp 100 -o sketches_1000/gather_15/{}.sig.zip_gather_k_15.csv sketches_1000/{}.sig.zip /data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_1000.sig.zip';

cd /data/shared_data/TwinsStudy/data/MZ
tail -n +2 ../../metadata/metadata_MZ.csv | rev | cut -d',' -f1 | rev | parallel -n 1 -P 5 'sourmash gather --protein -k 5 --estimate-ani-ci --threshold-bp 100 -o sketches_1000/gather_5/{}.sig.zip_gather_k_5.csv sketches_1000/{}.sig.zip /data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_1000.sig.zip';
tail -n +2 ../../metadata/metadata_MZ.csv | rev | cut -d',' -f1 | rev | parallel -n 1 -P 5 'sourmash gather --protein -k 7 --estimate-ani-ci --threshold-bp 100 -o sketches_1000/gather_7/{}.sig.zip_gather_k_7.csv sketches_1000/{}.sig.zip /data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_1000.sig.zip';
tail -n +2 ../../metadata/metadata_MZ.csv | rev | cut -d',' -f1 | rev | parallel -n 1 -P 5 'sourmash gather --protein -k 11 --estimate-ani-ci --threshold-bp 100 -o sketches_1000/gather_11/{}.sig.zip_gather_k_11.csv sketches_1000/{}.sig.zip /data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_1000.sig.zip';
tail -n +2 ../../metadata/metadata_MZ.csv | rev | cut -d',' -f1 | rev | parallel -n 1 -P 5 'sourmash gather --protein -k 15 --estimate-ani-ci --threshold-bp 100 -o sketches_1000/gather_15/{}.sig.zip_gather_k_15.csv sketches_1000/{}.sig.zip /data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_1000.sig.zip';
