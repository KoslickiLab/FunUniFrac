#!/usr/bin/env bash
set -e
set -u
set -o pipefail
#REF=/data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_10.sig.zip
REFDIR=/data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/
cd /data/shared_data/TwinsStudy/data/DZ
nohup sourmash multigather --query sketches/DZ_merged.sig.zip --protein -k 5 --threshold-bp 100 --db ${REFDIR}/KOs_sketched_scaled_10_k_5.sbt.zip --outdir gather_5 &
nohup sourmash multigather --query sketches/DZ_merged.sig.zip --protein -k 7 --threshold-bp 100 --db ${REFDIR}/KOs_sketched_scaled_10_k_7.sbt.zip --outdir gather_7 &
nohup sourmash multigather --query sketches/DZ_merged.sig.zip --protein -k 11 --threshold-bp 100 --db ${REFDIR}/KOs_sketched_scaled_10_k_11.sbt.zip --outdir gather_11 &
nohup sourmash multigather --query sketches/DZ_merged.sig.zip --protein -k 15 --threshold-bp 100 --db ${REFDIR}/KOs_sketched_scaled_10_k_15.sbt.zip --outdir gather_15 &
#nohup sourmash multigather --query sketches/DZ_merged.sig.zip --protein -k 5 --threshold-bp 100 --db ${REF} --outdir gather_5 &
#nohup sourmash multigather --query sketches/DZ_merged.sig.zip --protein -k 7 --threshold-bp 100 --db ${REF} --outdir gather_7 &
#nohup sourmash multigather --query sketches/DZ_merged.sig.zip --protein -k 11 --threshold-bp 100 --db ${REF} --outdir gather_11 &
#nohup sourmash multigather --query sketches/DZ_merged.sig.zip --protein -k 15 --threshold-bp 100 --db ${REF} --outdir gather_15 &


#cat ../metadata/metadata_DZ.csv | cut -d',' -f1 | rev | parallel -P 10 'sourmash gather --protein -k 5 --estimate-ani-ci --threshold-bp 100 -o ${line}.sig.zip_gather_k_5.csv ${line}.sig.zip ${REF}'
#cat ../metadata/metadata_DZ.csv | cut -d',' -f1 | rev | parallel -P 10 'sourmash gather --protein -k 7 --estimate-ani-ci --threshold-bp 100 -o ${line}.sig.zip_gather_k_7.csv ${line}.sig.zip ${REF}'
#cat ../metadata/metadata_DZ.csv | cut -d',' -f1 | rev | parallel -P 10 'sourmash gather --protein -k 11 --estimate-ani-ci --threshold-bp 100 -o ${line}.sig.zip_gather_k_11.csv ${line}.sig.zip ${REF}'
#cat ../metadata/metadata_DZ.csv | cut -d',' -f1 | rev | parallel -P 10 'sourmash gather --protein -k 15 --estimate-ani-ci --threshold-bp 100 -o ${line}.sig.zip_gather_k_15.csv ${line}.sig.zip ${REF}'

cd /data/shared_data/TwinsStudy/data/MZ
#nohup sourmash multigather --query sketches/MZ_merged.sig.zip --protein -k 5 --threshold-bp 100 --db ${REFDIR}/KOs_sketched_scaled_10_k_5.sbt.zip --outdir gather_5 &
#nohup sourmash multigather --query sketches/MZ_merged.sig.zip --protein -k 7 --threshold-bp 100 --db ${REFDIR}/KOs_sketched_scaled_10_k_7.sbt.zip --outdir gather_7 &
#nohup sourmash multigather --query sketches/MZ_merged.sig.zip --protein -k 11 --threshold-bp 100 --db ${REFDIR}/KOs_sketched_scaled_10_k_11.sbt.zip --outdir gather_11 &
#nohup sourmash multigather --query sketches/MZ_merged.sig.zip --protein -k 15 --threshold-bp 100 --db ${REFDIR}/KOs_sketched_scaled_10_k_15.sbt.zip --outdir gather_15 &

#nohup sourmash multigather --query sketches/MZ_merged.sig.zip --protein -k 5 --threshold-bp 100 --db ${REF} --outdir gather_5 &
#nohup sourmash multigather --query sketches/MZ_merged.sig.zip --protein -k 7 --threshold-bp 100 --db ${REF} --outdir gather_7 &
#nohup sourmash multigather --query sketches/MZ_merged.sig.zip --protein -k 11 --threshold-bp 100 --db ${REF} --outdir gather_11 &
#nohup sourmash multigather --query sketches/MZ_merged.sig.zip --protein -k 15 --threshold-bp 100 --db ${REF} --outdir gather_15 &


#cat ../metadata/metadata_MZ.csv | cut -d',' -f1 | rev | parallel -P 10 'sourmash gather --protein -k 5 --estimate-ani-ci --threshold-bp 100 -o ${line}.sig.zip_gather_k_5.csv ${line}.sig.zip ${REF}'
#cat ../metadata/metadata_MZ.csv | cut -d',' -f1 | rev | parallel -P 10 'sourmash gather --protein -k 7 --estimate-ani-ci --threshold-bp 100 -o ${line}.sig.zip_gather_k_7.csv ${line}.sig.zip ${REF}'
#cat ../metadata/metadata_MZ.csv | cut -d',' -f1 | rev | parallel -P 10 'sourmash gather --protein -k 11 --estimate-ani-ci --threshold-bp 100 -o ${line}.sig.zip_gather_k_11.csv ${line}.sig.zip ${REF}'
#cat ../metadata/metadata_MZ.csv | cut -d',' -f1 | rev | parallel -P 10 'sourmash gather --protein -k 15 --estimate-ani-ci --threshold-bp 100 -o ${line}.sig.zip_gather_k_15.csv ${line}.sig.zip ${REF}'
