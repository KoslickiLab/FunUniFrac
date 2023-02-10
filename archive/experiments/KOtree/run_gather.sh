#!/usr/bin/env bash
set -e
set -u
while getopts f:o: flag
do
    case "${flag}" in
        f) FILE=${OPTARG};;
        o) OUTDIR=${OPTARG};;
        *) echo "Usage: $0 -f <fasta/q file> -o <outdir>" >&2
           exit 1
           ;;
    esac
done
# get basename of file
BASENAME=$(basename $FILE)
SIGFILE=${OUTDIR}/${BASENAME}.sig.zip
REF=/data/shared_data/KEGG_data/sourmash_sketches/output_KOs/KOs_sketched/KOs_sketched_scaled_10.sig.zip
REFNAME=$(basename $REF)
sourmash sketch translate -f -p scaled=1000,k=5 ${FILE} -o ${SIGFILE}
sourmash gather --protein -k 5 --estimate-ani-ci --threshold-bp 500 ${SIGFILE} ${REF} -o ${OUTDIR}/${BASENAME}_${REFNAME}_gather_k_5.csv
