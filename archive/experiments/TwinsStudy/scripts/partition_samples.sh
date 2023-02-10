#!/usr/bin/env bash
set -e
set -u
set -o pipefail
cd /data/shared_data/TwinsStudy/data
# iterate over the metadata files
for line in `cat ../metadata/metadata_DZ.csv | rev | cut -d',' -f1 | rev`; do
	fileName=${line}.1.fq.gz
	if [ -f $fileName ]; then
		mv $fileName DZ/
	fi
	fileName=${line}.2.fq.gz
        if [ -f $fileName ]; then
                mv $fileName DZ/
        fi
done
# Then do it for the MZ
for line in `cat ../metadata/metadata_MZ.csv | rev | cut -d',' -f1 | rev`; do
        fileName=${line}.1.fq.gz
        if [ -f $fileName ]; then
                mv $fileName MZ/
        fi
        fileName=${line}.2.fq.gz
        if [ -f $fileName ]; then
                mv $fileName MZ/
        fi
done
