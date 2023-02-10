#!/usr/bin/env bash
set -e
set -u
set -o pipefail

cd /data/shared_data/TwinsStudy/data/DZ/
for line in `cat ../../metadata/metadata_DZ.csv | rev | cut -d',' -f1 | rev`; do
	nohup sourmash sketch translate -f -p k=5,k=7,k=11,k=15,abund,scaled=1000 --merge ${line}.sig.zip -o sketches_1000/${line}.sig.zip ${line}* &> nohup_1000.log &
done
cd /data/shared_data/TwinsStudy/data/MZ/
for line in `cat ../../metadata/metadata_MZ.csv | rev | cut -d',' -f1 | rev`; do
        nohup sourmash sketch translate -f -p k=5,k=7,k=11,k=15,abund,scaled=1000 --merge ${line}.sig.zip -o sketches_1000/${line}.sig.zip ${line}* &> nohup_1000.log &
done
