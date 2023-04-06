#!/usr/bin/env bash
set -e
set -u
set -o pipefail
echo "This script was written after I did the following commands separately"
echo "This is the best that I could recall what I did, so you may want to run these commands separately"
exit 1
# download the data
wget "https://qiita.ucsd.edu/public_download/?data=raw&study_id=13984" -O raw_data.zip
unzip raw_data.zip
# remove the zip file
rm raw_data.zip
# get the metadata
wget "https://qiita.ucsd.edu/public_download/?data=sample_information&study_id=13984" -O metadata.zip
unzip metadata.zip
cd metadata/templates
cp *.txt ../../
cd ../..
# get a list of all the files
find ~+ -type f -name "*.fastq.gz" > ../absolute_files.txt
# use sourmash to sketch them
nohup sourmash sketch translate -p k=5,k=7,abund,scaled=100 -o sketches/all_sketches_k_5_7_scaled_100.sig.zip --from-file absolute_files.txt &
