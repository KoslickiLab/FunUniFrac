#!/usr/bin/env bash
set -e
set -u
set -o pipefail
cd /data/shared_data/TwinsStudy/data/DZ
nohup sourmash sig cat *.rmhost.sig.zip -o DZ_merged.sig.zip &

cd /data/shared_data/TwinsStudy/data/MZ
nohup sourmash sig cat *.rmhost.sig.zip -o MZ_merged.sig.zip &
