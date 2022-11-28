#!/usr/bin/env bash
set -e
set -u
set -o pipefail

cd /data/shared_data/TwinsStudy/data/DZ/sketches
sourmash index DZ_merged_k_5.sbt.zip DZ_merged.sig.zip -k 5 --protein
sourmash index DZ_merged_k_7.sbt.zip DZ_merged.sig.zip -k 7 --protein
sourmash index DZ_merged_k_11.sbt.zip DZ_merged.sig.zip -k 11 --protein
sourmash index DZ_merged_k_15.sbt.zip DZ_merged.sig.zip -k 15 --protein
cd /data/shared_data/TwinsStudy/data/MZ/
sourmash index MZ_merged_k_5.sbt.zip MZ_merged.sig.zip -k 5 --protein
sourmash index MZ_merged_k_7.sbt.zip MZ_merged.sig.zip -k 7 --protein
sourmash index MZ_merged_k_11.sbt.zip MZ_merged.sig.zip -k 11 --protein
sourmash index MZ_merged_k_15.sbt.zip MZ_merged.sig.zip -k 15 --protein
