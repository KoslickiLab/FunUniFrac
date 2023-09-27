#!/bin/bash

EDGE_LIST='kegg_ko00001_no_edge_lengths.txt'
PW_LABEL='KOs_sketched_scaled_100_k_11_pw_dist.labels.txt'
FILTERED_FILE='kegg_ko00001_no_edge_lengths_filtered.txt'
touch $FILTERED_FILE	
echo "#parent	child	edge_length" > $FILTERED_FILE
#grep -v K[0-9] $EDGE_LIST >> $FILTERED_FILE
while read -r line
do
	if [[ ! $line == *K[0-9]* ]]; then
		echo "$line" >> $FILTERED_FILE
	fi
	#grep "$line" "$EDGE_LIST" >> "$FILTERED_FILE"
done <$EDGE_LIST

while read -r line
do
	grep "$line" "$EDGE_LIST" >> "$FILTERED_FILE"
done<$PW_LABEL
