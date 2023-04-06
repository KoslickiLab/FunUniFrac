#!/usr/bin/env bash
set -e
set -u
set -o pipefail
cd /data/shared_data/TwinsStudy/data
OutDir=/data/shared_data/TwinsStudy/output/
FULoc=/data/dmk333/Repositories/FunUniFrac/scripts/
for scale in "1000" "10000"; do
for treeType in "AAI" "motifs"; do
	EdgeList=kegg_ko_edge_df_br_ko00001.txt_${treeType}_lengths_n_50_f_10_r_100.txt
	for dirName in "DZ" "MZ" "merged"; do
		for kSize in 5 7 11 15; do
			GatherDir=/data/shared_data/TwinsStudy/data/${dirName}/sketches_${scale}/gather_${kSize}
			OutFile="${dirName}_pw_fu_${treeType}_scale_${scale}_k_${kSize}_f_unique_weighted.np"
			${FULoc}./make_all_pw_fununifrac.py -e ${EdgeList} -d ${GatherDir} -o ${OutDir}/${OutFile} -fp "*.csv" -b ko00001 --diffab -v -f --L2 &
			#${FULoc}./make_all_pw_fununifrac.py -e ${EdgeList} -d ${GatherDir} -o ${OutDir}/${OutFile} -fp "*.csv" -b ko00001 -v -f --L2 &
			OutFile="${dirName}_pw_fu_unweighted_${treeType}_scale_${scale}_k_${kSize}_f_unique_weighted.np"
			${FULoc}./make_all_pw_fununifrac.py -e ${EdgeList} -d ${GatherDir} -o ${OutDir}/${OutFile} -fp "*.csv" -b ko00001 -v -f --unweighted --L2 &
		done
	done
done
done

# then do the uniform
treeType="uniform"
EdgeList=kegg_ko_edge_df_br_ko00001_uniform_lengths.txt
for scale in "1000" "10000"; do
        for dirName in "DZ" "MZ" "merged"; do
                for kSize in 5 7 11 15; do
                        GatherDir=/data/shared_data/TwinsStudy/data/${dirName}/sketches_${scale}/gather_${kSize}
                        OutFile="${dirName}_pw_fu_${treeType}_scale_${scale}_k_${kSize}_f_unique_weighted.np"
                        ${FULoc}./make_all_pw_fununifrac.py -e ${EdgeList} -d ${GatherDir} -o ${OutDir}/${OutFile} -fp "*.csv" -b ko00001 --diffab -v -f --L2 &
			${FULoc}./make_all_pw_fununifrac.py -e ${EdgeList} -d ${GatherDir} -o ${OutDir}/${OutFile} -fp "*.csv" -b ko00001 -v -f --L2 &
			OutFile="${dirName}_pw_fu_unweighted_${treeType}_scale_${scale}_k_${kSize}_f_unique_weighted.np"
			#${FULoc}./make_all_pw_fununifrac.py -e ${EdgeList} -d ${GatherDir} -o ${OutDir}/${OutFile} -fp "*.csv" -b ko00001 --unweighted -v -f --L2 &
                done
        done
done

