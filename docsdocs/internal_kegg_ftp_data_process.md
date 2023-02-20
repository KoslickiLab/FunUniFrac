## KEGG FTP data extraction note

Shaopeng Liu (sml6467@psu.edu)

Last update: 2023.2.19



#### Data location (CPU server):

```
# raw data
cd /data/shared_data/KEGG_FTP/kegg

# temp outdata
cd /data/shared_data/temp_shaopeng_process_kegg
```



#### File notes:

1. "links/genes_ko.list.gz" and "ko.tar.gz:ko/ko_genes.list" are the same file (col switched), I used the latter one



#### Check gene data

1. there is no aggregated fna data, so need to merge from each single organism file
2. there is an aggregated faa file in `genes/fasta`, which is the ref-db for KOALA. However, I didn't use it because it has no matched fna. I'll merge both fna and faa for consistency consideration. But will compare with it then.
3. ignore the 4 empty organism folders in `kegg/genes/organisms_new`
4. based on readme, there should be a `organisms/org/` folder which contains the information, but didn't find it
5. NOT all genes are associated with a KO
6. Besides genes without KO, the record number for gene and pep can also be different for non-ko genes
7. Some records also have no description (though it could be found from the nuc/pep files), but I keep consistant with the gff file
8. some genes appear in more than 1 KO, will do double count. E.g. `https://www.kegg.jp/entry/psp:PSPPH_2593`





```
target_dir=/data/shared_data/temp_shaopeng_process_kegg
kegg_dir=/data/shared_data/KEGG_FTP/kegg
single_organism_dir=/data/shared_data/KEGG_FTP/kegg/genes/organisms
py_code=/data/shared_data/temp_shaopeng_process_kegg/script/build_single_organism_ko_gene_map.py
virus_dir=/data/shared_data/KEGG_FTP/kegg/genes/viruses
addendum_dir=/data/shared_data/KEGG_FTP/kegg/genes/addendum

# get ko-gene match files only
cd ${target_dir}
tar -xvf ${kegg_dir}/genes/ko.tar.gz -C ${target_dir} ko/ko_genes.list
cut -f 1 ./ko/ko_genes.list | sort | uniq -c | awk '{print $2"\t"$1}' > count_genes_per_ko.txt
mv ko/ko_genes.list .
rmdir ko

# merge fna and faa files from single record
cd ${target_dir}
mkdir -p single_organism_gene_ko_map
cd single_organism_gene_ko_map

# required python modules: Bio, gzip
python ${py_code} -i ${single_organism_dir} -o ${PWD}

# also add virus and addendum genes
mkdir -p temp_virus_addendum
cd temp_virus_addendum
ln -s ${virus_dir} .
ln -s ${addendum_dir} .
python ${py_code} -i ${PWD} -o ${PWD}

# addendum will fail because it has no DNA sequence, need manually check
# manually run the "collection_single_folder" function
# 1. set df_nuc = df_pep.copy(), then change all values to "No_seq" after merge. 
# 2. delete output fna files for addendum

rm addendum viruses
cat running_log.txt >> ../running_log.txt && rm running_log.txt
mv *.csv ..
mv *.fna ..
mv *.faa ..
cd ..
rmdir temp_virus_addendum


# merge files
mkdir -p all_record_csv
mv cleaned_organism_*.csv ./all_record_csv

mkdir -p KO_gene_dna
mv KO_genes_*.fna ./KO_gene_dna

mkdir -p KO_gene_aa
mv KO_genes_*.faa ./KO_gene_aa

mkdir -p no_ko_gene_dna
mv no_ko_genes_*.fna ./no_ko_gene_dna

mkdir -p no_ko_gene_aa
mv no_ko_genes_*.faa ./no_ko_gene_aa

cat ./KO_gene_dna/KO_genes_*.fna > ../merged_KO_genes.fna
cat ./KO_gene_aa/KO_genes_*.faa > ../merged_KO_genes.faa
cd ..

# count records number and compare to the ref file
grep "^>" merged_KO_genes.faa | wc -l
grep "^>" merged_KO_genes.fna | wc -l
wc -l ko_genes.list

# check the mismatches
# grep "^>" merged_KO_genes.faa | cut -d"|" -f 1 | sed 's/>//g'  > temp_all_gene_hits.txt
# cut -f 2 ko_genes.list > temp_expected_hits.txt
### grep too slow, use py set difference
cut -d":" -f 1 linkfile_extra_ko_gene.txt | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2,2rn > dist_count_extra_link_file.txt
```



#### Compare results

| Filename            | Record number     | Expected number from KEGG  | Old number |
| ------------------- | ----------------- | -------------------------- | ---------- |
| merged_KO_genes.fna | 22551531 (unique) | 22573108 (dedup: 22551679) | 14250334   |
| merged_KO_genes.faa | 22555590 (unique) |                            | 13888067   |
| KO number           | 25416             | 25485                      | 14606      |

1. 22551679 - 22529393 < 47198 different links -> indicates that there are mismatches for both, possibly due to non-update of the record link file

2. manually checking of the missings turned out to be mismatches between the link file and the kff file for organism, e.g.
   ```
   # check this distribution file
   /data/shared_data/temp_shaopeng_process_kegg/explore_record_dif/dist_count_extra_link_file.txt
   
   #### several organisms have lots of mismatch
   mtr     11745   # no KO record in organism mtr
   tva     5248 # no KO record in organism tva
   vp      264 # no organism vp
   nvi     77 # no KO record in organism nvi
   
   
   
   
   #### there are many organisms with few mismatchs, e.g.
   # in link file:
   ko:K18967	woc:BA177_15925
   
   # in kff file (NO KO information)
   woc:BA177_15925 CDS     383     1152    complement(3515487..3516638)            ANO52474                hypothetical protein
   ```
   
   
