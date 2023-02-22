## KEGG FTP data extraction note

Shaopeng Liu (sml6467@psu.edu)

Last update: 2023.2.21



Data information:

Location: CPU server `/data/shared_data/cleaned_KEGG_FTP_2022_NOV`

| File                         | Description                                              |
| ---------------------------- | -------------------------------------------------------- |
| **sourmash_sketch/**         | KO AA sketch for scale=100,1000 and k=5,7,11,15,25       |
| SBT/                         | Sequence bloom tree based on sourmash sketches above     |
| **merged_KO_genes.faa**      | gene-level AA sequences with KO information              |
| merged_KO_genes.fna          | gene-level DNA sequences with KO information             |
| script                       | reproducible code                                        |
| single_organism_gene_ko_map/ | organism-level gene metagdata from KEGG FTP database     |
| **PartitionByKO/**           | KO-level AA fasta files (input data for sourmash sketch) |
| count_genes_per_ko.txt       | Count of genes in each KO                                |
| ko_genes.list                | single KO-gene map (or edge) from KEGG FTP databse       |





### Contents:

- [Gene data note](#note)
- [Code: clean and merge all files](#clean)
- [KEGG data summary](#summary)
- [Script code](#code) (note the folder name has been changed)
  1. build_single_organism_ko_gene_map.py
  2. partition_aa_per_ko.sh
  3. sketch_aa_per_ko.sh
  4. build_sbt_lca.sh





#### Data location (CPU server):

```
# raw data
cd /data/shared_data/KEGG_FTP/kegg

# temp outdata
cd /data/shared_data/temp_shaopeng_process_kegg
```



#### File notes:

1. "links/genes_ko.list.gz" and "ko.tar.gz:ko/ko_genes.list" are the same file (col switched), I used the latter one



#### Check gene data <a name="note"></a>

1. there is no aggregated fna data, so need to merge from each single organism file
2. there is an aggregated faa file in `genes/fasta`, which is the ref-db for KOALA. However, I didn't use it because it has no matched fna. I'll merge both fna and faa for consistency consideration. But will compare with it then.
3. ignore the 4 empty organism folders in `kegg/genes/organisms_new`
4. based on readme, there should be a `organisms/org/` folder which contains the information, but didn't find it
5. NOT all genes are associated with a KO
6. Besides genes without KO, the record number for gene and pep can also be different for non-ko genes
7. Some records also have no description (though it could be found from the nuc/pep files), but I keep consistant with the gff file
8. some genes appear in more than 1 KO, will do double count. E.g. `https://www.kegg.jp/entry/psp:PSPPH_2593`



#### Clean all files  <a name="clean"></a>

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



#### KEGG data summary <a name="summary"></a>

| Filename            | Record number     | Expected number from KEGG  | Old number |
| ------------------- | ----------------- | -------------------------- | ---------- |
| merged_KO_genes.fna | 22551531 (unique) | 22573108 (dedup: 22551679) | 14250334   |
| merged_KO_genes.faa | 22555590 (unique) |                            | 13888067   |
| KO number           | 25416             | 25485                      | 14606      |

manually checking of the missings turned out to be mismatches between the link file and the kff file for organism, e.g.
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





#### Script code <a name="py"></a>

1. build_single_organism_ko_gene_map.py
   ```
   import pandas as pd
   import os, argparse, glob
   from Bio import SeqIO
   import gzip
   
   
   def collection_single_folder(f_in):
   	"""
   	Collect gene-ko map in a given folder.
   	:param f_in:
   	:param f_out:
   	:return: a merged df with index: ko_gene_id, columns: [koid, desc, ntseq, aaseq]
   	"""
   	# load files
   	temp_gff = glob.glob(f_in + "/*.kff.gz")[0]
   	temp_nuc = glob.glob(f_in + "/*.nuc.gz")[0]
   	temp_pep = glob.glob(f_in + "/*.pep.gz")[0]
   	# read metadata
   	temp_meta = pd.read_csv(temp_gff, sep='\t', compression='gzip', header=None, index_col='kegg_gene_id',
   							usecols=[0, 9, 10], names=['kegg_gene_id', 'koid', 'desc'])
   	# remove records with NO KO information
   	#temp_meta.dropna(subset=['koid'], inplace=True)
   	
   	# fillna with a special case:
   	temp_meta['koid'] = temp_meta['koid'].fillna("No_ko")
   	temp_meta['desc'] = temp_meta['desc'].fillna("No_desc")
   	
   	# load fna file
   	temp_nuc_dict = {}
   	with gzip.open(temp_nuc, "rt") as handle:
   		for fasta_record in SeqIO.parse(handle, "fasta"):
   			temp_nuc_dict[fasta_record.id] = str(fasta_record.seq.upper())
   	df_nuc = pd.DataFrame(temp_nuc_dict.items(), columns=['kegg_gene_id', 'ntseq'])
   	df_nuc.index = df_nuc.pop('kegg_gene_id')
   	# load faa file
   	temp_pep_dict = {}
   	with gzip.open(temp_pep, "rt") as handle:
   		for fasta_record in SeqIO.parse(handle, "fasta"):
   			temp_pep_dict[fasta_record.id] = str(fasta_record.seq.upper())
   	df_pep = pd.DataFrame(temp_pep_dict.items(), columns=['kegg_gene_id', 'aaseq'])
   	df_pep.index = df_pep.pop('kegg_gene_id')
   	# merge together
   	out_df = pd.concat([temp_meta, df_nuc, df_pep], axis=1)
   	# in case no seq available
   	out_df.fillna("No_seq", inplace=True)
   	
   	return out_df
   		
   def write_merged_df_to_2_fasta(input_df, f_out, out_name):
   	"""
   	Write a merged organism record to 2 fasta file: 1 for KO-gene map, 1 for non-KO genes
   	:param input_df:
   	:param f_out:
   	:param out_name:
   	:return: No return value
   	"""
   	# write this df
   	input_df.to_csv(os.path.join(f_out, "cleaned_organism_"+out_name+".csv"))
   	
   	# write non-ko fasta
   	sub_df = input_df[input_df['koid'] == "No_ko"]
   	# loop through index to write out
   	for item in sub_df.index:
   		temp_row = sub_df.loc[item]
   		temp_record = "|".join([">"+item, temp_row['desc'], temp_row['koid']])
   		temp_fna = temp_row['ntseq']
   		temp_faa = temp_row['aaseq']
   		# write to file
   		with open(os.path.join(f_out, "no_ko_genes_" + out_name + ".fna"), 'a') as f:
   			f.write(temp_record + "\n" + temp_fna + "\n")
   		with open(os.path.join(f_out, "no_ko_genes_" + out_name + ".faa"), 'a') as f:
   			f.write(temp_record + "\n" + temp_faa + "\n")
   			
   	# write KO-gene fasta
   	sub_df = input_df[input_df['koid'] != "No_ko"]
   	# loop through index to write out
   	for item in sub_df.index:
   		temp_row = sub_df.loc[item]
   		temp_record = "|".join([">" + item, temp_row['desc'], temp_row['koid']])
   		temp_fna = temp_row['ntseq']
   		temp_faa = temp_row['aaseq']
   		# write to file
   		with open(os.path.join(f_out, "KO_genes_" + out_name + ".fna"), 'a') as f:
   			f.write(temp_record + "\n" + temp_fna + "\n")
   		with open(os.path.join(f_out, "KO_genes_" + out_name + ".faa"), 'a') as f:
   			f.write(temp_record + "\n" + temp_faa + "\n")
   		
   def local_variable():
   	"""
   	For test purpose
   	:return:
   	"""
   	os.chdir("./clean_kegg_db")
   	input_dir = os.path.abspath("./input_dir")
   	out_dir = os.getcwd()
   	logfile = "running_log.txt"
   
   
   if __name__ == '__main__':
   	parser = argparse.ArgumentParser(description="Build single file gene-ko map",
   									 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
   	parser.add_argument('-i', '--input_dir', help="Path to input folder of single organisms", required=True)
   	parser.add_argument('-o', '--output_dir', help="Path to output folder", required=True)
   	parser.add_argument('-l', '--logfile', help="Log file to write", default="running_log.txt")
   
   	args = parser.parse_args()
   	input_dir = args.input_dir
   	out_dir = args.output_dir
   	logfile = args.logfile
   
   	for record in os.listdir(input_dir):
   		temp_dir = os.path.join(input_dir, record)
   		# add a dir filter (in case of unexpected single file)
   		if os.path.isdir(temp_dir):
   			name = record
   			try:
   				temp_df = collection_single_folder(f_in=temp_dir)
   			except Exception:
   				print("Failed on %s" % name)
   				with open(os.path.join(out_dir, logfile), 'a') as f:
   					f.write("Fail\t%s\n" % name)
   			else:
   				with open(os.path.join(out_dir, logfile), 'a') as f:
   					f.write("Success\t%s\n" % name)
   				# write all ko records to a new fasta file
   				write_merged_df_to_2_fasta(input_df=temp_df, f_out=out_dir, out_name=name)
   ```

   

2. Partition in parallel
   ```
   #!/bin/bash
   
   # this is in the CPU server
   
   # active conda inside script
   temp=$(which conda)
   conda_path=$(echo ${temp%/*bin/conda})
   . ${conda_path}/etc/profile.d/conda.sh
   if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
           . ${conda_path}/etc/profile.d/conda.sh
   else
           echo "ERROR: conda path can't be corrected identified!!!"
           exit 1
   fi
   unset conda_path
   
   conda activate reproduce_kegg_profile
   
   
   
   target_dir=/data/shared_data/temp_shaopeng_process_kegg/PartitionByKO
   faa_file=/data/shared_data/temp_shaopeng_process_kegg/merged_KO_genes.faa
   
   cd ${target_dir}
   mkdir -p aa_per_KO
   
   # get KO list
   grep '>' ${faa_file} | rev | cut -d'|' -f1 | rev | sort -u > ${target_dir}/KO_list.txt
   
   # delete the double/multi count KOs records (both will be counted, no need to extra combination category)
   # may contain "|" or ","
   sed -i '/K.*\/K.*/d' ${target_dir}/KO_list.txt
   sed -i '/K.*,K.*/d' ${target_dir}/KO_list.txt
   
   # partition faa file using parallel
   distribute_single_ko () {
    ko=$1
    full_file=$2
    out_dir=$3
    echo "Processing ${ko}"
    grep -A 1 --no-group-separator -Fw ${ko} ${full_file} > ${out_dir}/${ko}.faa
   }
   
   export -f distribute_single_ko
   cat ${target_dir}/KO_list.txt | parallel -j 51 distribute_single_ko {} ${faa_file} ${target_dir}/aa_per_KO
   ```



3. sketch_aa_per_ko.sh
   ```
   #!/bin/bash
   
   # this is in the CPU server
   
   # active conda inside script
   temp=$(which conda)
   conda_path=$(echo ${temp%/*bin/conda})
   if [ -f ${conda_path}/etc/profile.d/conda.sh ]; then
           . ${conda_path}/etc/profile.d/conda.sh
   else
           echo "ERROR: conda path can't be corrected identified!!!"
           exit 1
   fi
   unset conda_path
   conda activate reproduce_kegg_profile
   
   target_dir=/data/shared_data/temp_shaopeng_process_kegg/sourmash_sketch
   faa_dir=/data/shared_data/temp_shaopeng_process_kegg/PartitionByKO/aa_per_KO
   
   # generate input csv file for batch usage
   cd ${faa_dir}
   echo "name,genome_filename,protein_filename" > ../KO_faa_for_sketch.csv
   for file in $(ls -1); do
    name=$(echo ${file%.faa})
    path=$(readlink -f ${file})
    echo "${name},,${path}" >> ../KO_faa_for_sketch.csv
   done
   
   inputcsv=$(readlink -f ../KO_faa_for_sketch.csv)
   
   cd ${target_dir}
   for scaleFactor in 100 1000; do
    for k_value in 5 7 11 15 25; do
     /usr/bin/time -av -o runlog_scale_${scaleFactor}_k_${k_value} sourmash sketch fromfile -p protein,k=${k_value},abund,scaled=${scaleFactor} -o ${target_dir}/KO_AA_sketch_scaled_${scaleFactor}_k_${k_value}.sig.zip ${inputcsv}
    done
   done
   
   date
   ```



4. build_sbt_lca.sh

   ```
   conda activate reproduce_kegg_profile
   
   sketch_folder=/data/shared_data/temp_shaopeng_process_kegg/sourmash_sketch
   sbt_folder=/data/shared_data/temp_shaopeng_process_kegg/SBT
   
   cd ${sbt_folder}
   for sig_db in $(ls ${sketch_folder}/KO_AA_sketch_*.sig.zip); do
    # didn't check parameter here as they are in the same type when creating them
    echo "Processing $sig_db"
    out_name=$(echo ${sig_db##*/})
    out_name=$(echo ${out_name%.sig.zip})
    echo "Out name is ${out_name}"
    /usr/bin/time -av -o runlog_SBT_${out_name} sourmash index ${out_name}.sbt.zip ${sig_db}
    unset sig_db out_name
   done
   
   date
   echo "done"
   ```
   
   
