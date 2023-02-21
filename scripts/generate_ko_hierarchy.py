#!/usr/bin/env python
import os, sys
import argparse
import pandas as pd
import requests
import re

def organize_hierarchy(hiearchy_json, regex='^\w?\w?\d{5} ', prefix='', stop_level=None):

    def _iterate_multidimensional(res, hiearchy_json, res_list):
        if isinstance(hiearchy_json,dict):
            for k,v in hiearchy_json.items():
                if k == 'name':
                    temp_name = re.sub(' \[.*\]','',hiearchy_json['name'])
                    res += f"|{temp_name}"
                    if re.search(regex, hiearchy_json['name']) is not None:
                        res_list += [res]
                elif k == 'children':
                    for elem in hiearchy_json['children']:
                        _iterate_multidimensional(res, elem, res_list)
        else:
            # FIXME: self is not defined
            raise ValueError(f"{hiearchy_json} is not dictionary")
            

    res_list = []
    _iterate_multidimensional('', hiearchy_json, res_list)
    if stop_level is None:
        res_dict = {prefix+string.split('|')[-1].split(' ')[0]:'|'.join(string.split('|')[1:]) for string in res_list}
        return res_dict
    else:
        res_dict = {prefix+string.split('|')[-1].split(' ')[0].split('\t')[0]:'|'.join(string.split('|')[1:(stop_level+1)]) for string in res_list}
        return res_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", required=False, type=str, default='./', help="Path to which the output result is saved")
    parser.add_argument("--brite_list", required=False, type=str, nargs="*", default=['ko00001'], help="A list of BRITE IDs (eg. ko00001) that are associated with KO to extract their subtree. Default: ['ko00001']")
    parser.add_argument("--saved_kos", required=False, type=str, default=None, help="Path to a file that contains a list of ko IDs reserved in a hiearchical substree")
    args = parser.parse_args()
    
    outdir = args.outdir
    brite_list = args.brite_list
    saved_kos = args.saved_kos

    KEGG_api_link = 'http://rest.kegg.jp'
            
    # get brite table
    link = f"{KEGG_api_link}/list/brite"
    res = requests.get(link)
    if res.status_code == 200:
        brite_table = pd.DataFrame([x.split('\t') for x in res.text.split('\n') if x.split('\t')[0]])
        brite_table.columns = ['kegg_brite_id','desc']
    else:
        print(f"Error: Fail to download KEGG brite information from {link}", flush=True)
        exit()

    ## check parameter
    if len(brite_list) == 0:
        print(f"No brite ID is given via the paramter 'brite_list'", flush=True)
        exit()

    diff_brites = set(brite_list).difference(set(brite_table['kegg_brite_id'].str.replace('br:','')))
    if len(diff_brites) > 0:
        print(f"The given brite ids '{list(diff_brites)}' is not in the allowable brite list:\n{list(brite_table['kegg_brite_id'])}", flush=True)
        exit()

    # set up identifier mapping
    id_mapping = {f"{re.sub('^[a-z]*:[a-z]*','',x[0])} {x[1]}":x[0].split(':')[1] for x in brite_table.to_numpy()}

    ## get brite hierarchical substree and process hierarchy
    print("Get brite hierarchical substree and process hierarchy")
    brite_list = ['br:'+str(x) for x in brite_list]
    for index, brite_id in enumerate(brite_list):
        substree_dict = dict()
        print(f"Processing brite id {brite_id}", flush=True)
        check_link = f"{KEGG_api_link}/get/{brite_id}"
        res = requests.get(check_link)
        if res.status_code == 200:
            m = re.search('K\d{5}', res.text)
        else:
            print(F"ERROR: Fail to download KEGG brite information from {check_link} and thus skip it.", flush=True)
            continue
        if m is None:
            print(f"WARNING: Brite ID {brite_id} doesn't contain KO ids and thus skip it.", flush=True)
            continue
        link = f"{KEGG_api_link}/get/{brite_id}/json"
        json_res = requests.get(link)
        if json_res.status_code == 200:
            temp_dict = organize_hierarchy(json_res.json(), regex='K\d{5}')
            for key in temp_dict:
                if key in substree_dict:
                    substree_dict[key] += [temp_dict[key]]
                else:
                    substree_dict[key] = [temp_dict[key]]
        else:
            print(f"Error: Fail to download brite information from {link}")

        if saved_kos is not None:
            saved_kos = pd.read_csv(saved_kos, sep='\t', header=None)
            saved_kos = list(saved_kos[0])
            diff_kos = list(set(saved_kos).difference(set(substree_dict.keys())))
            if len(diff_kos) > 0:
                print(f"The following KOs in your given KO list are not found in hierarchical tree:\n {diff_kos}", flush=True)
            substree_dict = {ko_id:hierarchy_list for ko_id, hierarchy_list in substree_dict.items() if ko_id in saved_kos}
            outfile = f"kegg_ko_edge_df_filtered_{brite_id.replace('br:','')}.txt"
        else:
            outfile = f"kegg_ko_edge_df_{brite_id.replace('br:','')}.txt"

        # convert hierarchy to edge list
        ko_edge_list = []
        brite_root_id_list = []
        for ko_id, hierarchy_list in substree_dict.items():
            for item in hierarchy_list:
                temp_list = [id_mapping.get(x,x) for x in item.split('|')[:-1]]
                brite_root_id_list += [temp_list[0]]
                temp_list += [ko_id]
                ko_edge_list += [(temp_list[index-1], temp_list[index]) for index in range(1, len(temp_list))]
                
        kegg_ko_edge_df = pd.DataFrame(ko_edge_list)
        kegg_ko_edge_df.columns = ['parent','child']
        kegg_ko_edge_df = kegg_ko_edge_df.drop_duplicates().reset_index(drop=True)
        kegg_ko_edge_df.to_csv(os.path.join(outdir, outfile), sep='\t', index=None)
