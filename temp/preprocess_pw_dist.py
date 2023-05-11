from Kegg_tree import KeggTree, get_KeggTree_from_edgelist
import time

pw_dist_file = 'KOs_sketched_scaled_10_compare_5'
label_file = 'labels_clean.txt'
edge_list = 'kegg_ko00001_no_edge_lengths.txt'

kegg_tree = get_KeggTree_from_edgelist(edge_list, edge_length=False)
start_time = time.time()
kegg_tree.preprocess_pw_dist(pw_dist_file, label_file)
kegg_tree.write_pw_dist('pw_dist_dict.txt')
end_time = time.time()
print(f"time to process is {(end_time-start_time)/60} seconds")