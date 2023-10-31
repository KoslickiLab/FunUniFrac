import pandas as pd
import os
import glob
import numpy as np

def merge_files(file_dir, file_pattern, abund_col, ko_col):
    file_pattern = os.path.join(file_dir, file_pattern)
    fun_files = glob.glob(file_pattern)
    KOs = set()
    for f in fun_files:
        df = pd.read_csv(f)
        KOs.update(df[ko_col])
    KO_list = list(KOs)
    combined_df = pd.DataFrame(index=KO_list)
    KO_index_dict = {ko:i for i, ko in enumerate(KO_list)}
    for f in fun_files:
        sample_name = os.path.basename(f).split('.')[0]
        df = pd.read_csv(f)
        vector = np.zeros(len(KO_list))
        for i, row in df.iterrows():
            ko = row[ko_col]
            index_in_all_kos = KO_index_dict[ko]
            vector[index_in_all_kos] = row[abund_col]
        combined_df[sample_name] = vector
    return combined_df

combined_df = merge_files('/home/grads/wjw5274/FunUniFrac/fununifrac/reproducibility/data/skin_vs_gut', '*.csv', 'f_unique_weighted', 'name')
print(combined_df)
