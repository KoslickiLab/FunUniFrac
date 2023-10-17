import argparse
import pandas as pd
import random
import numpy as np

ROW_NUM = 1000
COL_NUM = 100

#get hold of a template
def parsearg():
    parser = argparse.ArgumentParser(description="This script generates simulated data")
    parser.add_argument("-f", "--template_file", type=str, help="Template file")
    parser.add_argument("-o", "--out_dir", type=str, help="Output direcctory")
    return parser.parse_args()


def main():
    args = parsearg()

    df = pd.read_csv(args.template_file, index_col='name')
    print(f"No. of samples: {len(df.columns)}")
    print(f"No. of KO: {len(df.index)}")
    df = df[df.columns[:COL_NUM]]
    df = df[:ROW_NUM]
    meta_dict = dict()
    for i in range(100):
        for percentage in [0.5, 0.75, 0.9]:
            partition = int(ROW_NUM * percentage)
            file_name = f"{args.out_dir}/sim_sample_{partition}_i.csv"
            for col in df.columns[:50]:
                vector = np.zeros(ROW_NUM)
                vector[:partition] =[random.random() for i in range(partition)]
                df[col] = vector
                meta_dict[col] = "environment 1"
            for col in df.columns[50:]:
                vector = np.zeros(ROW_NUM)
                vector[ROW_NUM-partition:] = [random.random() for i in range(partition)]
                df[col] = vector
                meta_dict[col] = "environment 2"
            df.to_csv(file_name)
    meta_df = pd.DataFrame(meta_dict.items(), columns=['sample', 'env'])
    meta_df.to_csv(f"{args.out_dir}/simulated_metadata.csv")



if __name__ == '__main__':
    main()
