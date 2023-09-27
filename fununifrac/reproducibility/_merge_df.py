#merge multiple csv files with the same column names into one
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Combine infinite amount of files with the same column names into one")
parser.add_argument('-f', '--files', nargs='+')
parser.add_argument('-o', '--output', help="Output file name.")

def merge_files(*files):
    df = pd.DataFrame()
    for file in files:
        new_df = pd.read_csv(file)
        df = pd.concat([df, new_df], ignore_index=True)
    return df

if __name__ == '__main__':
    args = parser.parse_args()
    df = merge_files(*args.files)
    df.to_csv(args.output)
