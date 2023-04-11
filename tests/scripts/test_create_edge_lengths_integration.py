import subprocess
import numpy as np
from scipy import sparse
import pandas as pd
import data


def test_small_edge_lengths():
    """
    Uses a complete binary tree with 4 leaf nodes, all branch lengths set to 1
    :return: None
    """
    edge_list = "small_edge_list.txt"
    distances_file = "small_pairwise_distances.npy"
    out_file = "test_output/small_edge_list_with_lengths.txt"
    brite = "ko00001"
    A_file = f"{brite}_{distances_file}_A.npz"
    cmd = f"python ../scripts/create_edge_lengths.py -e {edge_list} -d {distances_file} -o {out_file} -b {brite} -A {A_file} -n 10 -f 2 -r 100 --force --distance"
    subprocess.run(cmd, shell=True, check=True, executable='/usr/bin/bash')
    
    # check that the output file is correct
    df = pd.read_csv(data.get_data_abspath(out_file), sep="\t", header=0)
    # get the edge lengths
    edge_lengths = df["edge_length"].values
    # check if the edge lengths are all close to 1
    assert np.allclose(edge_lengths, np.ones_like(edge_lengths), atol=1e-2)
