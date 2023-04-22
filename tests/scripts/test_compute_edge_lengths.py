import subprocess
import numpy as np
from scipy import sparse
import pandas as pd
import data
from src.utility import constant


def test_small_edge_lengths():
    """
    Uses a complete binary tree with 4 leaf nodes, all branch lengths set to 1
    :return: None
    """
    edge_list = data.get_data_abspath("small_edge_list.txt")
    distances_file = data.get_data_abspath("small_pairwise_distances.npy")
    out_dir = data.get_data_abspath("test_output")
    out_id = "test"
    brite = "ko00001"
    cmd = f"python ../fununifrac/compute_edges.py -e {edge_list} -d {distances_file} -o {out_dir} -i {out_id} -b {brite} -n 10 -f 2 -r 100 --distance"
    subprocess.run(cmd, shell=True, check=True)
    
    # check that the output file is correct
    out_file = f"{out_dir}/{constant.EDGE_LIST_OUT__FILE_NAME.format(out_id)}"
    df = pd.read_csv(out_file, sep="\t", header=0)
    # get the edge lengths
    edge_lengths = df["edge_length"].values
    # check if the edge lengths are all close to 1
    assert np.allclose(edge_lengths, np.ones_like(edge_lengths), atol=1e-2)
