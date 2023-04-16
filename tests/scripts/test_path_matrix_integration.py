import subprocess
import numpy as np
import os
from scipy import sparse
import data
import src.utility.constant as constant


def test_small_edge_list():
    """
    Uses a complete binary tree with 4 leaf nodes, all branch lengths set to 1
    :return: None
    """
    edge_file = "small_edge_list.txt"
    distances_file = "small_pairwise_distances.npy"
    brite = "ko00001"
    out_dir = "test_output"
    cmd = f"python ../scripts/create_edge_matrix.py -e {edge_file} -d {distances_file} -o {out_dir} -b {brite}"
    subprocess.run(cmd, shell=True, check=True)
    
    # check that the output file exists
    out_file = data.get_data_abspath(f"{out_dir}/{constant.A_MATRIX__FILE_NAME}") 
    # check that the output file is correct
    A = sparse.load_npz(out_file)
    assert A.shape == (4**2, 6+1)
    # use the hand-calculated A matrix to check that the output is correct
    A_correct = np.array([[0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 1, 1, 0, 0],
                          [0, 1, 1, 1, 0, 1, 0],
                          [0, 1, 1, 1, 0, 0, 1],
                          [0, 0, 0, 1, 1, 0, 0],
                          [0, 0, 0, 0, 0, 0, 0],
                          [0, 1, 1, 0, 1, 1, 0],
                          [0, 1, 1, 0, 1, 0, 1],
                          [0, 1, 1, 1, 0, 1, 0],
                          [0, 1, 1, 0, 1, 1, 0],
                          [0, 0, 0, 0, 0, 0, 0],
                          [0, 0, 0, 0, 0, 1, 1],
                          [0, 1, 1, 1, 0, 0, 1],
                          [0, 1, 1, 0, 1, 0, 1],
                          [0, 0, 0, 0, 0, 1, 1],
                          [0, 0, 0, 0, 0, 0, 0]])
    assert np.allclose(A.toarray(), A_correct)
