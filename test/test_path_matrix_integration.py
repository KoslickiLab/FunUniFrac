import subprocess
from os.path import exists
import tempfile
import numpy as np
import os
from scipy import sparse

def test_small_edge_list():
    """
    Uses a complete binary tree with 4 leaf nodes, all branch lengths set to 1
    :return: None
    """
    edge_list = "test_data/small_edge_list.txt"
    distances_file = "test_data/small_pairwise_distances.npy"
    out_dir = "."
    brite = "ko00001"
    cmd = f"python ../scripts/graph_to_path_matrix.py -e {edge_list} -d {distances_file} -o {out_dir} -b {brite}"
    subprocess.run(cmd, shell=True, check=True)
    out_file = f"test_data/{brite}_{os.path.basename(distances_file)}_A.npz"
    # check that the output file exists
    assert exists(out_file)
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
