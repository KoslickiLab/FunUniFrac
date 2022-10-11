import subprocess
from os.path import exists
import tempfile
import numpy as np
import os
from scipy import sparse

def make_test_data():
    distances_file = "test_data/small_pairwise_distances.npy"
    labels_file = "test_data/small_pairwise_distances.npy.labels.txt"
    with open(labels_file, 'w') as fid:
        fid.write("d\n")
        fid.write("e\n")
        fid.write("f\n")
        fid.write("g\n")
        fid.write("h\n")
    mat = np.array([[0, 2, 4, 4, 7],
                    [2, 0, 4, 4, 7],
                    [4, 4, 0, 2, 7],
                    [4, 4, 2, 0, 7],
                    [7, 7, 7, 7, 0]])
    np.save(distances_file, mat)


def test_small_edge_list():
    """
    Uses a complete binary tree with 4 leaf nodes, all branch lengths set to 1
    :return: None
    """
    edge_list = "test_data/small_edge_list.txt"
    distances_file = "test_data/small_pairwise_distances.npy"
    out_dir = "test_data"
    brite = "ko00001"
    out_file = os.path.join(out_dir, f"{brite}_{os.path.basename(distances_file)}_A.npz")
    # delete the output file if it exists
    if exists(out_file):
        os.remove(out_file)
    cmd = f"python ../scripts/graph_to_path_matrix.py -e {edge_list} -d {distances_file} -o {out_dir} -b {brite}"
    subprocess.run(cmd, shell=True, check=True)
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
