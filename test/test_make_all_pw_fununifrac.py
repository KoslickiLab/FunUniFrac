import subprocess
from os.path import exists
import tempfile
import numpy as np
import os
from scipy import sparse
import pandas as pd


def test_small_edge_lengths():
    """
    Uses a complete binary tree with 4 leaf nodes, all branch lengths set to 1
    :return: None
    """
    edge_list = "test_data/small_edge_list_with_lengths.txt"
    directory = "test_data"
    profile1 = "test_data/small_sim_10_KOs_gather.csv"
    profile2 = "test_data/small_sim2_10_KOs_gather.csv"
    out_file = "test_data/pairwise_dists.npy"
    brite = "ko00001"
    cmd = f"python ../scripts/make_all_pw_fununifrac.py -e {edge_list} -d " \
          f"{directory} -o {out_file} --force -b {brite}"
    res = subprocess.run(cmd, shell=True, check=True)
    # check that the output file exists
    assert exists(out_file)
    # check that the output file is correct
    pw_dists = np.load(out_file)
    # This has been computed by hand
    known_pw_dists = np.array([[0.0, 10/21], [10/21, 0.0]])
    assert np.allclose(pw_dists, known_pw_dists, atol=1e-3)
    # check that the basis file has been created and is correct
    basis_file = f"{out_file}.basis.txt"
    assert exists(basis_file)
    # open the basis file
    with open(basis_file, "r") as f:
        basis = os.path.basename(f.readline().strip())
        assert basis == os.path.basename(profile2)
        basis = os.path.basename(f.readline().strip())
        assert basis == os.path.basename(profile1)
        # check that the end of the file has been reached
        assert f.readline() == ""

    # now also check the diffabs
    cmd = f"python ../scripts/make_all_pw_fununifrac.py -e {edge_list} -d " \
              f"{directory} -o {out_file} --force -b {brite} --diffab"
    diffab_out_file = "test_data/pairwise_dists.npy.diffab.npy"
    res = subprocess.run(cmd, shell=True, check=True)
    # check that the output file exists
    assert exists(out_file)
    # check that the output file is correct
    pw_dists = np.load(out_file)
    # This has been computed by hand
    known_pw_dists = np.array([[0.0, 10/21], [10/21, 0.0]])
    assert np.allclose(pw_dists, known_pw_dists, atol=1e-3)
    # check that the basis file has been created and is correct
    basis_file = f"{out_file}.basis.txt"
    assert exists(basis_file)
    # open the basis file
    with open(basis_file, "r") as f:
        basis = os.path.basename(f.readline().strip())
        assert basis == os.path.basename(profile2)
        basis = os.path.basename(f.readline().strip())
        assert basis == os.path.basename(profile1)
        # check that the end of the file has been reached
        assert f.readline() == ""
    # Check that the diffabs file is correct
    assert exists(diffab_out_file)
    diffabs = np.load(diffab_out_file)
    known_diffab = np.array([-1/14, -1/21, -5/42, 5/42, 0, 5/42, 0])
    zeros = np.zeros_like(known_diffab)
    assert np.allclose(diffabs[0, 1, :], known_diffab, atol=1e-3)
    assert np.allclose(diffabs[1, 0, :], known_diffab, atol=1e-3)
    assert np.allclose(diffabs[0, 0, :], zeros, atol=1e-3)
    assert np.allclose(diffabs[1, 1, :], zeros, atol=1e-3)
