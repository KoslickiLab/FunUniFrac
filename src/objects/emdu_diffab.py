import dataclasses
import os.path
import numpy as np

    
def convert_diffab_array_to_df(self, diffabs, nodes_in_order, file_basis):
        """
        Converts the differential abundance array to a 3D dataframe.
        DO NOT USE EXCEPT FOR VERY SMALL DATA SETS

        :param diffabs: 3D numpy array: [i,j,:] is the differential abundance vector between sample i and j
        :param nodes_in_order: list of nodes in post-ish order
        :param file_basis: basis for the file names
        :return: pandas dataframe
        """
        # FIXME: this is extremely slow
        file_basis = [os.path.basename(x) for x in file_basis]
        num_samples = diffabs.shape[0]
        if num_samples != len(file_basis):
            raise ValueError('Number of samples does not match number of file names')
        if num_samples != diffabs.shape[1]:
            raise ValueError('Differential abundance array is not square')
        if diffabs[0, 0, :].shape[0] != len(nodes_in_order):
            raise ValueError('Differential abundance array does not match nodes in order')
        diffabs_dict = {}
        for i, file in enumerate(file_basis):
            diffabs_dict[file] = {}
            for j, file2 in enumerate(file_basis):
                diffabs_dict[file][file2] = dict()
                for k, node in enumerate(nodes_in_order):
                    diffabs_dict[file][file2][node] = diffabs[i, j, k]
        diffabs_df = pd.DataFrame.from_dict(diffabs_dict)
        return diffabs_df


class DiffabArrayIndexer:
    """
    This class is used to index the differential abundance array. It is used to make the indexing more readable.
    """
    def __init__(self, diffabs, nodes_in_order, file_basis, EMDU_index_2_node):
        self.diffabs = diffabs
        self.node_to_index = {node: i for i, node in enumerate(nodes_in_order)}
        self.file_to_index = {file: i for i, file in enumerate(file_basis)}
        # reverse the EMDU_index_2_node dictionary
        self.node_to_EMDU_index = {v: k for k, v in EMDU_index_2_node.items()}

    def get_diffab(self, files1, files2):
        """
        Get the differential abundance vector between two sets of files. Sets can be singletons.

        :param files1: list of files, or single file
        :param files2: list of files, or single file
        :return: numpy array
        """
        if isinstance(files1, str):
            files1 = [files1]
        if isinstance(files2, str):
            files2 = [files2]
        indices1 = np.array([self.file_to_index[file] for file in files1])
        indices2 = np.array([self.file_to_index[file] for file in files2])
        #return self.diffabs[indices1, indices2, :]
        return self.diffabs[np.ix_(indices1, indices2, np.arange(self.diffabs.shape[2]))]

    def get_diffab_for_node(self, files1, files2, nodes):
        """
        Get the differential abundance vector between two sets of files. Sets can be singletons. So can nodes.

        :param files1: list of files, or single file
        :param files2: list of files, or single file
        :param nodes: node(s) to get the differential abundance for. Eg. K01234
        :return: numpy array
        """
        if isinstance(files1, str):
            files1 = [files1]
        if isinstance(files2, str):
            files2 = [files2]
        if isinstance(nodes, str) or isinstance(nodes, int):
            nodes = [nodes]
        indices1 = [self.file_to_index[file] for file in files1]
        indices2 = [self.file_to_index[file] for file in files2]
        indices3 = [self.node_to_index[self.node_to_EMDU_index[node]] for node in nodes]
        #return self.diffabs[np.ix_(indices1, indices2, indices3)]
        return self.diffabs[indices1, :, :][:, indices2, :][:, :, indices3]


