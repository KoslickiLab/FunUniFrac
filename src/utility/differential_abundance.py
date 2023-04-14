import matplotlib.pyplot as plt
import numpy as np
import os.path
import pandas as pd


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


def plot_diffab(self, nodes_in_order, diffab, P_label, Q_label, plot_zeros=True, thresh=0):
        """
        plot_diffab(nodes_in_order, diffab, P_label, Q_label)
        Plots the differential abundance vector.

        :param nodes_in_order: list returned from parse_envs
        :param diffab: differential abundance vector (returned from one flavor of EMDUnifrac)
        :param P_label: label corresponding to the sample name for P (e.g. when calling EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q))
        :param Q_label: label corresponding to the sample name for P (e.g. when calling EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q))
        :param plot_zeros: flag (either True or False) that specifies if the zero locations should be plotted. Warning, if your tree is large and plot_zeros=True, this can cause a crash.
        :param thresh: only plot those parts of the diffab vector that are above thresh, specify everything else as zero
        :return: None (makes plot)
        """
        x = range(len(nodes_in_order))
        y = np.zeros(len(nodes_in_order))
        keys = diffab.keys()
        # for i in range(len(nodes_in_order)):
        #	for key in keys:
        #		if key[0] == i:
        #			y[i] = diffab[key]
        # Much faster way to compute this
        for key in keys:
            y[key[0]] = diffab[key]

        pos_loc = [x[i] for i in range(len(y)) if y[i] > thresh]
        neg_loc = [x[i] for i in range(len(y)) if y[i] < -thresh]
        zero_loc = [x[i] for i in range(len(y)) if -thresh <= y[i] <= thresh]
        if not pos_loc:
            raise Exception('Threshold too high! Please lower and try again.')
        if not neg_loc:
            raise Exception('Threshold too high! Please lower and try again.')

        pos_val = [y[i] for i in range(len(y)) if y[i] > thresh]
        neg_val = [y[i] for i in range(len(y)) if y[i] < -thresh]
        zero_val = [y[i] for i in range(len(y)) if -thresh <= y[i] <= thresh]

        # The following is to get the indicies in order. Basically, I iterate down both pos_loc and neg_loc simultaneously
        # and create new lists (pos_loc_adj and neg_loc_adj) that are in the same order as pos_loc and neg_loc, but whose
        # union of indicies is equal to range(len(pos_loc + neg_loc)). Simply to make things pretty
        if plot_zeros:
            pos_loc_adj = pos_loc
            neg_loc_adj = neg_loc
            zero_loc_adj = zero_loc
        else:
            pos_loc_adj = []
            neg_loc_adj = []
            tick_names = []

            # rename the indicies so they are increasing by 1
            pos_ind = 0
            neg_ind = 0
            it = 0
            while pos_ind < len(pos_loc) or neg_ind < len(neg_loc):
                if pos_ind >= len(pos_loc):
                    neg_loc_adj.append(it)
                    tick_names.append(nodes_in_order[neg_loc[neg_ind]])
                    it += 1
                    neg_ind += 1
                elif neg_ind >= len(neg_loc):
                    pos_loc_adj.append(it)
                    tick_names.append(nodes_in_order[pos_loc[pos_ind]])
                    it += 1
                    pos_ind += 1
                elif pos_loc[pos_ind] < neg_loc[neg_ind]:
                    pos_loc_adj.append(it)
                    tick_names.append(nodes_in_order[pos_loc[pos_ind]])
                    it += 1
                    pos_ind += 1
                elif pos_loc[pos_ind] > neg_loc[neg_ind]:
                    neg_loc_adj.append(it)
                    tick_names.append(nodes_in_order[neg_loc[neg_ind]])
                    it += 1
                    neg_ind += 1
                else:
                    print('Something went wrong')
                    break

        fig, ax = plt.subplots()

        markerline, stemlines, baseline = ax.stem(neg_loc_adj, neg_val)
        plt.setp(baseline, 'color', 'k', 'linewidth', 1)
        plt.setp(markerline, 'color', 'r')
        stemlines.set_linewidth(3)
        stemlines.set_color('r')
        #for i in range(len(neg_loc)):
        #    plt.setp(stemlines[i], 'linewidth', 3)
        #    plt.setp(stemlines[i], 'color', 'r')

        markerline, stemlines, baseline = ax.stem(pos_loc_adj, pos_val)
        plt.setp(baseline, 'color', 'k', 'linewidth', 1)
        plt.setp(markerline, 'color', 'b')
        stemlines.set_linewidth(3)
        stemlines.set_color('b')
        #for i in range(len(pos_loc)):
        #    plt.setp(stemlines[i], 'linewidth', 3)
        #    plt.setp(stemlines[i], 'color', 'b')

        if plot_zeros:
            markerline, stemlines, baseline = ax.stem(zero_loc, zero_val)
            plt.setp(baseline, 'color', 'k', 'linewidth', 1)
            plt.setp(markerline, 'color', 'k')
            stemlines.set_linewidth(3)
            stemlines.set_color('k')
            #for i in range(len(zero_loc)):
            #    plt.setp(stemlines[i], 'linewidth', 3)
            #    plt.setp(stemlines[i], 'color', 'k')

        plt.ylabel('DiffAbund', fontsize=16)
        plt.gcf().subplots_adjust(right=0.93, left=0.15)
        # If you want the zeros plotted, label EVERYTHING, otherwise just label the things that are there...
        if plot_zeros:
            plt.xticks(x, nodes_in_order, rotation='vertical', fontsize=14)
        else:
            # tick_names = [nodes_in_order[i] for i in pos_loc] + [nodes_in_order[i] for i in neg_loc]  # Don't need this with new code
            plt.xticks(range(len(pos_loc_adj + neg_loc_adj)), tick_names, rotation='vertical', fontsize=14)

        plt.subplots_adjust(bottom=0.3, top=.93)
        plt.text(plt.xticks()[0][-1] + 0.1, max(pos_val), P_label, rotation=90, horizontalalignment='center',
                verticalalignment='top', multialignment='center', color='b', fontsize=14)
        plt.text(plt.xticks()[0][-1] + 0.1, min(neg_val), Q_label, rotation=90, horizontalalignment='center',
                verticalalignment='bottom', multialignment='center', color='r', fontsize=14)
        plt.show()

    