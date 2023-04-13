import matplotlib.pyplot as plt
import numpy as np


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

    