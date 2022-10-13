# Helper class for the the EMDUniFrac associated functions
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import copy
epsilon = sys.float_info.epsilon
from matplotlib import pyplot, transforms

def functional_profile_to_EMDU_vector(functional_profile_file, EMDU_index_2_node, abundance_key='median_abund',
                                      normalize=True):
    """
    This function will take a sourmash functional profile and convert it to a form that can be used by EMDUniFrac

    :param functional_profile_file: csv file output from `sourmash gather`
    :param EMDU_index_2_node: dictionary that translates between the indices used by EMDUniFrac and the actual graph
    node names
    :param abundance_key: key in the functional profile that contains the abundance information
    :param normalize: whether to normalize the abundance information
    :return: vector P
    """
    # reverse the dictionary for convenience
    node_2_EMDU_index = {v: k for k, v in EMDU_index_2_node.items()}
    # import the functional profile
    df = pd.read_csv(functional_profile_file)
    # get the functional profile as a vector
    P = np.zeros(len(EMDU_index_2_node))
    for i, row in df.iterrows():
        try:
            abundance = float(row[abundance_key])
        except ValueError:
            raise ValueError(
                f"The abundance key {abundance_key} gives a value that is not a number: {row[abundance_key]}")
        try:
            if row['name'].startswith('ko:'):
                ko = row['name'].split(':')[-1]  # ko:K12345 case
            else:
                ko = row['name'].split('|')[-1].split(':')[-1]  # abc|xyz|lmnop|ko:K12345 case
            if ko not in node_2_EMDU_index:
                print(f"Warning: {ko} not found in EMDU index, skipping.")
            else:
                P[node_2_EMDU_index[ko]] = abundance
        except:
            raise Exception(f"Could not parse the name {row['name']} in the functional profile {functional_profile_file}")
    if P.sum() == 0:
        raise Exception(f"Functional profile {functional_profile_file} is empty! Perhaps try a different abundance key? "
                        f"You used {abundance_key}.")
    if normalize:
        P = P / P.sum()
    return P


def EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, P, Q):
    """
    (Z, F, diffab) = EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, P, Q)
    This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
    and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
    Returns the weighted Unifrac distance Z, the flow F, and the differential abundance vector diffab.
    The flow F is a dictionary with keys of the form (i,j) where F[(i,j)] == num means that in the calculation of the
    Unifrac distance, a total mass of num was moved from the node nodes_in_order[i] to the node nodes_in_order[j].
    The differential abundance vector diffab is	a dictionary with tuple keys using elements of Tint and values
    diffab[(i, j)] equal to the signed difference of abundance between the two samples restricted to the sub-tree
    defined by nodes_in_order(i) and weighted by the edge length lint[(i,j)].
    """
    num_nodes = len(nodes_in_order)
    F = dict()
    G = dict()
    diffab = dict()
    Z = 0
    w = np.zeros(num_nodes)
    pos = dict()
    neg = dict()
    for i in range(num_nodes):
        pos[i] = set([])
        neg[i] = set([])
    for i in range(num_nodes):
        if P[i] > 0 and Q[i] > 0:
            F[(i, i)] = np.minimum(P[i], Q[i])
        G[(i, i)] = P[i] - Q[i]
        if P[i] > Q[i]:
            pos[i].add(i)
        elif P[i] < Q[i]:
            neg[i].add(i)
        posremove = set()
        negremove = set()
        for j in pos[i]:
            for k in neg[i]:
                if (j not in posremove) and (k not in negremove):
                    val = np.minimum(G[(i, j)], -G[(i, k)])
                    if val > 0:
                        F[(j, k)] = np.minimum(G[(i, j)], -G[(i, k)])
                        G[(i, j)] = G[(i, j)] - val
                        G[(i, k)] = G[(i, k)] + val
                        Z = Z + (w[j] + w[k]) * val
                    if G[(i, j)] == 0:
                        posremove.add(j)
                    if G[(i, k)] == 0:
                        negremove.add(k)
        pos[i].difference_update(posremove)
        neg[i].difference_update(negremove)
        if i < num_nodes - 1:
            for j in pos[i].union(neg[i]):
                if (Tint[i], j) in G:
                    G[(Tint[i], j)] = G[(Tint[i], j)] + G[(i, j)]
                    diffab[(i, Tint[i])] = diffab[(i, Tint[i])] + G[(i, j)]  # Added to capture 'slack' in subtree JLMC
                else:
                    G[(Tint[i], j)] = G[(i, j)]
                    diffab[(i, Tint[i])] = G[(i, j)]  # Added to capture 'slack' in subtree JLMC
                w[j] = w[j] + lint[i, Tint[i]]
            if (i, Tint[i]) in diffab:
                # diffab[(i,Tint[i])] = lint[i,Tint[i]]*abs(diffab[(i,Tint[i])])  # Captures DiffAbund, scales by length of edge JLMC
                diffab[(i, Tint[i])] = lint[i, Tint[i]] * diffab[
                    (i, Tint[i])]  # Captures DiffAbund, scales by length of edge JLMC
            pos[Tint[i]] |= pos[i]
            neg[Tint[i]] |= neg[i]
    return (Z, F,
            diffab)  # The returned flow and diffab are on the basis nodes_in_order and is given in sparse matrix dictionary format. eg {(0,0):.5,(1,2):.5}


# This will return the EMDUnifrac distance only
def EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q):
    """
    (Z, diffab) = EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
    This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
    and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
    Returns the weighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
    F[(i,j)] == num means that in the calculation of the Unifrac distance, a total mass of num was moved from the node
    nodes_in_order[i] to the node nodes_in_order[j].
    """
    num_nodes = len(nodes_in_order)
    Z = 0
    diffab = dict()
    partial_sums = P - Q
    for i in range(num_nodes - 1):
        val = partial_sums[i]
        partial_sums[Tint[i]] += val
        if val != 0:
            diffab[(i, Tint[i])] = lint[i, Tint[i]] * val  # Captures diffab
        Z += lint[i, Tint[i]] * abs(val)
    return (Z, diffab)


# This will return the EMDUnifrac distance only
def EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q):
    """
    (Z, diffab) = EMDUnifrac_unweighted(Tint, lint, nodes_in_order, P, Q)
    This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
    and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
    Returns the unweighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
    F[(i,j)] == 1 means that in the calculation of the Unifrac distance, a total mass of 1 was moved from the node
    nodes_in_order[i] to the node nodes_in_order[j].
    """
    num_nodes = len(nodes_in_order)
    Z = 0
    diffab = dict()
    for i in range(num_nodes):
        if P[i] > 0:
            P[i] = 1
        if Q[i] > 0:
            Q[i] = 1
    partial_sums = P - Q
    for i in range(num_nodes - 1):
        val = partial_sums[i]
        partial_sums[Tint[i]] += val
        if val != 0:
            diffab[(i, Tint[i])] = lint[i, Tint[i]] * val  # Captures diffab
        Z += lint[i, Tint[i]] * abs(val)
    return Z, diffab


def EMDUnifrac_unweighted_flow(Tint, lint, nodes_in_order, P, Q):
    """
    (Z, F) = EMDUnifrac_unweighted_flow(Tint, lint, nodes_in_order, P, Q)
    This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
    and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
    Returns the unweighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
    F[(i,j)] == 1 means that in the calculation of the Unifrac distance, a total mass of 1 was moved from the node
    nodes_in_order[i] to the node nodes_in_order[j].
    """
    num_nodes = len(nodes_in_order)
    F = dict()
    G = dict()
    diffab = dict()
    Z = 0
    w = np.zeros(num_nodes)
    pos = dict()
    neg = dict()
    for i in range(num_nodes):
        pos[i] = set([])
        neg[i] = set([])
    for i in range(num_nodes):
        if P[i] > 0:
            P[i] = 1
        if Q[i] > 0:
            Q[i] = 1
        if P[i] > 0 and Q[i] > 0:
            F[(i, i)] = np.minimum(P[i], Q[i])
        G[(i, i)] = P[i] - Q[i]
        if P[i] > Q[i]:
            pos[i].add(i)
        elif P[i] < Q[i]:
            neg[i].add(i)
        posremove = set()
        negremove = set()
        for j in pos[i]:
            for k in neg[i]:
                if (j not in posremove) and (k not in negremove):
                    val = np.minimum(G[(i, j)], -G[(i, k)])
                    if val > 0:
                        F[(j, k)] = np.minimum(G[(i, j)], -G[(i, k)])
                        G[(i, j)] = G[(i, j)] - val
                        G[(i, k)] = G[(i, k)] + val
                        Z = Z + (w[j] + w[k]) * val
                    if G[(i, j)] == 0:
                        posremove.add(j)
                    if G[(i, k)] == 0:
                        negremove.add(k)
        pos[i].difference_update(posremove)
        neg[i].difference_update(negremove)
        if i < num_nodes - 1:
            for j in pos[i].union(neg[i]):
                if (Tint[i], j) in G:
                    G[(Tint[i], j)] = G[(Tint[i], j)] + G[(i, j)]
                    diffab[(i, Tint[i])] = diffab[(i, Tint[i])] + G[(i, j)]  # Added to capture 'slack' in subtree JLMC
                else:
                    G[(Tint[i], j)] = G[(i, j)]
                    diffab[(i, Tint[i])] = G[(i, j)]  # Added to capture 'slack' in subtree JLMC
                w[j] = w[j] + lint[i, Tint[i]]
            if (i, Tint[i]) in diffab:
                diffab[(i, Tint[i])] = lint[i, Tint[i]] * diffab[
                    (i, Tint[i])]  # Captures DiffAbund, scales by length of edge JLMC
            pos[Tint[i]] |= pos[i]
            neg[Tint[i]] |= neg[i]
    return (Z, F,
            diffab)  # The returned flow and diffab are on the basis nodes_in_order and is given in sparse matrix dictionary format. eg {(0,0):.5,(1,2):.5}


def plot_diffab(nodes_in_order, diffab, P_label, Q_label, plot_zeros=True, thresh=0):
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


def EMDUnifrac_weighted_plain(ancestors, edge_lengths, nodes_in_order, P, Q):
    """
    Z = EMDUnifrac_weighted(ancestors, edge_lengths, nodes_in_order, P, Q)
    This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
    and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
    Returns the weighted Unifrac distance Z.
    """
    num_nodes = len(nodes_in_order)
    Z = 0
    eps = 1e-8
    partial_sums = P - Q
    total_mass = 1
    for i in range(num_nodes - 1):
        val = partial_sums[i]
        if abs(val) > eps:  # if the value is big enough to care about (i.e. ignore zeros)
            # if np.sign(val) * np.sign(partial_sums[ancestors[i]]) == -1:  # if mass is getting destroyed
            # total_mass -= abs(val + partial_sums[ancestors[i]])  # keep track of how much was lost
            # print(total_mass)  # Can use this to break once the change is small enough
            partial_sums[ancestors[i]] += val
            Z += edge_lengths[i, ancestors[i]] * abs(val)
    return Z


def EMDUnifrac_group(ancestors, edge_lengths, nodes_in_order, rel_abund):
    eps = 1e-8
    num_nodes = len(nodes_in_order)
    num_samples = len(rel_abund)
    Z = np.zeros((num_samples, num_samples))
    partial_sums = np.zeros((num_samples, num_samples, num_nodes))
    for i in range(num_samples):
        for j in range(num_samples):
            partial_sums[i, j] = rel_abund[i] - rel_abund[j]
    for n in range(num_nodes - 1):
        for i in range(num_samples):
            for j in range(i, num_samples):
                val = partial_sums[i, j, n]
                if abs(val) > eps:  # only do the work if it's big enough
                    partial_sums[i, j, ancestors[n]] += val
                    Z[i, j] += edge_lengths[n, ancestors[n]] * abs(val)
    Z = Z + np.transpose(Z)
    return Z


def push_up_L2(P, Tint, lint, nodes_in_order):
    """
    Push the vector P up the tree, to prep the vector for L2 unifrac

    :param P: numpy vector
    :param Tint: dictionary of ancestors
    :param lint: dictionary of edge lengths
    :param nodes_in_order: list of nodes in post-ish order
    :return: numpy vector
    """
    P_pushed = copy.deepcopy(P)
    for i in range(len(nodes_in_order)-1):
        if lint[i, Tint[i]] == 0:
            lint[i, Tint[i]] = epsilon
        P_pushed[Tint[i]] += P_pushed[i] #push mass up
        P_pushed[i] *= np.sqrt(lint[i, Tint[i]])
    return P_pushed


def push_up_L1(P, Tint, lint, nodes_in_order):
    """
    Push the vector P up the tree, to prep the vector for L2 unifrac

    :param P: numpy vector
    :param Tint: dictionary of ancestors
    :param lint: dictionary of edge lengths
    :param nodes_in_order: list of nodes in post-ish order
    :return: numpy vector
    """
    P_pushed = copy.deepcopy(P)
    for i in range(len(nodes_in_order)-1):
        if lint[i, Tint[i]] == 0:
            lint[i, Tint[i]] = epsilon
        P_pushed[Tint[i]] += P_pushed[i] #push mass up
        P_pushed[i] *=lint[i, Tint[i]]
    return P_pushed

def EMD_L1_on_pushed(P, Q):
    """
    Computes the L1 earth movers distance on the **pushed up** vectors

    :param P: numpy vector
    :param Q: numpy vector
    :return: float
    """
    return np.sum(np.abs(P-Q))
