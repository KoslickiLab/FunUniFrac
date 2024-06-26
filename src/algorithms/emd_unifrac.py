# Helper class for the the EMDUniFrac associated functions
import abc
import numpy as np
import sys
import copy
from src.objects.func_tree import FuncTreeEmduInput
from src.objects.profile_vector import FuncProfileVector, get_L1, get_L1_diffab, get_L2, get_L2_diffab
from src.objects.emdu_out import UnifracDistance, EmdFlow, DifferentialAbundance
from typing import Tuple
epsilon = sys.float_info.epsilon
from itertools import combinations
from blist import blist
import sparse
    

class EarthMoverDistanceUniFracSolver:
    def __init__(self) -> None:
        pass    

    def push_up_L2(self, P: FuncProfileVector, input: FuncTreeEmduInput)->FuncProfileVector:
        P_pushed = copy.deepcopy(P)
        Tint, lint, nodes_in_order = input.Tint, input.lint, input.basis
        for i in range(len(nodes_in_order)-1):
            if lint[i, Tint[i]] == 0:
                lint[i, Tint[i]] = epsilon
            P_pushed[Tint[i]] += P_pushed[i] #push mass up
            P_pushed[i] *= np.sqrt(lint[i, Tint[i]])
        return P_pushed


    def push_up_L1(self, P: FuncProfileVector, input: FuncTreeEmduInput):
        P_pushed = copy.deepcopy(P)
        Tint, lint, nodes_in_order = input.Tint, input.lint, input.basis
        for i in range(len(nodes_in_order)-1):
            if lint[i, Tint[i]] == 0:
                lint[i, Tint[i]] = epsilon
            P_pushed[Tint[i]] += P_pushed[i] #push mass up
            P_pushed[i] *= lint[i, Tint[i]]
        return P_pushed


    def pairwise_computation(self, Ps: dict[str, FuncProfileVector], vector_names, input: FuncTreeEmduInput, make_diffab, is_L2):
        dists = np.zeros((len(vector_names), len(vector_names)))
        i_coords = blist([])
        j_coords = blist([])
        k_coords = blist([])
        data_vals = blist([])
        diffab_dims = (len(vector_names), len(vector_names), len(input.basis))
        for i, j in combinations(range(len(vector_names)), 2):
            if not make_diffab:
                if not is_L2:
                    dists[i, j] = dists[j, i] = get_L1(Ps[vector_names[i]], Ps[vector_names[j]])
                else:
                    dists[i, j] = dists[j, i] = get_L2(Ps[vector_names[i]], Ps[vector_names[j]])
            else:
                if not is_L2:
                    Z = get_L1(Ps[vector_names[i]], Ps[vector_names[j]])
                    diffab = get_L1_diffab(Ps[vector_names[i]], Ps[vector_names[j]])
                else:
                    Z = get_L2(Ps[vector_names[i]], Ps[vector_names[j]])
                    diffab = get_L2_diffab(Ps[vector_names[i]], Ps[vector_names[j]])
                dists[i, j] = dists[j, i] = Z
                nonzero_diffab_locs = np.nonzero(diffab)[0]
                i_coords.extend([i] * len(nonzero_diffab_locs))
                j_coords.extend([j] * len(nonzero_diffab_locs))
                k_coords.extend(nonzero_diffab_locs)
                data_vals.extend(diffab[nonzero_diffab_locs])
        if make_diffab:
            coords = [i_coords, j_coords, k_coords]
            diffabs_sparse = sparse.COO(coords, data_vals, shape=diffab_dims)
            # add the negative transpose ([[0, 1], [0, 0]] -> [[0, 1], [-1, 0]])
            diffabs_sparse = diffabs_sparse - diffabs_sparse.transpose(axes=(1, 0, 2))
        else:
            diffabs_sparse = None
        return dists, diffabs_sparse 

    ########################################################################################################################
    # Not used. Moved to pushup paradigm
    ########################################################################################################################
    def solve(self, 
              input: FuncTreeEmduInput, 
              P: FuncProfileVector, 
              Q: FuncProfileVector, 
              weighted: bool):
        """(Z, diffab) = EMDUnifrac_weighted(Tint, lint, nodes_in_order, P, Q)
        This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
        and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
        Returns the weighted Unifrac distance Z and the flow F. The flow F is a dictionary with keys of the form (i,j) where
        F[(i,j)] == num means that in the calculation of the Unifrac distance, a total mass of num was moved from the node
        nodes_in_order[i] to the node nodes_in_order[j].
        """
        Tint, lint, nodes_in_order = input.Tint, input.lint, input.basis
        num_nodes = len(nodes_in_order)
        Z: UnifracDistance = 0
        diffab: DifferentialAbundance = dict()
        
        if not weighted:
            for i in range(num_nodes):
                ################################################################################
                # Unweighted
                ################################################################################
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

    ########################################################################################################################
    # Not used. Moved to pushup paradigm
    ########################################################################################################################
    def solve_with_flow(self, 
                        input: FuncTreeEmduInput, 
                        P: FuncProfileVector, 
                        Q: FuncProfileVector, 
                        weighted: bool):
        """ (Z, F, diffab) = EMDUnifrac_weighted_flow(Tint, lint, nodes_in_order, P, Q)
        This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
        and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
        Returns the weighted Unifrac distance Z, the flow F, and the differential abundance vector diffab.
        The flow F is a dictionary with keys of the form (i,j) where F[(i,j)] == num means that in the calculation of the
        Unifrac distance, a total mass of num was moved from the node nodes_in_order[i] to the node nodes_in_order[j].
        The differential abundance vector diffab is	a dictionary with tuple keys using elements of Tint and values
        diffab[(i, j)] equal to the signed difference of abundance between the two samples restricted to the sub-tree
        defined by nodes_in_order(i) and weighted by the edge length lint[(i,j)].
        """
        Tint, lint, nodes_in_order = input.Tint, input.lint, input.basis
        num_nodes = len(nodes_in_order)
        F: EmdFlow = dict()
        G = dict()
        Z: UnifracDistance = 0
        diffab: DifferentialAbundance = dict()
        w = np.zeros(num_nodes)
        pos = dict()
        neg = dict()
        for i in range(num_nodes):
            pos[i] = set([])
            neg[i] = set([])
        for i in range(num_nodes):
            ################################################################################
            # Unweighted
            ################################################################################
            if not weighted:
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
        return (Z, F, diffab)
    # @abc.abstractmethod
    # def EMDUnifrac_weighted_plain(ancestors, 
    #                               edge_lengths, 
    #                               nodes_in_order, 
    #                               P, 
    #                               Q):
    #     """
    #     Z = EMDUnifrac_weighted(ancestors, edge_lengths, nodes_in_order, P, Q)
    #     This function takes the ancestor dictionary Tint, the lengths dictionary lint, the basis nodes_in_order
    #     and two probability vectors P and Q (typically P = envs_prob_dict[samples[i]], Q = envs_prob_dict[samples[j]]).
    #     Returns the weighted Unifrac distance Z.
    #     """
    #     raise NotImplementedError()
    
    # @abc.abstractmethod
    # def EMDUnifrac_group(ancestors, 
    #                      edge_lengths, 
    #                      nodes_in_order, 
    #                      rel_abund):
    #     raise NotImplementedError()

    # def EMDUnifrac_weighted_plain(self, ancestors, edge_lengths, nodes_in_order, P, Q):
    #     num_nodes = len(nodes_in_order)
    #     Z = 0
    #     eps = 1e-8
    #     partial_sums = P - Q
    #     total_mass = 1
    #     for i in range(num_nodes - 1):
    #         val = partial_sums[i]
    #         if abs(val) > eps:  # if the value is big enough to care about (i.e. ignore zeros)
    #             # if np.sign(val) * np.sign(partial_sums[ancestors[i]]) == -1:  # if mass is getting destroyed
    #             # total_mass -= abs(val + partial_sums[ancestors[i]])  # keep track of how much was lost
    #             # print(total_mass)  # Can use this to break once the change is small enough
    #             partial_sums[ancestors[i]] += val
    #             Z += edge_lengths[i, ancestors[i]] * abs(val)
    #     return Z

    # def EMDUnifrac_group(self, ancestors, edge_lengths, nodes_in_order, rel_abund):
    #     eps = 1e-8
    #     num_nodes = len(nodes_in_order)
    #     num_samples = len(rel_abund)
    #     Z = np.zeros((num_samples, num_samples))
    #     partial_sums = np.zeros((num_samples, num_samples, num_nodes))
    #     for i in range(num_samples):
    #         for j in range(num_samples):
    #             partial_sums[i, j] = rel_abund[i] - rel_abund[j]
    #     for n in range(num_nodes - 1):
    #         for i in range(num_samples):
    #             for j in range(i, num_samples):
    #                 val = partial_sums[i, j, n]
    #                 if abs(val) > eps:  # only do the work if it's big enough
    #                     partial_sums[i, j, ancestors[n]] += val
    #                     Z[i, j] += edge_lengths[n, ancestors[n]] * abs(val)
    #     Z = Z + np.transpose(Z)
    #     return Z