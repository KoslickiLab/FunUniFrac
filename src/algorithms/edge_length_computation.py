#!/usr/bin/env python
import networkx as nx
import numpy as np
from scipy import sparse
import multiprocessing
try:
    from blist import blist
except ModuleNotFoundError:
    print("Warning: Could not import blist. Please install blist to speed up the path matrix calculation.")
from src.algorithms.lp_edge_length import get_descendant
from src.base.func_tree import FuncTree
from src.base.pairwise_dist import PairwiseDistance
from itertools import repeat
from scipy.optimize import lsq_linear


class EdgeLengthSolver:
    def __init__(self) -> None:
        pass

    def shortest_path_parallel(self, args):
        def map_func(node_i, G):
            if node_i in G:
                return {node_i: nx.single_source_dijkstra_path(G, node_i)}
            else:
                return {}
        return map_func(*args)

    def get_A_matrix(self, tree: FuncTree, pairwise_distances):
        ########################################################################################################################
        # Let's do the following: since I've already computed all pairwise distances, we can just make a large
        # least squares problem fitting the tree distances to the pairwise distances
        # Let's get the matrix describing which edges are traversed between all pairs of nodes
        # This is a sparse matrix, so we'll need to use scipy.sparse
        tree.check_subtree_valid()
        G = tree.current_subtree
        G_undirected = tree.current_subtree_undirected 
        basis = tree.basis
        basis_index = tree.basis_index

        # Let's go with a csr_array to store the (pairwise_distance, edge_list) matrix
        try:
            d = blist([])
            row_inds = blist([])
            col_inds = blist([])
        except NameError:
            d = []
            row_inds = []
            col_inds = []

        ########
        # for all pairs of edges, find which of the two connected nodes is the descendant of the other
        edge_2_descendant = {}
        edges = list(G.edges())
        num_edges = len(edges)
        for i, (v1, v2) in enumerate(edges):
            if i % 100 == 0:
                print(f"descendant iteration: {i}/{num_edges}")
            desc = get_descendant(G, v1, v2)
            edge_2_descendant[(v1, v2)] = desc
            edge_2_descendant[(v2, v1)] = desc

        # iterate over the pairs of nodes for which we have pairwise distances
        # In parallel, get all shortest paths
        print("Getting all shortest paths...")
        num_processes = multiprocessing.cpu_count() // 2
        pool = multiprocessing.Pool(num_processes)
        paths_list = pool.imap(self.shortest_path_parallel, zip(pairwise_distances, repeat(G_undirected)), chunksize=max(1, len(pairwise_distances) // num_processes))
        # The results should be ordered the same as the pairwise_dist_KOs
        print("Done getting all shortest paths")

        row_ind = -1
        for i, node_i in enumerate(pairwise_distances):
            # Some KOs may not be in the subtree selected, so append a row of zeros for those (i.e. don't add to the data list).
            # That won't affect the least squares fit
            if node_i in G_undirected:
                paths = next(paths_list)
            else:
                paths = dict()
            for node_j in pairwise_distances:
                row_ind += 1
                # if the nodes are the same, skip since this row of the matrix is all zeros
                if node_i != node_j:
                    # get the shortest path between the two nodes
                    try:
                        path = paths[node_i][node_j]
                        # get the index of each path element in the basis
                        path_tuples = [(path[i], path[i+1]) for i in range(len(path)-1)]
                        # for each tuple, find which is the descendant
                        path_edges = [edge_2_descendant[(x[0], x[1])] for x in path_tuples]
                        path_indices = [basis_index[node] for node in path_edges]
                        # set the corresponding entries in the sparse matrix to 1
                        for path_index in path_indices:
                            d.append(1)
                            row_inds.append(row_ind)
                            col_inds.append(path_index)
                    except KeyError:
                        # if there is no path between the two nodes, skip
                        pass
                if row_ind % 100000 == 0:
                    print(f"Finished {row_ind}/{len(pairwise_distances)**2} rows")
        A = sparse.csr_matrix((d, (row_inds, col_inds)), shape=(len(pairwise_distances)**2, len(basis)))
        return basis, A


    def least_square_parallel(self, args):
        def map_func(itr, A, y, factor, reg_factor):
            num_rows = int(factor * A.shape[1])
            row_indices = np.random.choice(A.shape[0], num_rows, replace=False)
            A_small = A[row_indices, :]
            y_small = y[row_indices]
            # append a row of 1's to A_small
            A_small = sparse.vstack([reg_factor * A_small, sparse.csr_matrix(np.ones(A_small.shape[1]))])
            # append a 0 to y_small
            y_small = np.append(reg_factor * y_small, 0)
            # Use lsq_linear to solve the NNLS problem
            res = lsq_linear(A_small, y_small, bounds=(0, 1), verbose=2)
            x = res.x
            return x
        return map_func(*args)


    def compute_edges(self, A, basis, edge_list, pairwise_distances: PairwiseDistance,
                      num_iter, factor, reg_factor, isdistance):
        # create the y vector of all pairwise distances
        y = pairwise_distances.get_pairwise_vector(isdistance)
        
        if num_iter < 1:
            raise ValueError('Number of iterations must be at least 1')
        if factor < 1:
            raise ValueError('Factor must be at least 1')
        
        if A.shape[0] != y.shape[0]:
            raise ValueError(f"The A matrix has {A.shape[0]} rows, but the y vector has {y.shape[0]} elements. "
                            f"Something is wrong.")
        
        num_threads = 1  # numpy apparently uses all PHYSICAL cores, twice that for hyperthreading
        pool = multiprocessing.Pool(num_threads)
        xs = np.array(pool.map(self.least_square_parallel, zip(range(num_iter), repeat(A), repeat(y), repeat(factor), repeat(reg_factor)), chunksize=num_iter // num_threads))

        # take the average of the solutions
        x = np.mean(xs, axis=0)

        df = edge_list
        # add a new column for the edge lengths
        df['edge_length'] = np.nan
        # iterate through the basis and add the edge lengths to the dataframe
        for i, tail in enumerate(basis):
            df.loc[df['child'] == tail, 'edge_length'] = x[i]
        # remove all the rows that aren't in the basis, so only focusing on the subtree defined by the given brite
        df = df[df['child'].isin(basis)]
        return df
