#!/usr/bin/env python
import networkx as nx
from scipy import sparse
import multiprocessing
try:
    from blist import blist
except ModuleNotFoundError:
    print("Warning: Could not import blist. Please install blist to speed up the path matrix calculation.")
from src.algorithms.lp_edge_length import get_descendant
from src.utils.func_tree import FuncTree
from itertools import repeat


class EdgeLengthSolver:
    def __init__(self) -> None:
        pass

    def map_star(self, args):
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
        paths_list = pool.imap(self.map_star, zip(pairwise_distances, repeat(G_undirected)), chunksize=max(1, len(pairwise_distances) // num_processes))
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




# def main(args):
#         edge_list = args.edge_list
#         out_dir = args.out_dir
#         distances_file = args.distances
#         brite = args.brite_id
#         # check that the files exist
#         # if not exists(edge_list):
#         #     raise FileNotFoundError(f"Could not find {edge_list}")
#         # if not exists(distances_file):
#         #     raise FileNotFoundError(f"Could not find {distances_file}")
#         # if not exists(f"{distances_file}.labels.txt"):
#         #     raise FileNotFoundError(f"Could not find {distances_file}.labels.txt in the same directory as {distances_file}")
#         edge_list = data.get_data_abspath(edge_list)
#         distances_file = data.get_data_abspath(distances_file)
#         distances_labels_file = data.get_data_abspath(f"{distances_file}.labels.txt")
#         if not exists(out_dir):
#             os.mkdir(out_dir)

#         # check if brite is legit
#         if brite not in kegg_db.instance.brites:
#             raise ValueError(f"{brite} is not a valid BRITE ID. Choices are: {kegg_db.instance.brites}")

#         matrix_name = f"{brite}_{os.path.basename(distances_file)}_A.npz"
#         basis_name = f"{brite}_{os.path.basename(distances_file)}_column_basis.txt"

#         # import pairwise distances
#         pairwise_dist = np.load(distances_file)
#         # import label names
#         pairwise_dist_KOs, pairwise_dist_KO_index = get_KO_labels_and_index(distances_labels_file, basis=basis)

#         ########################################################################################################################
#         # Let's do the following: since I've already computed all pairwise distances, we can just make a large
#         # least squares problem fitting the tree distances to the pairwise distances
#         # Let's get the matrix describing which edges are traversed between all pairs of nodes
#         # This is a sparse matrix, so we'll need to use scipy.sparse

#         # read in the edge list
#         G = nx.read_edgelist(edge_list, delimiter='\t', nodetype=str, create_using=nx.DiGraph)

#         # Save the basis
#         with open(f"{os.path.join(out_dir, basis_name)}", 'w') as f:
#             f.write('\n'.join(basis))

#         # save the sparse matrix
#         print(f"Saving sparse matrix to {os.path.join(out_dir, matrix_name)}")
#         sparse.save_npz(os.path.join(out_dir, matrix_name), A)
