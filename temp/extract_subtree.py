# Open input file and read edge list
with open('../data/kegg_ko_edge_df_br_ko00001.txt_AAI_lengths_n_50_f_10_r_100.txt', 'r') as f:
    edges = [tuple(map(str, line.split())) for line in f.readlines()]

edges.remove(('#parent', 'child', 'edge_length'))

# Define a function to extract a subtree of the specified size
def extract_subtree_of_size(size, edges):
    def dfs(node, visited):
        # Returns the size of the subtree rooted at the given node, and the edges in the subtree
        size = 1
        subtree = []
        visited.add(node)
        for edge in edges:
            if edge[0] == node and edge[1] not in visited:
                child_size, child_subtree = dfs(edge[1], visited)
                size += child_size
                if size == target_size:
                    return size, subtree + [edge] + child_subtree
                elif size > target_size:
                    return size - child_size, subtree
                subtree += [edge] + child_subtree
        return size, subtree

    # Compute the target size of the output subtree
    target_size = size

    # Try extracting a subtree from each leaf node of the input tree
    best_size = 0
    best_subtree = None
    for leaf in range(1, len(edges) + 2):
        visited = set()
        size, subtree = dfs(leaf, visited)
        if size == target_size and (best_subtree is None or size > best_size):
            best_size = size
            best_subtree = subtree

    return best_subtree

# Extract a subtree of the specified size
subtree = extract_subtree_of_size(30, edges)
print(subtree)

# Write the subtree to an output file in edge list format
with open('./sub_graph_size_30_chatgpt.txt', 'w+') as f:
    for edge in subtree:
        f.write(f'{edge[0]} {edge[1]}\n')
