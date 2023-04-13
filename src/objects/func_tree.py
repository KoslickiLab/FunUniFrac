import networkx as nx


class FuncTree:
    def __init__(self, G) -> None:
        self._G :nx.DiGraph = G
        self._subtree = None

    @property
    def main_tree(self)->nx.DiGraph:
        return self._G
        
    @property
    def current_subtree(self)->nx.DiGraph:
        self.check_subtree_valid()
        return self._subtree
    
    @property
    def basis(self):
        self.check_subtree_valid()
        # set the basis for the tree, which is an ordering of the edges. I'll identify each edge by its terminal node
        basis = [x for x in self._subtree.nodes()]
        return basis
    
    @property
    def basis_index(self):
        self.check_subtree_valid()
        basis_index = {node: i for i, node in enumerate(self.basis)}
        return basis_index
    
    def check_subtree_valid(self):
        if self._subtree is None:
            raise Exception("subtree has not been made")
    
    def set_subtree(self, classification):
        G = self._G
        # get the descendants of the brite
        descendants = FuncTree.get_descendants(G, classification)
        # add the root node back in (this will only have affect if you are using multiple brites. TODO: not yet implemented)
        descendants.add('root')
        # select the subgraph from the brite to the leaves
        subgraph = G.subgraph(descendants)
        self._subtree = subgraph
        return subgraph

    @classmethod
    def get_descendants(cls, G, node, leaf_only = False):
        """
        Return all descendants of a node, including the node itself.

        :param G: networkx graph
        :param node: name of a node
        :return: set of nodes
        """
            
        descendants = set()
        if leaf_only:
            for n in nx.descendants(G, node):
                if G.out_degree(n) == 0:
                    descendants.add(n)
        else:
            descendants.add(node)
            for n in nx.descendants(G, node):
                descendants.add(n)
        return descendants
    
    @classmethod
    def get_descendant(cls, G: nx.DiGraph, v1, v2):
        """
        of the two nodes v1 and v2, ASSUMED TO BE ADJACENT, find out which one is the descendant of the other

        :param graph: networkx graph, directed
        :param v1: node name
        :param v2: node name
        :return: descendant node name
        """
        if v1 in G.predecessors(v2):
            return v2
        elif v2 in G.predecessors(v1):
            return v1
        else:
            raise ValueError(f"Nodes {v1} and {v2} are not adjacent")

    @classmethod
    def get_root_of_tree(cls, G:nx.DiGraph):
        """
        Returns the root node of a directed tree

        :param G: directed tree
        :return: root node name
        """
        roots = [n for n, d in G.in_degree() if d == 0]
        if len(roots) > 1:
            raise Exception(f"The graph has multiple roots: {roots}")
        return roots[0]


    @classmethod
    def infer_edge_len_property(cls, G:nx.DiGraph):
        """
        From a networkx graph, find out the name of the edge length property

        :param G: graph
        :return: name of the edge length property
        """
        edge_properties = list(list(G.edges(data=True))[0][-1].keys())
        if len(edge_properties) > 1:
            raise ValueError(f'I found multiple edge properties: {edge_properties}. I don\'t know which one to use for '
                            f'edge lengths.')
        else:
            edge_len_property = edge_properties[0]
        return edge_len_property

    
    @classmethod
    def get_leaf_nodes(cls, G:nx.DiGraph):
        """
        Return all leaf nodes of a graph

        :param G: graph
        :return: set of leaf nodes
        """
        leaf_nodes = set()
        for n in G.nodes():
            if G.out_degree(n) == 0:
                leaf_nodes.add(n)
        leaf_nodes = sorted(list(leaf_nodes))
        return leaf_nodes

