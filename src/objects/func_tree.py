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
    
    def make_subtree(self, classification):
        G = self._G
        # get the descendants of the brite
        descendants = self._get_descendants(G, classification)
        # add the root node back in (this will only have affect if you are using multiple brites. TODO: not yet implemented)
        descendants.add('root')
        # select the subgraph from the brite to the leaves
        subgraph = G.subgraph(descendants)
        self._subtree = subgraph
        return subgraph

    def _get_descendants(self, G, node):
        """
        Return all descendants of a node, including the node itself.

        :param G: networkx graph
        :param node: name of a node
        :return: set of nodes
        """
        descendants = set()
        descendants.add(node)
        for n in nx.descendants(G, node):
            descendants.add(n)
        return descendants

