import numpy as np


class FuncTreeLeafPairwiseDistances:
    def __init__(self, 
                 pairwise_dists,
                 pairwise_dist_labels,
                 pairwise_dist_indices) -> None:
        self.dists = pairwise_dists
        self.labels = pairwise_dist_labels
        self.indices = pairwise_dist_indices
        
    def filter_labels(self, basis):
        """
        Make pairwise dist filtering explicit. 
        Given the pairwise KO distance matrix, return a list of labels and a
        dictionary mapping labels to indices. Only include the labels in the basis if specified.
        """
        # remove those that are not in the basis
        labels = [x for x in self.labels if x in basis]
        # also remove them from the pairwise dist index
        indices = {x: self.indices[x] for x in self.indices}
        return labels, indices

    def get_pairwise_vector(self, basis=None, isdistance=True):
        if basis is None:
            labels, indices = self.labels, self.indices
        else:
            labels, indices = self.filter_labels(basis)
        # create the y vector of all pairwise distances
        y = []
        for func1 in labels:
            for func2 in labels:
                print(f"funcs: {func1}, {func2}")
                print(indices[func1], indices[func2])
                print(type(indices[func1]), type(indices[func2]))
                print(self.dists)
                print(self.dists[indices[func1], indices[func2]])
                y.append(self.dists[indices[func1], indices[func2]])
        y = np.array(y)
        # by default, the values are: 0 = most dissimilar, 1 = most similar, so to convert to a distance, we subtract from 1
        if not isdistance:
            y = 1 - y
        return y
        