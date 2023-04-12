import numpy as np

class PairwiseDistance:
    def __init__(self, 
                 pairwise_dists,
                 pairwise_dist_labels,
                 pairwise_dist_indices) -> None:
        self.dists = pairwise_dists
        self.labels = pairwise_dist_labels
        self.indices = pairwise_dist_indices
        
    def get_pairwise_vector(self, isdistance=True):
        # create the y vector of all pairwise distances
        y = []
        for func1 in self.labels:
            for func2 in self.labels:
                y.append(self.dists[self.indices[func1], self.indices[func2]])
        y = np.array(y)
        # by default, the values are: 0 = most dissimilar, 1 = most similar, so to convert to a distance, we subtract from 1
        if not isdistance:
            y = 1 - y
        return y
        