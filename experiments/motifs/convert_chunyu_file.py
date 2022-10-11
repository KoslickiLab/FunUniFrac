import numpy as np
import pickle
res = pickle.load(open('final_jindex_results.pkl', 'rb'))
jac_dicts = {}
basis_set = set()
for x,y,d in res:
    jac_dicts[x,y] = d
    jac_dicts[y,x] = d
    basis_set.add(x)

basis = list(basis_set)
basis.sort()
with open('final_jindex_results_dists_pw_mat.labels.txt', 'w') as f:
    for x in basis:
        f.write(f"{x}\n")


dists = np.zeros((len(basis), len(basis)))
for i,x in enumerate(basis):
    for j,y in enumerate(basis):
        if i == j:
            dists[i,j] = 1
        else:
            dists[i,j] = jac_dicts[x,y]

np.save('final_jindex_results_dists_pw_mat.npy', dists)

