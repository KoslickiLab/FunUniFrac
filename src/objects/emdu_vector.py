import numpy as np
from typing import NewType

ProfileVector = NewType('ProfileVector', np.ndarray)
    
def get_L1(P: ProfileVector, Q: ProfileVector):
    Z = np.sum(np.abs(P-Q))
    return Z

def get_L1_diffab(P: ProfileVector, Q: ProfileVector):
    diffab = P-Q
    return diffab

def get_L2(P: ProfileVector, Q: ProfileVector):
    Z = np.sum(np.abs(P-Q)**2)
    return Z

def get_L2_diffab(P: ProfileVector, Q: ProfileVector):
    diffab = (P-Q)**2
    return diffab

