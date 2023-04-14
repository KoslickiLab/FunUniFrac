import numpy as np
from typing import NewType

########################################################################
# 
########################################################################
FuncProfileVector = NewType('FuncProfileVector', np.ndarray)
    
def get_L1(P: FuncProfileVector, Q: FuncProfileVector):
    Z = np.sum(np.abs(P-Q))
    return Z

def get_L1_diffab(P: FuncProfileVector, Q: FuncProfileVector):
    diffab = P-Q
    return diffab

def get_L2(P: FuncProfileVector, Q: FuncProfileVector):
    Z = np.sum(np.abs(P-Q)**2)
    return Z

def get_L2_diffab(P: FuncProfileVector, Q: FuncProfileVector):
    diffab = (P-Q)**2
    return diffab

