import numpy as np
from typing import NewType, Dict

########################################################################
# 
########################################################################
UnifracDistance = NewType('UnifracDistance', float)
Flow = NewType('Flow', Dict[(int,int), float])
DifferentialAbundance = NewType('DifferentialAbundance', np.ndarray)
