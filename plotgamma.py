#!/bin/env python3
# Usage: ./plotgamma.py <gamma.log>

import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import sys

# import my plotting code
from plot import *

show = True # show plots?
dfpath = "gamma.log"
if len(sys.argv) > 1:
    dfpath = sys.argv[1]
    
data = np.loadtxt(dfpath)
data = data[np.lexsort((data[:, 1], data[:, 0]))].T

sizes = []
u, ui, uc = np.unique(data[0], axis=0, return_counts=True, return_index=True)
labels=[]

for i in range(len(u)):
    gamma = np.array([np.array(data[1, ui[i]:ui[i]+uc[i]]), np.array(data[2, ui[i]:ui[i]+uc[i]])]).T
    sizes.append(gamma)
    labels.append("$N$ = %i" % u[i])

linescatter(sizes, titles=["Relaxation parameter in different system sizes", "$\gamma$", "Iteration number"], labels=labels, fpath="gamma.pdf", show=show, log=[False, True])

