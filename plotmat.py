#!/bin/env python3
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import sys

# plot matrix a
def matplot(a, titles=["", "", ""], fpath="", show=True):
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    plt.rcParams.update({'font.size': 22})
    f = plt.figure(figsize=(10, 10))
 
    plt.matshow(a, f)
        
    plt.title(titles[0])
    plt.xlabel(titles[1])
    plt.ylabel(titles[2])
    plt.colorbar()
    
    if len(fpath) > 0:
        f.savefig(fpath, bbox_inches='tight') # save to file

    if show:
        plt.show()
        plt.close()

show = True # show plots?
dfpath1 = "out.dat"
ndata = len(sys.argv) - 1
if len(sys.argv) > 1:
    dfpath1 = sys.argv[1]
data1 = np.loadtxt(dfpath1)

if len(sys.argv) > 2:
    dfpath2 = sys.argv[2]
    data2 = np.loadtxt(dfpath2)
    
# plot matrix from file
if ndata > 1:
    matplot(data1, ["Matrix 1", "$x$", "$y$"], fpath="matrix1.pdf")
    matplot(data2, ["Matrix 2", "$x$", "$y$"], fpath="matrix2.pdf")
    matplot(np.abs(data1-data2), ["Difference", "$x$", "$y$"], fpath="diff.pdf")
else:
    matplot(data1, ["Matrix", "$x$", "$y$"], fpath="matrix.pdf")
