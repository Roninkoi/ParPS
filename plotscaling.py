#!/bin/env python3
# Usage: ./plotscaling.py <file1.log> <file2.log> ...

import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import sys

from scipy.optimize import curve_fit

# import my plotting code
from plot import *

def inv(x, *p):
    a, b, c = p
    return a*1/x**b+c

show = True # show plots?
dfpath = "scaling.log"
if len(sys.argv) > 1:
    dfpath = sys.argv[1]
    data = np.loadtxt(dfpath)
    data = data[np.lexsort((data[:, 1], data[:, 0]))].T
elif len(sys.argv) > 2:
    datas = []
    for i in range(1, len(sys.argv)):
        dfpath = sys.argv[i]
        data = np.loadtxt(dfpath)
        data = data[np.lexsort((data[:, 1], data[:, 0]))].T
        datas.append(data)
    data = np.mean(data, axis=0)
else:
    data = np.loadtxt(dfpath)
    data = data[np.lexsort((data[:, 1], data[:, 0]))].T

sizes = []
u, ui, uc = np.unique(data[0], axis=0, return_counts=True, return_index=True)
labels=[]
styles=[]

for i in range(len(u)):
    scaling = np.array([np.array(data[1, ui[i]:ui[i]+uc[i]]), np.array(data[2, ui[i]:ui[i]+uc[i]]), ]).T
    sizes.append(scaling)
    labels.append("N = %i" % u[i])
    styles.append(0)

for i in range(len(u)):
    coeff, _ = curve_fit(inv, sizes[i][:, 0], sizes[i][:, 1], p0=[1., 1., 0.])
    x = np.linspace(np.min(sizes[i][:, 0]), np.max(sizes[i][:, 0]), 1000)
    sizes.append(np.array([x, inv(x, *coeff)]).T)
    print("fit", u[i], coeff)
    labels.append("N = %i fit" % u[i])
    styles.append(4)

linescatter(sizes, titles=["Scaling of wall time with number of processors", "$N_p$", "$t$ (s)"], labels=labels, fpath="scaling.pdf", colornum=len(u), styles=styles, show=show, log=[False, True])

