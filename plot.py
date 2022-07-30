# Roni Koitermaa 2022
# Plotting code using Matplotlib
# Line/scatter plots, subplots, histogram plot, bar plot

import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# line + scatter plot of multiple data sets a
# with different labels and styles
def linescatter(a, titles=["", "", ""], labels=[], styles=[], fpath="", show=True, log=[False, False], colornum=0):
    plt.style.use('default')
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    plt.rcParams.update({'font.size': 22})
    #colors = plt.cm.hsv(np.linspace(0.66, 0.0, len(a)))
    #colors = plt.cm.plasma(np.linspace(0, 0.95, len(a)))
    if colornum == 0:
        colornum = len(labels)
    if len(a) > 5:
        colors = plt.cm.plasma(np.linspace(0, 0.9, colornum))
        plt.rcParams['axes.prop_cycle'] = plt.cycler(color=colors)
    f = plt.figure(figsize=(10, 10))
    ax = plt.gca()
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.5g'))

    if log[0]:
        ax.set_xscale('log')
    if log[1]:
        ax.set_yscale('log')

    i = 0
    for ai in a:
        if len(ai) <= 0:
            continue
        
        alabel = str(i)
        amarker = ""
        alinestyle = "-"
        amarkersize = 1.0
        errfill = False
        
        if len(labels) > 0:
            if i < len(labels):
                alabel = labels[i]
            else:
                alabel = labels[-1]
            
        if len(styles) > 0:
            if styles[i] == 1:
                amarker = "o"
                alinestyle = ""
                amarkersize = 3.0
            if styles[i] == 4:
                alinestyle = "dashed"
            if styles[i] == 10:
                errfill = True

        if errfill and len(ai[0]) > 2:
            plt.plot(ai[:, 0], ai[:, 1], label=alabel, marker=amarker, markersize=amarkersize, linestyle=alinestyle)
            plt.fill_between(a[0][:, 0], a[0][:, 1]-a[0][:, 2], a[0][:, 1]+a[0][:, 2], alpha=0.33)
        elif len(ai[0]) > 2:
            plt.errorbar(ai[:, 0], ai[:, 1], yerr=ai[:, 2], label=alabel, marker=amarker, markersize=amarkersize, linestyle=alinestyle)
        elif len(ai[0]) > 3:
            plt.errorbar(ai[:, 0], ai[:, 1], yerr=ai[:, 3], xerr=ai[:, 2], label=alabel, marker=amarker, markersize=amarkersize, linestyle=alinestyle)
        else:
            plt.plot(ai[:, 0], ai[:, 1], label=alabel, marker=amarker, markersize=amarkersize, linestyle=alinestyle)
        i += 1
        
    plt.title(titles[0])
    plt.xlabel(titles[1])
    plt.ylabel(titles[2])
    
    plt.grid(True)
    if len(labels) > 0 or len(a) > 1:
        plt.legend(fontsize='xx-small')

    if len(fpath) > 0:
        f.savefig(fpath, bbox_inches='tight', dpi=300) # save to file

    if show:
        plt.show()
        plt.close()

# line plot with two subplots
def linesub2(a, titles=["", "", ""], labels=["", ""], fpath="", show=True):
    plt.style.use('default')
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    plt.rcParams.update({'font.size': 22})
    f, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(10, 10))

    ax1.plot(a[0][:, 0], a[0][:, 1], label=labels[0])
    ax2.plot(a[1][:, 0], a[1][:, 1], color="tab:orange", label=labels[1])
    
    if (len(a[0][0]) > 2):
        ax1.fill_between(a[0][:, 0], a[0][:, 1]-a[0][:, 2], a[0][:, 1]+a[0][:, 2], alpha=0.33, facecolor="tab:blue")
    if (len(a[1][0]) > 2):
        ax2.fill_between(a[1][:, 0], a[1][:, 1]-a[1][:, 2], a[1][:, 1]+a[1][:, 2], alpha=0.33, facecolor="tab:orange")

    ax1.set_title(titles[0])
    plt.xlabel(titles[1])
    ax1.set_ylabel(titles[2])
    if len(titles) > 3:
        ax2.set_ylabel(titles[3])
    else:
        ax2.set_ylabel(titles[2])
        
    ax1.grid(True)
    ax2.grid(True)
    if len(labels[0]) > 0:
        ax1.legend(fontsize='xx-small')
    if len(labels[1]) > 0:
        ax2.legend(fontsize='xx-small')
        
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.5g'))
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.5g'))

    if len(fpath) > 0:
        f.savefig(fpath, bbox_inches='tight', dpi=300)

    if show:
        plt.show()
        plt.close()

# line plot with three subplots
def linesub3(a, titles=["", "", ""], labels=["", "", ""], fpath="", show=True):
    plt.style.use('default')
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    plt.rcParams.update({'font.size': 22})
    f, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, figsize=(10, 10))
    
    ax1.plot(a[0][:, 0], a[0][:, 1], color="tab:blue", label=labels[0])
    ax2.plot(a[1][:, 0], a[1][:, 1], color="tab:orange", label=labels[1])
    ax3.plot(a[2][:, 0], a[2][:, 1], color="tab:green", label=labels[2])

    if (len(a[0][0]) > 2):
        ax1.fill_between(a[0][:, 0], a[0][:, 1]-a[0][:, 2], a[0][:, 1]+a[0][:, 2], alpha=0.33, facecolor="tab:blue")
    if (len(a[1][0]) > 2):
        ax2.fill_between(a[1][:, 0], a[1][:, 1]-a[1][:, 2], a[1][:, 1]+a[1][:, 2], alpha=0.33, facecolor="tab:orange")
    if (len(a[2][0]) > 2):
        ax3.fill_between(a[2][:, 0], a[2][:, 1]-a[2][:, 2], a[2][:, 1]+a[2][:, 2], alpha=0.33, facecolor="tab:green")
    
    ax1.set_title(titles[0])
    plt.xlabel(titles[1])
    ax1.set_ylabel(titles[2])
    if len(titles) > 4:
        ax2.set_ylabel(titles[3])
        ax3.set_ylabel(titles[4])
    else:
        ax2.set_ylabel(titles[2])
        ax3.set_ylabel(titles[2])
        
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    if len(labels[0]) > 0:
        ax1.legend(fontsize='xx-small')
    if len(labels[1]) > 0:
        ax2.legend(fontsize='xx-small')
    if len(labels[2]) > 0:
        ax3.legend(fontsize='xx-small')
    
    ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.5g'))
    ax2.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.5g'))
    ax3.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.5g'))

    if len(fpath) > 0:
        f.savefig(fpath, bbox_inches='tight', dpi=300)

    if show:
        plt.show()
        plt.close()
        
# histogram plot
def histplt(x, nbins, theor=None, titles=["Histogram", "x", "N"], label="", fpath="", show=True):
    plt.style.use('default')
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    plt.rcParams.update({'font.size': 22})
    f = plt.figure(figsize=(10, 10))
    # generate (normalized) histogram of data x
    plt.hist(x, nbins, ec="tab:blue", label=label) #, density=True
    if theor is not None:
        tx = np.linspace(np.min(x), np.max(x), 1000)
        plt.plot(tx, theor(tx), label="Theory")
        plt.legend(fontsize='xx-small')

    plt.title(titles[0])
    plt.xlabel(titles[1])
    plt.ylabel(titles[2])
    
    #plt.grid(True)

    if len(fpath) > 0:
        f.savefig(fpath, bbox_inches='tight', dpi=300) # save to file

    if show:
        plt.show()
        plt.close()
        
# bar plot
def barplt(a, fit, titles=["", "", ""], labels="", al=1.0, fpath="", show=True):
    plt.style.use('default')
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    plt.rcParams.update({'font.size': 22})
    f = plt.figure(figsize=(10, 10))
    
    i = 0
    for ai in a:
        if len(ai) <= 0:
            continue
        
        alabel = str(i)
        
        if len(labels) > 0:
            alabel = labels[i]
            
        if len(fit) > 0:
            fi = fit[i]
            plt.plot(fi[:, 0], fi[:, 1], label=alabel+" fit")

        eps=0.0
        plt.bar(ai[:, 0]-eps, ai[:, 2], (ai[:, 1]-ai[:,0])+2*eps, label=alabel, alpha=al, align="edge")
        i += 1

    plt.title(titles[0])
    plt.xlabel(titles[1])
    plt.ylabel(titles[2])

    #plt.grid(True)

    if len(labels) > 0 or len(a) > 1:
        plt.legend(fontsize='xx-small')

    if len(fpath) > 0:
        f.savefig(fpath, bbox_inches='tight', dpi=300) # save to file

    if show:
        plt.show()
        plt.close()
        
# plot matrix a
def matplot(a, titles=["", "", ""], fpath="", show=True):
    plt.style.use('default')
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    plt.rcParams.update({'font.size': 22})
    f = plt.figure(figsize=(10, 10))
    #ax.set_aspect('equal')
    plt.matshow(a, f)
        
    plt.title(titles[0])
    plt.xlabel(titles[1])
    plt.ylabel(titles[2])
    plt.colorbar()
    
    if len(fpath) > 0:
        f.savefig(fpath, bbox_inches='tight', dpi=300) # save to file

    if show:
        plt.show()
        plt.close()

