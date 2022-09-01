import sys
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt
import discopy.util as util

def plotCheckpoint(file):
    
    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim, dat = util.loadCheckpoint(file)
    pars = util.loadPars(file)
    opts = util.loadOpts(file)

    print("   Plotting...")
    nq = prim.shape[1]
    nc = opts['NUM_C']

    if nc <= 5:
        fig, ax = plt.subplots(2,3,figsize=(14,9))
        util.plotAx(ax[0,0], r, prim[:,0], "log", "log", r"$r$", r"$\rho$", 
                    'k+')
        util.plotAx(ax[0,1], r, prim[:,1], "log", "log", r"$r$", r"$P$", 
                    'k+')
        util.plotAx(ax[1,0], r, prim[:,2], "log", "linear", r"$r$", r"$u_1$",
                    'k+')
        util.plotAx(ax[1,1], r, prim[:,3], "log", "linear", r"$r$", 
                    r"$u_2$",'k+')
        util.plotAx(ax[1,2], r, prim[:,4], "log", "linear", r"$r$", r"$u_3$", 
                    'k+')
        if nq > 5:
            util.plotAx(ax[0,2], r, prim[:,5], "log", "linear", r"$r$", 
                    r"$q$", 'k+')
    else:
        fig, ax = plt.subplots(3,3,figsize=(14,12))
        util.plotAx(ax[0,0], r, prim[:,0], "linear", "linear", r"$r$", r"$\rho$", 
                    'k+')
        util.plotAx(ax[0,1], r, prim[:,1], "linear", "linear", r"$r$", r"$P$", 
                    'k+')
        util.plotAx(ax[1,0], r, prim[:,2], "linear", "linear", r"$r$", r"$u_1$",
                    'k+')
        util.plotAx(ax[1,1], r, prim[:,3], "linear", "linear", r"$r$", 
                    r"$u_2$",'k+')
        util.plotAx(ax[1,2], r, prim[:,4], "linear", "linear", r"$r$", r"$u_3$", 
                    'k+')
        if nq > 8:
            util.plotAx(ax[0,2], r, prim[:,8], "linear", "linear", r"$r$", 
                    r"$q$", 'k+')
        util.plotAx(ax[2,0], r, prim[:,5], "linear", "linear", r"$r$", r"$B_1$",
                    'k+')
        util.plotAx(ax[2,1], r, prim[:,6], "linear", "linear", r"$r$", 
                    r"$B_2$",'k+')
        util.plotAx(ax[2,2], r, prim[:,7], "linear", "linear", r"$r$", r"$B_z$", 
                    'k+')

    title = "DISCO t = {0:.3g}".format(t)

    #fig.suptitle(title, fontsize=18)

    plt.tight_layout()

    name = file.split('/')[-1].split('.')[0].split('_')[-1]
    plotname = "plot_r_{0:s}.png".format(name)
    
    print("   Saving {0:s}...".format(plotname))
    fig.savefig(plotname)

    plt.close(fig)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Makes plots of Disco prims as a function of r.")
        print("usage: python plotDiscoR.py <checkpoint.h5 ...>")
        sys.exit()

    files = sys.argv[1:]
    for f in files:
        plotCheckpoint(f)
