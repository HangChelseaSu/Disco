import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.plot as plot

def plotDiff(f, file0):

    pars = util.loadPars(f)
    opts = util.loadOpts(f)

    print("Loading " + file0)
    t0, x10, x20, x30, prim0, dat0 = util.loadCheckpoint(file0)
    print("Loading " + f)
    t, x1, x2, x3, prim, dat = util.loadCheckpoint(f)

    primSlice0 = dat0[2]
    primSlice = dat[2]
    diff = primSlice - primSlice0

    diff_eq = prim - prim0

    names, texnames, num_c, num_n = util.getVarNames(f)
    
    name = '.'.join(f.split('/')[-1].split('.')[:-1])

    for q in range(diff.shape[2]):
        figname = "diff_phi0_{0:s}_{1:s}.png".format(name, names[q])
        fig, ax = plt.subplots(1,1, figsize=(12,9))
        plot.plotPhiSlice(fig, ax, dat[0], dat[1], diff[:,:,q], texnames[q],
                        pars, opts)
        print("Saving " + figname)
        fig.savefig(figname)
        plt.close(fig)
        
        figname = "diff_eq_{0:s}_{1:s}.png".format(name, names[q])
        fig, ax = plt.subplots(1,1, figsize=(12,9))
        plot.plotZSlice(fig, ax, dat[0], dat[3], x1, diff_eq[:,q], x3.mean(),
                        texnames[q], pars, opts)
        print("Saving " + figname)
        fig.savefig(figname)
        plt.close(fig)

    

if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("Need at least 2 checkpoints bub!")
        sys.exit()

    filenames = sys.argv[1:-1]
    file0 = sys.argv[-1]

    for f in filenames:
        plotDiff(f, file0)

