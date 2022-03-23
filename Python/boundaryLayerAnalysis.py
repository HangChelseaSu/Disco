import os
import sys
import math
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.geom as geom
import discopy.plot as plot
import argparse as ag
from multiprocessing import Pool

try:
  import cmocean as cmo
  vrcmap = plt.get_cmap('cmo.delta')
  vortcmap = plt.get_cmap('cmo.balance_r')
except ModuleNotFoundError:
  vrcmap = plt.get_cmap('PRGn')
  vortcmap = plt.get_cmap('RdBu')


def analyzeSingle(filename):

    opts = util.loadOpts(filename)
    pars = util.loadPars(filename)
    print("Loading " + filename)
    t, x1, x2, x3, prim, dat = util.loadCheckpoint(filename)
    name = filename.split('/')[-1].split('.')[0].split('_')[-1]

    R = np.unique(x1)

    rjph = dat[0]
    piph = dat[3]

    nx = pars['Num_R']
    gam = pars['Adiabatic_Index']

    rho = prim[:, 0]
    P = prim[:, 1]
    v1 = prim[:, 2]
    v2 = prim[:, 3]
    v3 = prim[:, 4]
    q = prim[:, 5]

    x, y, z = geom.getXYZ(x1, x2, x3, opts, pars)
    vx, vy, vz = geom.getVXYZ(x1, x2, x3, v1, v2, v3, opts)

    print("Plotting") 
    dphi = 2.0*(piph-x2)

    fig, ax = plt.subplots(1, 1)

    plot.plotZSlice(fig, ax, rjph, piph, x1, rho, x3.mean(), r"$\rho$",
                    pars, opts, vmax = 1.2)
    print("Saving") 
    #fig.savefig('bl_deltaRho.png')
    plotname = "plot_BL_{0:s}_{1:s}.png".format(name, 'rho')
    fig.savefig(plotname, dpi = 400)
    plt.close(fig)

    print("Plotting vr")
    fig, ax = plt.subplots(1, 1)

    vrp = v1 * x1 * np.sqrt(rho)

    v1ext = np.max([-1*vrp.min(), vrp.max()])
    slt = 0.005 * v1ext

    plot.plotZSlice(fig, ax, rjph, piph, x1, vrp, x3.mean(), r"$r\sqrt{\Sigma}v_r$",
                    pars, opts, symlog=True, symlthresh=slt, vmax = v1ext, vmin = -1*v1ext, cmap=vrcmap )

    print("Saving")
    plotname = "plot_BL_{0:s}_{1:s}.png".format(name, 'vr')
    fig.savefig(plotname, dpi = 1000)
    plt.close(fig)

    #zoom-in
    fig, ax = plt.subplots(1, 1)

    vrp = v1 * x1 * np.sqrt(rho)

    v1ext = np.max([-1*vrp.min(), vrp.max()])
    slt = 0.005 * v1ext

    plot.plotZSlice(fig, ax, rjph, piph, x1, vrp, x3.mean(), r"$r\sqrt{\Sigma}v_r$",
                    pars, opts, symlog=True, symlthresh=slt, vmax = v1ext, vmin = -1*v1ext, cmap=vrcmap, square=True, rmax=1.0+(1.0-np.min(x1)) )

    print("Saving zoom")
    plotname = "plot_BL_{0:s}_{1:s}.png".format(name, 'vr_zoom')
    fig.savefig(plotname, dpi = 1000)
    plt.close(fig)

    print("Calculating div")
    divV = geom.calculateCurlV(x1, x2, x3, v1, v2, v3, dat, opts, pars)

    interior = (x1 > pars['R_Min']) & (x1 < pars['R_Max'] - 0.1)
    minDiv = divV[interior].min()
    maxDiv = divV[interior].max()
    if np.abs(maxDiv) < np.abs(minDiv): maxDiv = -1.0*minDiv
    else : minDiv = -1 * maxDiv

    print("Plotting")
    fig, ax = plt.subplots(1, 1)
    plot.plotZSlice(fig, ax, rjph, piph, x1, divV, x3.mean(), r"$(\nabla\times v)_z$", pars, opts, symlog=True, symlthresh=0.005*maxDiv, vmin=minDiv, vmax=maxDiv, cmap=vortcmap)
    print("Saving")
    plotname = "plot_BL_{0:s}_{1:s}.png".format(name, 'vorticity')
    fig.savefig(plotname, dpi = 1000)
    plt.close(fig)

    #zoom-in
    fig, ax = plt.subplots(1, 1)
    plot.plotZSlice(fig, ax, rjph, piph, x1, divV, x3.mean(), r"$(\nabla\times v)_z$", pars, opts, symlog=True, symlthresh=0.005*maxDiv, vmin=minDiv, vmax=maxDiv, cmap=vortcmap, square=True, rmax=1.0 + (1.0-np.min(x1)))

    print("Saving zoom")
    plotname = "plot_BL_{0:s}_{1:s}.png".format(name, 'vorticity_zoom')
    fig.savefig(plotname, dpi = 1000)
    plt.close(fig)


    return t, nx

def analyze(filenames):

    N = len(filenames)

    t = np.empty(N)
    nx = np.empty(N)

    for i, f in enumerate(filenames):
        dat = analyzeSingle(f)
        t[i] = dat[0]
        nx[i] = dat[1]

def plot_routine(f):
  analyzeSingle(f)

def process_images(files, ncpu):
  with Pool(ncpu) as pool:
    pool.starmap(plot_routine, zip(files))


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need a checkpoint dude!")
        sys.exit()

    parser = ag.ArgumentParser(description="Create 2D plots of Disco variables. Parallelized version.")
    parser.add_argument('checkpoints', nargs='+',
                            help="Checkpoint (.h5) files to plot.")
    parser.add_argument('-n', '--ncpu', type=int, nargs='?', action='store', const= os.cpu_count() - 1, default=False,
                            help="Turns on parallel processing. If number of cores not given, then it will default to N-1 cores.")

    args = parser.parse_args()
    files = args.checkpoints
    ncpu = args.ncpu

    if ncpu < 2:
      analyze(files)
    else:
      process_images(files, ncpu)
