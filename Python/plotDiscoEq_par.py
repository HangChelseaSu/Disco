import os
import sys
import argparse as ag
import h5py as h5
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.plot as plot
import discopy.geom as geom
from itertools import repeat
from multiprocessing import Pool

try:
  import cmocean as cmo
  symcmap = plt.get_cmap('cmo.balance')
except ModuleNotFoundError:
  symcmap = plt.get_cmap('RdBu')

def plotCheckpoint(file, vars=None, logvars=None, noGhost=False, om=None,
                    bounds=None, rmax=None, planets=False, k=None, symlogvars=None, slt=None):
    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim, dat = util.loadCheckpoint(file)
    rjph = dat[0]
    zkph = dat[1]
    primPhi0 = dat[2]
    piph = dat[3]
    pars = util.loadPars(file)
    opts = util.loadOpts(file)

    if k is None:
        k = int(zkph.shape[0]/2-1)

    zind = (z>zkph[k]) * (z<zkph[k+1])
    r = r[zind]
    phi = phi[zind]
    prim = prim[zind,:]
    primPhi0 = primPhi0[k,:]
    piph = piph[zind]


    if planets:
        planetDat = dat[4]
    else:
        planetDat = None

    if om is not None:
        phi1 = phi - om*t
        piph1 = piph - om*t
        if planetDat is not None:
            planetDat[:,4] -= om*t
    else:
        phi1 = phi
        piph1 = piph


    varnames, vartex, num_c, num_n = util.getVarNames(file)
    #nq = prim.shape[1]
    nq = num_c + num_n

    #Zs = np.unique(z)
    #z_eq = Zs[len(Zs)/2]
    #eq_ind = (z==z_eq)
    Z = z[zind].mean()
    title = "DISCO t = {0:.1f}".format(t/(2.0*np.pi))
    name = file.split('/')[-1].split('.')[0].split('_')[-1]

    if vars is None:
        vars = range(nq)

    if logvars is None:
        logvars = []

    if symlogvars is None:
        symlogvars = []

    for q in range(nq):
        if q in vars or q in logvars or q in symlogvars:
            print("   Plotting...")

            if bounds is not None:
                vmin, vmax = bounds[q]
            elif q >= num_c:
                vmin, vmax = 0.0, 1.0
            else:
                vmin, vmax = None, None

            if q in vars:

                fig, ax = plt.subplots(1,1, figsize=(12,9))

                plot.plotZSlice(fig, ax, rjph, piph1, r, prim[:,q], Z, vartex[q],
                                pars, opts, vmin=vmin, vmax=vmax, rmax=rmax, 
                                planets=planetDat)
                fig.suptitle(title, fontsize=24)
                plotname = "plot_eq_{0:s}_lin_{1:s}.png".format(name, varnames[q])
                
                print("   Saving {0:s}...".format(plotname))
                fig.savefig(plotname, dpi=300)
                plt.close(fig)

            if q in logvars:
                fig, ax = plt.subplots(1,1, figsize=(12,9))

                plot.plotZSlice(fig, ax, rjph, piph1, r, prim[:,q], Z, vartex[q],
                                pars, opts, vmin=vmin, vmax=vmax, rmax=rmax, 
                                planets=planetDat, log=True)
                fig.suptitle(title, fontsize=24)
                plotname = "plot_eq_{0:s}_log_{1:s}.png".format(name, varnames[q])

                print("   Saving {0:s}...".format(plotname))
                fig.savefig(plotname, dpi=300)
                plt.close(fig)

            if q in symlogvars:
                fig, ax = plt.subplots(1,1, figsize=(12,9))

                plot.plotZSlice(fig, ax, rjph, piph1, r, prim[:,q], Z, vartex[q],
                                pars, opts, vmin=vmin, vmax=vmax, rmax=rmax, 
                                planets=planetDat, symlog=True, symlthresh=slt, cmap=symcmap )
                fig.suptitle(title, fontsize=24)
                plotname = "plot_eq_{0:s}_symlog_{1:s}.png".format(name, varnames[q])

                print("   Saving {0:s}...".format(plotname))
                fig.savefig(plotname, dpi=300)
                plt.close(fig)

    if 10 in vars or 10 in logvars or 10 in symlogvars:
            v1 = prim[:,2]
            v2 = prim[:,3]
            v3 = prim[:,4]

            curl = geom.calculateZCurlV(r, phi, z, v1, v2, v3, dat, opts, pars)
            fig, ax = plt.subplots(1,1, figsize=(12,9))
            curl = np.sqrt(curl*curl)
            vmax = np.max(curl)
            vmin = np.min(curl)
            slt = vmin * 0.005*(vmax - vmin)
            vmin = -1.0*np.max(curl)

            plot.plotZSlice(fig, ax, rjph, piph1, r, curl, Z, vartex[q],
                                pars, opts, vmin=vmin, vmax=vmax, rmax=rmax, 
                                planets=planetDat, symlog=True, symlthresh=slt, cmap=symcmap )
            fig.suptitle(title, fontsize=24)
            plotname = "plot_eq_{0:s}_symlog_curl.png".format(name, varnames[q])

            print("   Saving {0:s}...".format(plotname))
            fig.savefig(plotname, dpi=300)
            plt.close(fig)

    
def getBounds(use_bounds, names, files):

    num_q = len(names)

    if use_bounds is not None:
        if use_bounds is True:
            bounds = plot.calcBounds(files)
        else:
            if os.path.isfile(use_bounds):
                bounds = plot.readBoundsFile(use_bounds, num_q)
            else:
                bounds = plot.calcBounds(files)
                plot.writeBoundsFile(use_bounds, names, bounds)
    else:
        bounds = None

    return bounds

def process_images(files, Vars, logvars, bounds, om, rmax, noghost, planets, ncpu, symlogvars, slt):
    
    with Pool(ncpu) as pool:
        pool.starmap(plot_routine, zip(files, repeat(Vars), repeat(logvars), repeat(bounds), 
			repeat(om), repeat(rmax), repeat(noghost), repeat(planets), repeat(symlogvars), repeat(slt)))

def plot_routine(f, Vars, logvars, bounds, om, rmax, noghost, planets, symlogvars, slt):
	
	plotCheckpoint(f, vars=Vars, logvars=logvars, bounds=bounds, om=om, 
                rmax=rmax, noGhost=noghost, planets=planets, symlogvars=symlogvars, slt=slt)

if __name__ == "__main__":

    parser = ag.ArgumentParser(description="Create 2D plots of Disco variables. Parallelized version.")
    parser.add_argument('checkpoints', nargs='+', 
                            help="Checkpoint (.h5) files to plot.")
    parser.add_argument('-v', '--vars', nargs='+', type=int,
                            help="Variables to plot.")
    parser.add_argument('-l', '--logvars', nargs='+', type=int,
                            help="Variables to plot logscale.")
    parser.add_argument('-s', '--symlogvars', nargs='+', type=int,
                            help="Variables to plot symlogscale.")
    parser.add_argument('-p', '--planets', action='store_true',
                            help="Plot planets.")
    parser.add_argument('-b', '--bounds', nargs='?', const=True,
                            help="Use global max/min for bounds. Optional argument BOUNDS is a file. If it exists, it will be read for parameter bounds. If it does not exist the global max/min will be calculated and saved to the file.")
    parser.add_argument('-r', '--rmax', type=float, 
                            help="Set plot limits to RMAX.")
    parser.add_argument('-o', '--omega', type=float, 
                            help="Rotate frame at rate OMEGA.")
    parser.add_argument('-st', '--symlogthresh', type=float, 
                            help="symlog threshold")
    parser.add_argument('-n', '--ncpu', type=int, nargs='?', action='store', const= os.cpu_count() - 1, default=False,
                            help="Turns on parallel processing. If number of cores not given, then it will default to N-1 cores.")
    parser.add_argument('--noghost', action='store_true', 
                            help="Do not plot ghost zones.")

    args = parser.parse_args()

    vars = args.vars
    logvars = args.logvars
    slogvars = args.symlogvars
    om = args.omega
    slt = args.symlogthresh
    rmax = args.rmax
    use_bounds = args.bounds
    planets = args.planets
    noghost = args.noghost
    ncpu = args.ncpu

    if ncpu < 1:
      ncpu = 1

    files = args.checkpoints

    names, texnames, num_c, num_n = util.getVarNames(files[0])

    bounds = getBounds(use_bounds, names, files)
    process_images(files, vars, logvars, bounds, om, rmax, noghost, planets, ncpu, slogvars, slt)
