import os
import sys
import argparse as ag
import h5py as h5
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.plot as plot
from itertools import repeat
from multiprocessing import Pool

def plotCheckpoint(file, vars=None, logvars=None, noGhost=False, om=None,
                    bounds=None, rmax=None, planets=False, k=None):
    
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
    title = "DISCO t = {0:.1f}".format(t)
    name = file.split('/')[-1].split('.')[0].split('_')[-1]

    if vars is None:
        vars = range(nq)

    if logvars is None:
        logvars = []

    for q in range(nq):
        if q in vars or q in logvars:
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
                fig.savefig(plotname, dpi=200)
                plt.close(fig)

            if q in logvars:
                fig, ax = plt.subplots(1,1, figsize=(12,9))

                plot.plotZSlice(fig, ax, rjph, piph1, r, prim[:,q], Z, vartex[q],
                                pars, opts, vmin=vmin, vmax=vmax, rmax=rmax, 
                                planets=planetDat, log=True)
                fig.suptitle(title, fontsize=24)
                plotname = "plot_eq_{0:s}_log_{1:s}.png".format(name, varnames[q])

                print("   Saving {0:s}...".format(plotname))
                fig.savefig(plotname, dpi=200)
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

def process_images(files, Vars, logvars, bounds, om, rmax, noghost, planets, ncpu):
    
    with Pool(ncpu) as pool:
        pool.starmap(plot_routine, zip(files, repeat(Vars), repeat(logvars), repeat(bounds), 
			repeat(om), repeat(rmax), repeat(noghost), repeat(planets)))

def plot_routine(f, Vars, logvars, bounds, om, rmax, noghost, planets):
	
	plotCheckpoint(f, vars=Vars, logvars=logvars, bounds=bounds, om=om, 
                rmax=rmax, noGhost=noghost, planets=planets)

if __name__ == "__main__":

    parser = ag.ArgumentParser(description="Create 2D plots of Disco variables. Parallelized version.")
    parser.add_argument('checkpoints', nargs='+', 
                            help="Checkpoint (.h5) files to plot.")
    parser.add_argument('-v', '--vars', nargs='+', type=int,
                            help="Variables to plot.")
    parser.add_argument('-l', '--logvars', nargs='+', type=int,
                            help="Variables to plot logscale.")
    parser.add_argument('-p', '--planets', action='store_true',
                            help="Plot planets.")
    parser.add_argument('-b', '--bounds', nargs='?', const=True,
                            help="Use global max/min for bounds. Optional argument BOUNDS is a file. If it exists, it will be read for parameter bounds. If it does not exist the global max/min will be calculated and saved to the file.")
    parser.add_argument('-r', '--rmax', type=float, 
                            help="Set plot limits to RMAX.")
    parser.add_argument('-o', '--omega', type=float, 
                            help="Rotate frame at rate OMEGA.")
    parser.add_argument('-n', '--ncpu', type=int, nargs='?', action='store', const= os.cpu_count() - 1, default=False,
                            help="Turns on parallel processing. If number of cores not given, then it will default to N-1 cores.")
    parser.add_argument('--noghost', action='store_true', 
                            help="Do not plot ghost zones.")

    args = parser.parse_args()

    vars = args.vars
    logvars = args.logvars
    om = args.omega
    rmax = args.rmax
    use_bounds = args.bounds
    planets = args.planets
    noghost = args.noghost
    ncpu = args.ncpu
    
    if ncpu < 1:
        raise SystemExit("CPU only has one core! Turn off parallel processing flag.")
    
    files = args.checkpoints

    names, texnames, num_c, num_n = util.getVarNames(files[0])

    bounds = getBounds(use_bounds, names, files)
    
    if ncpu:
        process_images(files, vars, logvars, bounds, om, rmax, noghost, planets, ncpu)
    else:
        for f in files:
                plotCheckpoint(f, vars=vars, logvars=logvars, bounds=bounds, om=om, 
			rmax=rmax, noGhost=noghost, planets=planets)

