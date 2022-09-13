import os
import argparse as ag
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.plot as plot


def plotCheckpoint(file, vars=None, logvars=None, noGhost=False, bounds=None,
                   rmax=None, planets=False, dpi=300):

    print("Loading {0:s}...".format(file))

    t, r, phi, z, prim, dat = util.loadCheckpoint(file)
    rjph = dat[0]
    zkph = dat[1]
    primPhi0 = dat[2]
    pars = util.loadPars(file)
    opts = util.loadOpts(file)

    if planets:
        planetDat = dat[4]
    else:
        planetDat = None

    varnames, vartex, num_c, num_n = util.getVarNames(file)
    nq = num_c + num_n

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

                fig, ax = plt.subplots(1, 1, figsize=(9, 12))

                plot.plotPhiSlice(fig, ax, rjph, zkph, primPhi0[:, :, q],
                                  vartex[q], pars, opts,
                                  vmin=vmin, vmax=vmax, rmax=rmax,
                                  planets=planetDat)
                fig.suptitle(title, fontsize=24)
                plotname = "plot_phi0_{0:s}_lin_{1:s}.png".format(
                    name, varnames[q])

                print("   Saving {0:s}...".format(plotname))
                fig.savefig(plotname, dpi=dpi)
                plt.close(fig)

            if q in logvars:
                fig, ax = plt.subplots(1, 1, figsize=(9, 12))

                plot.plotPhiSlice(fig, ax, rjph, zkph, primPhi0[:, :, q],
                                  vartex[q], pars, opts,
                                  vmin=vmin, vmax=vmax, rmax=rmax,
                                  planets=planetDat, log=True)
                fig.suptitle(title, fontsize=24)
                plotname = "plot_phi0_{0:s}_log_{1:s}.png".format(
                    name, varnames[q])

                print("   Saving {0:s}...".format(plotname))
                fig.savefig(plotname, dpi=dpi)
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


if __name__ == "__main__":

    parser = ag.ArgumentParser(
        description="Create 2D plots of Disco variables along phi = 0.")
    parser.add_argument('checkpoints', nargs='+',
                        help="Checkpoint (.h5) files to plot.")
    parser.add_argument('-v', '--vars', nargs='+', type=int,
                        help="Variables to plot.")
    parser.add_argument('-l', '--logvars', nargs='+', type=int,
                        help="Variables to plot logscale.")
    parser.add_argument('-p', '--planets', action='store_true',
                        help="Plot planets.")
    parser.add_argument('-b', '--bounds', nargs='?', const=True,
                        help=("Use global max/min for bounds. Optional "
                              "argument BOUNDS is a file. If it exists, it "
                              "will be read for parameter bounds. If it does "
                              "not exist the global max/min will be "
                              "calculated and saved to the file."))
    parser.add_argument('-r', '--rmax', type=float,
                        help="Set plot limits to RMAX.")
    parser.add_argument('--noghost', action='store_true',
                        help="Do not plot ghost zones.")
    parser.add_argument('-d', '--dpi', type=int, default=300,
                            help="DPI of final image")

    args = parser.parse_args()

    vars = args.vars
    logvars = args.logvars
    rmax = args.rmax
    use_bounds = args.bounds
    planets = args.planets
    noghost = args.noghost
    dpi = args.dpi

    files = args.checkpoints

    names, texnames, num_c, num_n = util.getVarNames(files[0])

    bounds = getBounds(use_bounds, names, files)

    for f in files:
        plotCheckpoint(f, vars=vars, logvars=logvars, bounds=bounds,
                       rmax=rmax, noGhost=noghost, planets=planets,
                       dpi=dpi)
