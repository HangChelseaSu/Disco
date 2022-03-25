import numpy as np
import matplotlib as mpl
from . import util
from . import geom


def plotZSlice(fig, ax, rjph, piph, r, q, Z, label, pars, opts, vmin=None,
               vmax=None, noGhost=False, colorbar=True, xlabel=None,
               ylabel=None, log=False, rmax=None, planets=None, cmap=None, 
               symlog=False, symlthresh=None, square=None):


    phi_max = pars['Phi_Max']

    if cmap is None:
        try:
            cmap = mpl.cm.inferno
        except:
            cmap = mpl.cm.afmhot

    if vmin is None:
        vmin = q.min()
    if vmax is None:
        vmax = q.max()

    if symlthresh is None:
        symlthresh = np.percentile(np.abs(q), 0.5)

    Rs = np.unique(r)
    if noGhost:
        Rs = Rs[:-2]

    if rmax is None or rmax <= 0.0:
        rmax = np.max(Rs)
        lim_float = True
    else:
        lim_float = False

    if square is True and xlabel is None:
      xlabel = r'r'
    if square is True and ylabel is None:
      ylabel = r'$\phi$'

    if log:
        norm = mpl.colors.LogNorm(vmin, vmax)
    else:
        norm = mpl.colors.Normalize(vmin, vmax)
    if symlog:
        norm = mpl.colors.SymLogNorm(symlthresh, vmin = vmin, vmax = vmax, base=10)

    pimh_min = np.inf
    rfactor = np.sqrt(2.1)
    if square is True: 
      rfactor = 1.0


    for i, R in enumerate(Rs):
      if R < rmax*rfactor:
        ind = (r == R)
        imax = np.argmax(piph[ind])

        apiph = np.roll(piph[ind], -imax-1)
        aq = np.roll(q[ind], -imax-1)

        phif = np.empty(len(apiph)+1)
        phif[1:] = apiph[:]
        phif[0] = apiph[-1] - phi_max

        rf = rjph[i:i+2]

        # x = rf[:,None] * np.cos(phif)[None,:]
        # y = rf[:,None] * np.sin(phif)[None,:]

        if square is True:
          X, Y = np.meshgrid(rf, phif)

          C = ax.pcolormesh(X, Y, aq.reshape((len(aq),1)),
                cmap=cmap, vmin=vmin, vmax=vmax, norm=norm, shading='flat')
        else:
          x, y, z = geom.getXYZ(rf[:,None], phif[None,:], Z.mean(), opts, pars)

        C = ax.pcolormesh(x, y, aq[None,:], cmap=cmap, norm=norm,
                          vmin=vmin, vmax=vmax)

        if lim_float and rf.max() > rmax:
            rmax = rf.max()
        if phif.min() < pimh_min:
            pimh_min =  phif.min()

    if planets is not None:
        rpl = planets[:,3]
        ppl = planets[:,4]
        xpl = rpl * np.cos(ppl)
        ypl = rpl * np.sin(ppl)
        ax.plot(xpl, ypl, color='grey', ls='', marker='o', mew=0, ms=5)

    if square is True:
        ax.set_aspect('auto')
        ax.set_xlim(Rs.min(), rmax)
        ax.set_ylim(pimh_min, piph.max())
    else:
        ax.set_aspect('equal')
    
    if opts['GEOMETRY'] == 'cylindrical':
        ax.set_xlim(-rmax, rmax)
        ax.set_ylim(-rmax, rmax)
    elif opts['GEOMETRY'] == 'spherical':
        ax.set_xlim(-rmax, rmax)
        ax.set_ylim(-rmax, rmax)
    else:
        ax.set_xlim(rjph.min(), rjph.max())
        ax.set_ylim(pimh_min, piph.max())

    if colorbar:
        cb = fig.colorbar(C, ax=ax)
        cb.set_label(label, fontsize=24)

    if xlabel == None:
        xlabel = r'$x$'
    if ylabel == None:
        ylabel = r'$y$'

    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)


def plotPhiSlice(fig, ax, rjph, zkph, q, label, pars, opts, vmin=None,
                 vmax=None, noGhost=False, colorbar=True, xlabel=None,
                 ylabel=None, log=False, rmax=None, planets=None,
                 cmap=None, symlog=False, symlthresh=None):

    if cmap is None:
        try:
            cmap = mpl.cm.inferno
        except:
            cmap = mpl.cm.afmhot

    if vmin is None:
        vmin = q.min()
    if vmax is None:
        vmax = q.max()
    if symlthresh is None:
        symlthresh = np.percentile(np.abs(q), 0.5)

    Nr = rjph.shape[0]-1
    Nz = zkph.shape[0]-1

    if noGhost:
        Nr_loc = Nr
        Nz_loc = Nz
        j1 = 0
        k1 = 0

        if pars['NoBC_Rmin'] == 0:
            Nr_loc -= 2
            j1 = 2
        if pars['NoBC_Rmax'] == 0:
            Nr_loc -= 2
        if pars['NoBC_Zmin'] == 0:
            Nz_loc -= 2
            k1 = 2
        if pars['NoBC_Zmax'] == 0:
            Nz_loc -= 2

        dat = np.empty((Nz_loc, Nr_loc))
        for k in range(Nz_loc):
            dat[k, :] = q[k+k1, j1:j1+Nr_loc]
        rjph_loc = rjph[j1:j1+Nr_loc+1]
        zkph_loc = zkph[k1:k1+Nz_loc+1]
    else:
        dat = q
        rjph_loc = rjph
        zkph_loc = zkph

    if rmax is None or rmax <= 0.0:
        rmax = 0.0
        lim_float = True
    else:
        lim_float = False

    if log:
        norm = mpl.colors.LogNorm(vmin, vmax)
    else:
        norm = mpl.colors.Normalize(vmin, vmax)
    if symlog:
        norm = mpl.colors.SymLogNorm(symlthresh, vmin = vmin, vmax = vmax, base=10)

    x, y, z = geom.getXYZ(rjph_loc[None, :], 0.0, zkph_loc[:, None],
                          opts, pars)

    C = ax.pcolormesh(x, z, dat, cmap=cmap, norm=norm)

    if lim_float and rjph_loc.max() > rmax:
        rmax = rjph_loc.max()

    if planets is not None:
        rpl = planets[:, 3]
        ppl = planets[:, 4]
        xpl = rpl * np.cos(ppl)
        ypl = rpl * np.sin(ppl)
        ax.plot(xpl, ypl, color='grey', ls='', marker='o', mew=0, ms=5)

    ax.set_aspect('equal')
    # ax.set_xlim(rjph.min(), rmax)
    # ax.set_ylim(zkph.min(), zkph.max())
    ax.set_xlim(x.min(), x.max())
    ax.set_ylim(z.min(), z.max())

    if colorbar:
        cb = fig.colorbar(C)
        cb.set_label(label, fontsize=24)

    if xlabel is None:
        xlabel = r'$x$'
    if ylabel is None:
        ylabel = r'$z$'

    ax.set_xlabel(xlabel, fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)


def calcBounds(files):

    print("Checking bounds...")
    f = files[0]
    
    print("  Checking bounds for {0:s}...".format(f))
    prim = util.loadCheckpointPrims(f)

    num_q = prim.shape[1]
    bounds = [[prim[:,q].min(), prim[:,q].max()] for q in range(num_q)]
    bounds = np.array(bounds)

    for f in files[1:]:
        print("  Checking bounds for {0:s}...".format(f))
        prim = util.loadCheckpointPrims(f)
        for q in range(num_q):
            bounds[q,0] = min(bounds[q,0], prim[:,q].min())
            bounds[q,1] = max(bounds[q,1], prim[:,q].max())
    
    return bounds

def writeBoundsFile(filename, names, bounds):

    lines = ["{0:s} {1:.12g} {2:.12g}\n".format(names[q], bounds[q,0], bounds[q,1])
                for q in range(len(names))]

    f = open(filename, "w")
    for line in lines:
        f.write(line)

    f.close()

def readBoundsFile(filename, num_q):

    bounds = np.loadtxt(filename, usecols=[1,2])

    bounds = bounds[:num_q]

    return bounds
    
