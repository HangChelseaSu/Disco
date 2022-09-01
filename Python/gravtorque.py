import sys
import math
import argparse as ap
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.plot as plot
import discopy.geom as geom


def analysisSingle(filename, flux=True, fourier=True, summary=True):

    print("Loading {0}".format(filename))
    t, r, phi, z, prim, dat = util.loadCheckpoint(filename)

    rjph = dat[0]
    zkph = dat[1]
    primPhi0 = dat[2]
    piph = dat[3]
    planetDat = dat[4]
    nphi = dat[7]
    pars = util.loadPars(filename)
    opts = util.loadOpts(filename)
    name = filename.split('/')[-1].split('.')[0].split('_')[-1]


    k = int(zkph.shape[0]/2-1)
    ind = (z>zkph[k]) * (z<zkph[k+1])

    r = r[ind]
    phi = phi[ind]
    z = z[ind]
    prim = prim[ind,:]
    piph = piph[ind]

    primPhi0 = primPhi0[k,:]
    nphi = nphi[k,:]

    npl = planetDat.shape[0]

    x = r * np.cos(phi)
    y = r * np.sin(phi)

    Mp = planetDat[0,0]
    rp = planetDat[0,3]
    pp = planetDat[0,4]
    eps = planetDat[0,5]
    xp = rp * math.cos(pp)
    yp = rp * math.sin(pp)
    dx = x-xp
    dy = y-yp
    #rp = np.sqrt(dx*dx + dy*dy + eps*eps)
    rp = np.sqrt(dx*dx + dy*dy)

    TG = Mp*(dy*xp - dx*yp)/(rp*rp*rp)


    for p in range(1, npl):
        Mp = planetDat[p,0]
        rp = planetDat[p,3]
        pp = planetDat[p,4]
        eps = planetDat[p,5]
        xp = rp * math.cos(pp)
        yp = rp * math.sin(pp)
        dx = x-xp
        dy = y-yp
        #rp = np.sqrt(dx*dx + dy*dy + eps*eps)
        rp = np.sqrt(dx*dx + dy*dy)

        TG = Mp*(dy*xp - dx*yp)/(rp*rp*rp)

    TG = TG*prim[:,0]
    cmax = np.max([np.max(TG), -1*np.min(TG)])
    slt = cmax*0.005

    fig, ax = plt.subplots(1,1, figsize=(12,9))

    plot.plotZSlice(fig, ax, rjph, piph, r, TG, z.mean(), r'$\Gamma_g/\Delta V$',
   	  pars, opts, vmin=-cmax, vmax=cmax, rmax=4.0, planets=planetDat, symlog=True, symlthresh=slt, cmap=plt.get_cmap('PRGn') )

    dV = geom.getDV(dat, opts, pars)
    summ = np.sum(TG*dV)
    plotname = "plot_eq_{0:s}_symlog_{1:s}.png".format(name, "Tgrav")
    #plt.savefig(plotname)
    plt.close(fig)
    return t, summ


def analysisSumm(summaries):

    print("Summary Time series")
    time = []
    tq = []
    for i in range(len(summaries[:])):
      data = summaries[i]
      time.append(data[0])
      tq.append(data[1])
    plt.plot(time, tq)
    plt.xlabel('time')
    plt.ylabel('Gravitational Torque')
    #plt.savefig('timeseries.png')
    #plt.close(fig)
    print(tq)

if __name__ == "__main__":

    parser = ap.ArgumentParser(description="Perform deaddisk analysis")
    parser.add_argument('checkpoints', nargs='+', 
                            help="Checkpoint (.h5) files to plot.")

    args = parser.parse_args()
    files = args.checkpoints

    print(files)

    summaries = []
    for f in files:
        summ = analysisSingle(f)
        if summ[1] is not None:
            summaries.append(summ)

    if len(summaries) > 0:
        analysisSumm(summaries)
