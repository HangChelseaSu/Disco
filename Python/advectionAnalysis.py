import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.geom as geom


def calcAdvection(t, x, y, z, pars, tol=1.0e-6):

    rho = np.empty(x.shape)
    P = np.empty(x.shape)
    vx = np.empty(x.shape)
    vy = np.empty(x.shape)
    vz = np.empty(x.shape)

    vx0 = pars['Init_Par1']
    vy0 = pars['Init_Par2']
    vz0 = pars['Init_Par3']
    gam = pars['Adiabatic_Index']

    rho[:] = 1.0
    P[:] = 1.0/gam
    vx[:] = vx0
    vy[:] = vy0
    vz[:] = vz0

    return rho, P, vx, vy, vz


def analyzeSingle(filename):

    opts = util.loadOpts(filename)
    pars = util.loadPars(filename)
    print("Loading " + filename)
    t, x1, x2, x3, prim, dat = util.loadCheckpoint(filename)

    nx = pars['Num_R']
    gam = pars['Adiabatic_Index']

    rho = prim[:, 0]
    P = prim[:, 1]
    v1 = prim[:, 2]
    v2 = prim[:, 3]
    v3 = prim[:, 4]
    eS = P * np.power(rho, -gam)

    x, y, z = geom.getXYZ(x1, x2, x3, opts, pars)
    vx, vy, vz = geom.getVXYZ(x1, x2, x3, v1, v2, v3, opts)

    r = np.sqrt(x*x + y*y)
    cp = x / r
    sp = y / r
    vr = cp*vx + sp*vy
    om = (-sp*vx + cp*vy) / r

    gam = pars['Adiabatic_Index']

    rhoS, PS, vxS, vyS, vzS = calcAdvection(t, x, y, z, pars)
    eSS = PS * np.power(rhoS, -gam)
    vrS = cp*vxS + sp*vyS
    omS = (-sp*vxS + cp*vyS) / r

    dV = geom.getDV(dat, opts, pars)

    errRho = geom.integrate(np.fabs(rho-rhoS), dat, opts, pars, dV)
    errP = geom.integrate(np.fabs(P-PS), dat, opts, pars, dV)
    errVx = geom.integrate(np.fabs(vx-vxS), dat, opts, pars, dV)
    errVy = geom.integrate(np.fabs(vy-vyS), dat, opts, pars, dV)
    errVz = geom.integrate(np.fabs(vz-vzS), dat, opts, pars, dV)
    errVr = geom.integrate(np.fabs(vr-vrS), dat, opts, pars, dV)
    errOm = geom.integrate(np.fabs(om-omS), dat, opts, pars, dV)
    errS = geom.integrate(np.fabs(eS-eSS), dat, opts, pars, dV)

    return t, nx, errRho, errP, errVx, errVy, errVz, errVr, errOm, errS


def analyze(filenames):

    N = len(filenames)

    t = np.empty(N)
    nx = np.empty(N)
    errRho = np.empty(N)
    errP = np.empty(N)
    errVx = np.empty(N)
    errVy = np.empty(N)
    errVz = np.empty(N)
    errVr = np.empty(N)
    errOm = np.empty(N)
    errS = np.empty(N)

    for i, f in enumerate(filenames):
        dat = analyzeSingle(f)
        t[i] = dat[0]
        nx[i] = dat[1]
        errRho[i] = dat[2]
        errP[i] = dat[3]
        errVx[i] = dat[4]
        errVy[i] = dat[5]
        errVz[i] = dat[6]
        errVr[i] = dat[7]
        errOm[i] = dat[8]
        errS[i] = dat[9]

    makeErrPlot(t, nx, errRho, "advection_errRho", r"$L_1(\rho)$")
    makeErrPlot(t, nx, errP, "advection_errP", r"$L_1(P)$")
    makeErrPlot(t, nx, errVx, "advection_errVx", r"$L_1(v_x)$")
    makeErrPlot(t, nx, errVy, "advection_errVy", r"$L_1(v_y)$")
    makeErrPlot(t, nx, errVz, "advection_errVz", r"$L_1(v_z)$")
    makeErrPlot(t, nx, errVr, "advection_errVr", r"$L_1(v_r)$")
    makeErrPlot(t, nx, errOm, "advection_errOm", r"$L_1(\Omega)$")
    makeErrPlot(t, nx, errS, "advection_errS", r"$L_1(P/\rho^\gamma)$")


def makeErrPlot(t, nx, err, name, label):
    print(name)
    print(err)
    T = np.unique(t)
    NX = np.unique(nx)
    figname = name + "_nx.png"
    if len(NX) > 1:
        fig, ax = plt.subplots(1, 1)
        for tt in T:
            ind = t == tt
            ax.plot(nx[ind], err[ind], marker='+', ms=10, mew=1, ls='')
        nn = np.logspace(math.log10(NX[1]), math.log10(NX[-1]),
                         num=10, base=10.0)
        ax.plot(nn, err.max() * np.power(nn/nn[0], -2), ls='--', lw=2,
                color='grey')
        ax.set_xlabel(r"$n_x$")
        ax.set_ylabel(label)
        ax.set_xscale('log')
        if (err > 0).any():
            ax.set_yscale('log')
            ylim = ax.get_ylim()
            if ylim[1] < 10*ylim[0]:
                ax.set_ylim(0.25*ylim[0], 4*ylim[1])
        print("Saving " + figname)
        fig.savefig(figname)
        plt.close(fig)

    figname = name + "_t.png"
    if len(T) > 1:
        fig, ax = plt.subplots(1, 1)
        for nn in NX:
            ind = nn == nx
            ax.plot(t[ind], err[ind])
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(label)
        ax.set_xscale('log')
        if (err > 0).any():
            ax.set_yscale('log')
            ylim = ax.get_ylim()
            if ylim[1] < 10*ylim[0]:
                ax.set_ylim(0.25*ylim[0], 4*ylim[1])
        print("Saving " + figname)
        fig.savefig(figname)
        plt.close(fig)


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need a checkpoint dude!")
        sys.exit()

    filenames = sys.argv[1:]
    analyze(filenames)
