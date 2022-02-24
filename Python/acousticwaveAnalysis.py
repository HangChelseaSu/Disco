import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.geom as geom


def rhoProf(x, a, x0, L, rho0):

    rho = np.empty(x.shape)
    rho[:] = rho0
    inn = np.fabs(x-x0) < L
    xin = x[inn]
    f = (xin-x0)*(xin-x0)/(L*L) - 1.0
    rho[inn] *= 1.0 + a*f*f*f*f

    return rho


def calcAcousticWave(t, x, pars, tol=1.0e-6):
    rho0 = 1.0
    cs0 = 1.0
    gam = pars['Adiabatic_Index']
    mach = pars['Init_Par1']
    a = pars['Init_Par2']
    x0 = pars['Init_Par3']
    L = pars['Init_Par4']

    xR = x.max()
    xL = x.min()

    v0 = mach * cs0

    rhomax = rho0 * (1+a)
    csmax = cs0*np.power(rhomax/rho0, 0.5*(gam-1))
    vmax = v0 + 2*(csmax-cs0)/(gam-1)

    xa = x - (1.1*(vmax+csmax))*t
    xb = x - (0.9*(v0+cs0))*t
    err = np.empty(x.shape)
    err[:] = np.inf

    while np.fabs(err).min() > tol:
        xc = 0.5*(xa+xb)
        rho = rhoProf(xc, a, x0, L, rho0)
        cs = cs0*np.power(rho/rho0, 0.5*(gam-1))
        v = v0 + 2*(cs-cs0)/(gam-1)
        xf = xc + (v+cs)*t
        err = (xf-x) / (xR-xL)
        over = err > 0
        under = err < 0
        xb[over] = xc[over]
        xa[under] = xc[under]

    xc = 0.5*(xa+xb)
    rho = rhoProf(xc, a, x0, L, rho0)
    cs = cs0*np.power(rho/rho0, 0.5*(gam-1))
    v = v0 + 2*(cs-cs0)/(gam-1)

    return rho, cs, v


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
    cs = np.sqrt(gam*P/rho)

    x, y, z = geom.getXYZ(x1, x2, x3, opts, pars)
    vx, vy, vz = geom.getVXYZ(x1, x2, x3, v1, v2, v3, opts)

    gam = pars['Adiabatic_Index']
    phi0 = pars['Init_Par5'] * np.pi
    cos_theta0 = pars['Init_Par6']
    sin_theta0 = math.sqrt(1.0 - cos_theta0*cos_theta0)
    kx = math.cos(phi0) * sin_theta0
    ky = math.sin(phi0) * sin_theta0
    kz = cos_theta0

    X = kx*x + ky*y + kz*z
    V = kx*vx + ky*vy + kz*vz
    Jm = V + 2*cs/(gam-1)
    Jp = V - 2*cs/(gam-1)

    rhoS, csS, VS = calcAcousticWave(t, X, pars, 1.0e-13)
    PS = rhoS*csS*csS/gam
    eSS = PS * np.power(rhoS, -gam)
    JmS = VS + 2*csS/(gam-1)
    JpS = VS - 2*csS/(gam-1)

    dV = geom.getDV(dat, opts, pars)

    maskDistance = 0.2

    dVmask = dV.copy()
    # dVmask[(x1 < pars['R_Min']+maskDistance)] = 0.0
    dVmask[(x1 > pars['R_Max']-maskDistance)] = 0.0

    errRho = geom.integrate(np.fabs(rho-rhoS), dat, opts, pars, dVmask)
    errP = geom.integrate(np.fabs(P-PS), dat, opts, pars, dVmask)
    errVx = geom.integrate(np.fabs(vx-kx*VS), dat, opts, pars, dVmask)
    errVy = geom.integrate(np.fabs(vy-ky*VS), dat, opts, pars, dVmask)
    errVz = geom.integrate(np.fabs(vz-kz*VS), dat, opts, pars, dVmask)
    errS = geom.integrate(np.fabs(eS-eSS), dat, opts, pars, dVmask)
    errJp = geom.integrate(np.fabs(Jp-JpS), dat, opts, pars, dVmask)
    errJm = geom.integrate(np.fabs(Jm-JmS), dat, opts, pars, dVmask)

    return t, nx, errRho, errP, errVx, errVy, errVz, errS, errJp, errJm


def analyze(filenames):

    N = len(filenames)

    t = np.empty(N)
    nx = np.empty(N)
    errRho = np.empty(N)
    errVx = np.empty(N)
    errVy = np.empty(N)
    errVz = np.empty(N)
    errP = np.empty(N)
    errS = np.empty(N)
    errJp = np.empty(N)
    errJm = np.empty(N)

    for i, f in enumerate(filenames):
        dat = analyzeSingle(f)
        t[i] = dat[0]
        nx[i] = dat[1]
        errRho[i] = dat[2]
        errP[i] = dat[3]
        errVx[i] = dat[4]
        errVy[i] = dat[5]
        errVz[i] = dat[6]
        errS[i] = dat[7]
        errJp[i] = dat[8]
        errJm[i] = dat[9]

    makeErrPlot(t, nx, errRho, "acousticwave_errRho", r"$L_1(\rho)$")
    makeErrPlot(t, nx, errP, "acousticwave_errP", r"$L_1(P)$")
    makeErrPlot(t, nx, errVx, "acousticwave_errVx", r"$L_1(v_x)$")
    makeErrPlot(t, nx, errVy, "acousticwave_errVy", r"$L_1(v_y)$")
    makeErrPlot(t, nx, errVz, "acousticwave_errVz", r"$L_1(v_z)$")
    makeErrPlot(t, nx, errS, "acousticwave_errS", r"$L_1(P/\rho^\gamma)$")
    makeErrPlot(t, nx, errJp, "acousticwave_errJp", r"$L_1(J_{+})$")
    makeErrPlot(t, nx, errJm, "acousticwave_errJm", r"$L_1(J_{-})$")


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
        ax.plot(nn, err.max() * np.power(nn/nn[0], -1), ls=':', lw=2,
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
