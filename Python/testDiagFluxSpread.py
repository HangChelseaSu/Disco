from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.plot as plot
import discopy.geom as geom

def prim2cons(prim, r, phi, z, pars):

    cons = np.empty_like(prim)

    rho = prim[:, 0]
    P = prim[:, 1]
    vr = prim[:, 2]
    om = prim[:, 3]
    vz = prim[:, 4]

    gam = pars['Adiabatic_Index']

    en = 0.5*rho*(vr**2 + (r*om)**2 + vz**2) + P/(gam - 1)

    cons[:, 0] = rho
    cons[:, 1] = en
    cons[:, 2] = rho*vr
    cons[:, 3] = rho*r*r*om
    cons[:, 4] = rho*vz

    nq = prim.shape[1]

    for q in range(5, nq):
        cons[:, q] = rho * prim[:, q]

    return cons


def analyze(filename):

    pars = util.loadPars(filename)
    opts = util.loadOpts(filename)
    t, r, phi, z, prim, dat = util.loadCheckpoint(filename)

    t, rf, zf, f, fv, dt, rjph, zkph = util.loadFluxR(filename)

    dz = zkph[-1] - zkph[0]

    rf = rf[0, :]
    f = f[0, :] / dz
    fv = fv[0, :] / dz

    rho = prim[:, 0]
    P = prim[:, 1]
    vr = prim[:, 2]
    om = prim[:, 3]

    nu = pars['Viscosity']

    f_rho = 2*np.pi*r*rho*vr
    f_Sr = 2*np.pi*r*(rho*vr*vr + P)
    f_Sp = 2*np.pi*r*(rho*vr*r*r*om)
    fv_Sp = 2*np.pi*r*(1.5*nu * rho *r*om)

    fig, ax = plt.subplots(2, 2, figsize=(8, 6))
    ax[0, 0].plot(rf[1:-1], f[1:-1, 0])
    ax[0, 0].plot(rf[1:-1], fv[1:-1, 0])
    ax[0, 0].plot(r, f_rho, '.', color='C0')
    ax[0, 1].plot(rf[1:-1], f[1:-1, 1])
    ax[0, 1].plot(rf[1:-1], fv[1:-1, 1])
    ax[1, 0].plot(rf[1:-1], f[1:-1, 2])
    ax[1, 0].plot(rf[1:-1], fv[1:-1, 2])
    ax[1, 0].plot(r, f_Sr, '.', color='C0')
    ax[1, 1].plot(rf[1:-1], f[1:-1, 3])
    ax[1, 1].plot(rf[1:-1], fv[1:-1, 3])
    ax[1, 1].plot(r, f_Sp, '.', color='C0')
    ax[1, 1].plot(r, fv_Sp, '.', color='C1')

    ax[0, 0].set(xscale='log')
    ax[0, 1].set(xscale='log')
    ax[1, 0].set(xscale='log')
    ax[1, 1].set(xscale='log')

    name = (filename.stem).split('_')[-1]

    figname = "flux_comp_{0:s}.png".format(name)
    print("Saving", figname)
    fig.savefig(figname)
    plt.close(fig)

    cons = prim2cons(prim, r, phi, z, pars)
    dV = geom.getDV(dat, opts, pars)

    Nr = len(rjph)-1
    Nq = prim.shape[1]
    U = np.empty((Nr, Nq))
    R = np.empty(Nr)
    V = np.empty(Nr)

    for j in range(Nr):
        idx = (r>rjph[j]) & (r<rjph[j+1])
        R[j] = r[idx].mean()
        U[j, :] = (cons[idx, :]*dV[idx, None]).sum(axis=0)
        V[j] = dV[idx].sum()

    return t, R, U, rf, f, fv, V
    


if __name__ == "__main__":

    filenames = sys.argv[1:]

    n = len(filenames)

    t = []
    U = []
    f = []
    fv = []

    for filename in filenames:
        ti, R, Ui, rf, fi, fvi, V = analyze(Path(filename))
        t.append(ti)
        U.append(Ui)
        f.append(fi)
        fv.append(fvi)

    t = np.array(t)
    U = np.array(U)
    f = np.array(f)
    fv = np.array(fv)
    f[0, :] = 0.0
    fv[0, :] = 0.0

    dt = t[1:] - t[:-1]

    fig, ax = plt.subplots(4, 3, figsize=(12, 10))

    Nr = R.shape[0]

    j1 = Nr//4
    j2 = 3*Nr//4

    Ut = U[:, j1:j2, :].sum(axis=1)
    Ft = f[:, j1-1, :] - f[:, j2-1, :]
    FVt = fv[:, j1-1, :] - fv[:, j2-1, :]

    Fdt = (Ft + FVt)[1:, :] * dt[:, None]
    dU = Ut[1:, :] - Ut[:-1, :]

    print(U.shape, f.shape, fv.shape)
    print(Ut.shape, Ft.shape, FVt.shape)
    print(Fdt.shape, dU.shape)

    ax[0, 0].plot(t[1:], dU[:, 0], '.', color='C0')
    ax[0, 0].plot(t[1:], Fdt[:, 0], 'x', color='C0')
    ax[0, 1].plot(t[1:], dU[:, 1], '.', color='C0')
    ax[0, 1].plot(t[1:], Fdt[:, 1], 'x', color='C0')
    # if U.shape[2] > 5:
    #     ax[0, 2].plot(t[1:], dU[:, 5], '.', color='C0')
    #     ax[0, 2].plot(t[1:], Fdt[:, 5], 'x', color='C0')
    ax[0, 2].plot(R, V)
    ax[0, 2].plot(R[1:-1], 2*np.pi*(0.5*(rf[1:]+rf[:-1]))*(rf[1:]-rf[:-1]))
    ax[1, 0].plot(t[1:], dU[:, 0]-Fdt[:, 0], '+', color='C0')
    ax[1, 1].plot(t[1:], dU[:, 1]-Fdt[:, 1], '+', color='C0')
    # if U.shape[2] > 5:
    #     ax[0, 2].plot(t[1:], dU[:, 5], '.', color='C0')
    #     ax[0, 2].plot(t[1:], Fdt[:, 5], 'x', color='C0')
    ax[1, 2].plot(R[1:-1], V[1:-1] - 2*np.pi*R[1:-1]*(rf[1:]-rf[:-1]))
    ax[2, 0].plot(t[1:], dU[:, 2], '.', color='C0')
    ax[2, 0].plot(t[1:], Fdt[:, 2], 'x', color='C0')
    ax[2, 1].plot(t[1:], dU[:, 3], '.', color='C0')
    ax[2, 1].plot(t[1:], Fdt[:, 3], 'x', color='C0')
    ax[2, 2].plot(t[1:], dU[:, 4], '.', color='C0')
    ax[2, 2].plot(t[1:], Fdt[:, 4], 'x', color='C0')
    ax[3, 0].plot(t[1:], dU[:, 2] - Fdt[:, 2], '+', color='C0')
    ax[3, 1].plot(t[1:], dU[:, 3] - Fdt[:, 3], '+', color='C0')
    ax[3, 2].plot(t[1:], dU[:, 4] - Fdt[:, 4], '+', color='C0')

    figname = "flux_comp_dU.png"
    print("Saving", figname)
    fig.savefig(figname)
    plt.close(fig)



