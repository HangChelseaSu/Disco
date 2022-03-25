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

    if pars['Energy_Omega'] == 2:
        omk = np.sqrt(1/(r*r*r))
        en = 0.5*rho*(vr**2 + (r*(om-omk))**2 + vz**2) + P/(gam - 1)
    else:
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
    t, rs, zs, s, sg, sv, ss, sc, sd, dt, rjph, zkph = util.loadSource(filename)

    dz = zkph[-1] - zkph[0]

    rf = rf[0, :]
    f = f[0, :] / dz
    fv = fv[0, :] / dz
    s = s[0, :] / dz
    sg = sg[0, :] / dz
    sv = sv[0, :] / dz
    ss = ss[0, :] / dz
    sc = sc[0, :] / dz
    sd = sd[0, :] / dz

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

    return t, R, U, rf, f, fv, s, sg, sv, ss, sc, sd, V
    


if __name__ == "__main__":

    filenames = sys.argv[1:]

    n = len(filenames)

    t = []
    U = []
    f = []
    fv = []
    s = []
    sg = []
    sv = []
    ss = []
    sc = []
    sd = []

    for filename in filenames:
        ti, R, Ui, rf, fi, fvi, si, sgi, svi, ssi, sci, sdi, V\
                = analyze(Path(filename))
        t.append(ti)
        U.append(Ui)
        f.append(fi)
        fv.append(fvi)
        s.append(si)
        sg.append(sgi)
        sv.append(svi)
        ss.append(ssi)
        sc.append(sci)
        sd.append(sdi)

    t = np.array(t)
    U = np.array(U)
    f = np.array(f)
    fv = np.array(fv)
    s = np.array(s)
    sg = np.array(sg)
    sv = np.array(sv)
    ss = np.array(ss)
    sd = np.array(sd)
    sc = np.array(sc)

    f[0, :] = 0.0
    fv[0, :] = 0.0
    s[0, :] = 0.0
    sg[0, :] = 0.0
    sv[0, :] = 0.0
    ss[0, :] = 0.0
    sc[0, :] = 0.0
    sd[0, :] = 0.0

    dt = t[1:] - t[:-1]

    fig, ax = plt.subplots(4, 3, figsize=(12, 10))

    Nr = R.shape[0]

    j1 = Nr//4
    j2 = 3*Nr//4

    Ut = U[:, j1:j2, :].sum(axis=1)
    Ft = f[:, j1-1, :] - f[:, j2-1, :]
    FVt = fv[:, j1-1, :] - fv[:, j2-1, :]
    St = s[:, j1:j2, :].sum(axis=1)
    SGt = sg[:, j1:j2, :].sum(axis=1)
    SVt = sv[:, j1:j2, :].sum(axis=1)
    SSt = ss[:, j1:j2, :].sum(axis=1)
    SCt = sc[:, j1:j2, :].sum(axis=1)
    SDt = sd[:, j1:j2, :].sum(axis=1)

    Fdt = (Ft + FVt)[1:, :] * dt[:, None]
    Sdt = (St + SGt + SVt + SSt + SCt + SDt)[1:, :] * dt[:, None]
    dU = Ut[1:, :] - Ut[:-1, :]

    print(U.shape, f.shape, fv.shape)
    print(Ut.shape, Ft.shape, FVt.shape)
    print(Fdt.shape, dU.shape)

    ax[0, 0].plot(t[1:], dU[:, 0], '.', color='C0')
    ax[0, 0].plot(t[1:], Fdt[:, 0], 'x', color='C0')
    ax[0, 0].plot(t[1:], Sdt[:, 0], '+', color='C0')
    ax[0, 1].plot(t[1:], dU[:, 1], '.', color='C0')
    ax[0, 1].plot(t[1:], Fdt[:, 1], 'x', color='C0')
    ax[0, 1].plot(t[1:], Sdt[:, 1], '+', color='C0')
    # if U.shape[2] > 5:
    #     ax[0, 2].plot(t[1:], dU[:, 5], '.', color='C0')
    #     ax[0, 2].plot(t[1:], Fdt[:, 5], 'x', color='C0')
    ax[0, 2].plot(R, V)
    ax[0, 2].plot(R[1:-1], 2*np.pi*(0.5*(rf[1:]+rf[:-1]))*(rf[1:]-rf[:-1]))
    ax[1, 0].plot(t[1:], dU[:, 0]-Fdt[:, 0]-Sdt[:, 0], '.', color='C0')
    ax[1, 1].plot(t[1:], dU[:, 1]-Fdt[:, 1]-Sdt[:, 1], '.', color='C0')
    # if U.shape[2] > 5:
    #     ax[0, 2].plot(t[1:], dU[:, 5], '.', color='C0')
    #     ax[0, 2].plot(t[1:], Fdt[:, 5], 'x', color='C0')
    ax[1, 2].plot(R[1:-1], V[1:-1] - 2*np.pi*(0.5*(rf[1:]+rf[:-1]))*(rf[1:]-rf[:-1]))
    ax[2, 0].plot(t[1:], dU[:, 2], '.', color='C0')
    ax[2, 0].plot(t[1:], Fdt[:, 2], 'x', color='C0')
    ax[2, 0].plot(t[1:], Sdt[:, 2], '+', color='C0')
    ax[2, 1].plot(t[1:], dU[:, 3], '.', color='C0')
    ax[2, 1].plot(t[1:], Fdt[:, 3], 'x', color='C0')
    ax[2, 1].plot(t[1:], Sdt[:, 3], '+', color='C0')
    ax[2, 2].plot(t[1:], dU[:, 4], '.', color='C0')
    ax[2, 2].plot(t[1:], Fdt[:, 4], 'x', color='C0')
    ax[2, 2].plot(t[1:], Sdt[:, 4], '+', color='C0')
    ax[3, 0].plot(t[1:], dU[:, 2] - Fdt[:, 2] - Sdt[:, 2], '.', color='C0')
    ax[3, 1].plot(t[1:], dU[:, 3] - Fdt[:, 3] - Sdt[:, 3], '.', color='C0')
    ax[3, 2].plot(t[1:], dU[:, 4] - Fdt[:, 4] - Sdt[:, 4], '.', color='C0')

    figname = "flux_comp_dU.png"
    print("Saving", figname)
    fig.savefig(figname)
    plt.close(fig)



