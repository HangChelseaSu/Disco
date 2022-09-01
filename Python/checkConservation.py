from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.geom as geom

class Cons:

    def __init__(self, N=10):
        self.size = N
        self.N = 0
        self.t = np.zeros(N)
        self.M = np.zeros(N)
        self.MX = np.zeros((N, 3))
        self.P = np.zeros((N, 3))
        self.L = np.zeros(N)
        self.S = np.zeros(N)
        self.K = np.zeros(N)
        self.Eint = np.zeros(N)
        self.U = np.zeros(N)

    def push(self, t, M, MX, P, L, S, K, Eint, U):
        if self.N == self.size:
            self.size *= 2
            self.t.resize(self.size)
            self.M.resize(self.size)
            tempMX = self.MX.copy()
            tempP = self.P.copy()
            self.MX = np.zeros((self.size, 3))
            self.P = np.zeros((self.size, 3))
            self.MX[:self.N, :] = tempMX[:, :]
            self.P[:self.N, :] = tempP[:, :]
            self.L.resize(self.size)
            self.S.resize(self.size)
            self.K.resize(self.size)
            self.Eint.resize(self.size)
            self.U.resize(self.size)

        self.t[self.N] = t
        self.M[self.N] = M
        self.MX[self.N, :] = MX
        self.P[self.N, :] = P
        self.L[self.N] = L
        self.S[self.N] = S
        self.K[self.N] = K
        self.Eint[self.N] = Eint
        self.U[self.N] = U
        self.N += 1

    def finalize(self):
        self.size = self.N
        self.t.resize(self.N)
        self.M.resize(self.N)
        MXtemp = self.MX.copy()
        Ptemp = self.P.copy()
        self.MX = np.empty((self.N, 3))
        self.P = np.empty((self.N, 3))
        self.MX[:, :] = MXtemp[:self.N, :]
        self.P[:, :] = Ptemp[:self.N, :]
        self.L.resize(self.N)
        self.S.resize(self.N)
        self.K.resize(self.N)
        self.Eint.resize(self.N)
        self.U.resize(self.N)

    def integrate(self):
        self.M = self.M.cumsum()
        self.MX = self.MX.cumsum(axis=0)
        self.P = self.P.cumsum(axis=0)
        self.L = self.L.cumsum()
        self.S = self.S.cumsum()
        self.K = self.K.cumsum()
        self.U = self.U.cumsum()
        self.Eint = self.Eint.cumsum()

    def __add__(self, other):
        tot = Cons(self.N)
        tot.N = self.N
        tot.t = self.t
        tot.M = self.M + other.M
        tot.MX = self.MX + other.MX
        tot.P = self.P + other.P
        tot.L = self.L + other.L
        tot.S = self.S + other.S
        tot.K = self.K + other.K
        tot.Eint = self.Eint + other.Eint
        tot.U = self.U + other.U

        return tot

    @property
    def J(self):
        return self.L + self.S

    @property
    def E(self):
        return self.K + self.Eint + self.U

    @property
    def X(self):
        return self.MX / self.M

    @property
    def x(self):
        return self.MX[:, 0] / self.M

    @property
    def y(self):
        return self.MX[:, 1] / self.M

    @property
    def z(self):
        return self.MX[:, 2] / self.M

    @property
    def Px(self):
        return self.P[:, 0]

    @property
    def Py(self):
        return self.P[:, 1]

    @property
    def Pz(self):
        return self.P[:, 2]


def PhiGrav(r, phi, z, planet):
    
    rp = planet[3]
    phip = planet[4]
    m = planet[0]
    eps = float(planet[5])
    plType = int(planet[6])

    dr = np.sqrt(rp**2 + r**2 - 2*rp*r*np.cos(phip-phi) + z**2)
    dr = np.atleast_1d(dr)

    if plType == 3:
        
        eps *= 2.8
        u = dr / eps
        val = np.empty(dr.shape)

        inn = u < 0.5
        mid = (u >= 0.5) & (u < 1)
        out = u >= 1

        if inn.any():
            val[inn] = (16/3*u*u - 48/5*u**4 + 32/5*u**5 - 14/5)[inn]
        if mid.any():
            val[mid] = (1/(15*u) + 32/3*u**2 - 16*u**3 + 48/5*u**4
                        - 32/15*u**5 - 3.2)[mid]
        if out.any():
            val[out] = (-1/u)[out]

        Phi = m * val/eps
    else:
        Phi = -m / np.sqrt(dr*dr + eps*eps)

    return Phi


def calcGasCons(r, phi, z, prim, dat, pars, opts):
                                                
    Rmin = pars['R_Min']
    Rmax = pars['R_Max']
    Zmin = pars['Z_Min']
    Zmax = pars['Z_Max']
    gam = pars['Adiabatic_Index']

    inn = (r > Rmin) & (r < Rmax) & (z > Zmin) & (z < Zmax)

    dV = geom.getDV(dat, opts, pars)
    dV[~inn] = 0.0

    rho = prim[:, 0]
    P = prim[:, 1]
    vr = prim[:, 2]
    om = prim[:, 3]
    vz = prim[:, 4]

    cos = np.cos(phi)
    sin = np.sin(phi)

    vx = cos * vr - sin * r*om
    vy = sin * vr + cos * r*om

    x = r * cos
    y = r * sin

    v2 = vr*vr + r*r*om*om + vz*vz

    rhoe = P / (gam - 1)

    planetDat = dat[4]
    Npl = planetDat.shape[0]

    Phi = np.zeros(r.shape)
    for i in range(Npl):
        Phi += PhiGrav(r, phi, z, planetDat[i])

    M = geom.integrate(rho, dat, opts, pars, dV)
    Mx = geom.integrate(rho*x, dat, opts, pars, dV)
    My = geom.integrate(rho*y, dat, opts, pars, dV)
    Mz = geom.integrate(rho*z, dat, opts, pars, dV)
    PX = geom.integrate(rho*vx, dat, opts, pars, dV)
    PY = geom.integrate(rho*vy, dat, opts, pars, dV)
    PZ = geom.integrate(rho*vz, dat, opts, pars, dV)

    MX = np.array([Mx, My, Mz])
    P = np.array([PX, PY, PZ])

    L = geom.integrate(rho*r*r*om, dat, opts, pars, dV)
    K = geom.integrate(0.5*rho*v2, dat, opts, pars, dV)
    Eint = geom.integrate(rhoe, dat, opts, pars, dV)
    Ug = geom.integrate(rho*Phi, dat, opts, pars, dV)

    return M, MX, P, L, K, Eint, Ug

def calcGasConsFlux(r, phi, z, prim, rf, zf, zkph, fr, fvr, dt, dat,
                    pars, opts):

    M = 0.0
    MX = np.zeros(3)
    P = np.zeros(3)
    L = 0.0
    K = 0.0
    Eint = 0.0
    Ug = 0.0

    nq = fr.shape[-1]

    dz = zkph[-1] - zkph[0]
    fr = fr[0, :, :] / dz
    fvr = fvr[0, :, :] / dz

    if dt > 0.0:
        f = np.zeros(nq)
        if pars['NoBC_Rmin'] == 0:
            f -= fr[1, :] + fvr[1, :]
        if pars['NoBC_Rmax'] == 0:
            f += fr[-2, :] + fvr[-2, :]
        fdt = f * dt
        F_M = fdt[0]
        F_J = fdt[3]
        F_E = fdt[1]

        M -= F_M
        L -= F_J
        K -= 0.5*F_E
        Eint -= 0.5*F_E

    return M, MX, P, L, K, Eint, Ug

def calcGasConsSource(r, phi, z, prim, rs, zs, zkph, sH, sV, sC, sD, dt,
        dat, pars, opts):

    M = 0.0
    MX = np.zeros(3)
    P = np.zeros(3)
    L = 0.0
    K = 0.0
    Eint = 0.0
    Ug = 0.0

    Nr = sH.shape[1]
    nq = sH.shape[2]

    if dt > 0.0:
        dz = zkph[-1] - zkph[0]
        s = (sH + sV + sC + sD)[0, :, :]

        j1 = 2 if pars['NoBC_Rmin'] == 0 else 0
        j2 = Nr-2 if pars['NoBC_Rmax'] == 0 else Nr

        Sdt = s[j1:j2, :].sum(axis=0) * dt/dz
        
        S_M = Sdt[0]
        S_J = Sdt[3]
        S_E = Sdt[1]

        M += S_M
        L += S_J
        K += 0.5*S_E
        Eint -= 0.5*S_E

    return M, MX, P, L, K, Eint, Ug

def calcGasConsSourceGrav(r, phi, z, prim, rs, zs, zkph, sG, dt,
        dat, pars, opts):

    M = 0.0
    MX = np.zeros(3)
    P = np.zeros(3)
    L = 0.0
    K = 0.0
    Eint = 0.0
    Ug = 0.0

    Nr = sG.shape[1]

    if dt > 0.0:
        dz = zkph[-1] - zkph[0]

        j1 = 2 if pars['NoBC_Rmin'] == 0 else 0
        j2 = Nr-2 if pars['NoBC_Rmax'] == 0 else Nr

        Sdt = sG[0, j1:j2, :].sum(axis=0) * dt/dz
        
        S_M = Sdt[0]  # Should be zero!
        S_J = Sdt[3]
        S_E = Sdt[1]

        M += S_M
        L += S_J
        Ug += S_E

    return M, MX, P, L, K, Eint, Ug

def calcGasConsSourceSink(r, phi, z, prim, rs, zs, zkph, sS, dt,
        dat, pars, opts):

    M = 0.0
    MX = np.zeros(3)
    P = np.zeros(3)
    L = 0.0
    K = 0.0
    Eint = 0.0
    Ug = 0.0

    Nr = sS.shape[1]

    if dt > 0.0:
        dz = zkph[-1] - zkph[0]

        j1 = 2 if pars['NoBC_Rmin'] == 0 else 0
        j2 = Nr-2 if pars['NoBC_Rmax'] == 0 else Nr

        Sdt = sS[0, j1:j2, :].sum(axis=0) * dt/dz
        
        S_M = Sdt[0]
        S_J = Sdt[3]
        S_Pz = Sdt[4]
        S_E = Sdt[1]

        M += S_M
        L += S_J
        P[2] += S_Pz
        K += 0.5*S_E
        Eint -= 0.5*S_E

    return M, MX, P, L, K, Eint, Ug

def calcPlCons(dat, pars, opts):

    planetDat = dat[4]
    Npl = planetDat.shape[0]

    Mtot = 0.0
    MXtot = np.zeros(3)
    Ptot = np.zeros(3)
    Ltot = 0.0
    Stot = 0.0
    Ktot = 0.0
    Einttot = 0.0
    Utot = 0.0

    for i in range(Npl):
        pl = planetDat[i]
        M, vr, om, r, phi, eps, plType = pl[:7]
        kin = pl[7:]
        Sz = kin[7]
        Eint = kin[8]

        cos = np.cos(phi)
        sin = np.sin(phi)

        x = np.array([r*cos, r*sin, 0.0])
        v = np.array([vr*cos - r*om*sin, vr*sin + r*om*cos, 0.0])

        Mtot += M
        MXtot += M*x
        Ptot += M*v
        Ltot += M*r*r*om
        Ktot += 0.5*M*(vr*vr + r*r*om*om)
        Stot += Sz
        Einttot += Eint
        
        for j in range(Npl):
            if j == i:
                continue
            Utot += 0.5*M*PhiGrav(r, phi, 0, planetDat[j])

    return Mtot, MXtot, Ptot, Ltot, Stot, Ktot, Einttot, Utot


def calcCons(filename):

    pars = util.loadPars(filename)
    opts = util.loadOpts(filename)
    t, r, phi, z, prim, dat = util.loadCheckpoint(filename)
    _, rf, zf, fr, fvr, dt, _, zkph = util.loadFluxR(filename)
    _, rs, zs, sH, sG, sV, sS, sC, sD, _, _, _ = util.loadSource(filename)

    Mg, Xg, Pg, Lg, Kg, Eintg, Ug = calcGasCons(r, phi, z, prim,
                                                dat, pars, opts)
    Mf, Xf, Pf, Lf, Kf, Eintf, Uf = calcGasConsFlux(r, phi, z, prim,
                                                rf, zf, zkph, fr, fvr, dt,
                                                dat, pars, opts)
    Ms, Xs, Ps, Ls, Ks, Eints, Us = calcGasConsSource(r, phi, z, prim,
                                                rs, zs, zkph,
                                                sH, sV, sC, sD, dt,
                                                dat, pars, opts)
    Msg, Xsg, Psg, Lsg, Ksg, Eintsg, Usg = calcGasConsSourceGrav(r, phi, z,
                                                prim, rs, zs, zkph, sG, dt,
                                                dat, pars, opts)
    Mss, Xss, Pss, Lss, Kss, Eintss, Uss = calcGasConsSourceSink(r, phi, z,
                                                prim, rs, zs, zkph, sS, dt,
                                                dat, pars, opts)
    Mp, Xp, Pp, Lp, Sp, Kp, Eintp, Up = calcPlCons(dat, pars, opts)

    consDatGas = (Mg, Xg, Pg, Lg, 0, Kg, Eintg, Ug)
    consDatFlux = (Mf, Xf, Pf, Lf, 0, Kf, Eintf, Uf)
    consDatPl = (Mp, Xp, Pp, Lp, Sp, Kp, Eintp, Up)
    consDatSrc = (Ms, Xs, Ps, Ls, 0, Ks, Eints, Us)
    consDatSrcGrav = (Msg, Xsg, Psg, Lsg, 0, Ksg, Eintsg, Usg)
    consDatSrcSink = (Mss, Xss, Pss, Lss, 0, Kss, Eintss, Uss)


    return t, consDatGas, consDatFlux, consDatPl, consDatSrc, consDatSrcGrav,\
            consDatSrcSink



def runCheckpoints(checkpointFilenames, report):

    N = len(checkpointFilenames)

    t = np.empty(N)

    consGas = Cons()
    consFlux = Cons()
    consPl = Cons()
    consSrc = Cons()
    consSrcGrv = Cons()
    consSrcSnk = Cons()

    for i in range(N):
        print("Loading", checkpointFilenames[i])
        ti, consDatGas, consDatFlux, consDatPl, consDatSrc, consDatSrcGrav,\
                consDatSrcSink = calcCons(checkpointFilenames[i])
        consGas.push(ti, *consDatGas) 
        consFlux.push(ti, *consDatFlux) 
        consPl.push(ti, *consDatPl)
        consSrc.push(ti, *consDatSrc)
        consSrcGrv.push(ti, *consDatSrcGrav)
        consSrcSnk.push(ti, *consDatSrcSink)
        t[i] = ti
        print(consFlux.M[consFlux.N-1])

    consGas.finalize()
    consPl.finalize()
    consFlux.finalize()
    consSrc.finalize()
    consSrcGrv.finalize()
    consSrcSnk.finalize()

    consFlux.integrate()
    consSrc.integrate()
    consSrcGrv.integrate()
    consSrcSnk.integrate()

    #print(consFlux.M)
    #print(consFlux.J)

    print("Plotting")

    consTot = consGas + consPl

    fig, ax = plt.subplots(6, 3, figsize=(12, 12))

    ax[0, 0].plot(t, consTot.M, label='Total')
    ax[0, 0].plot(t, consGas.M, label='Gas')
    ax[0, 0].plot(t, consPl.M, label='Planet')
    ax[0, 0].plot(t, consFlux.M, label='Flux')
    ax[0, 0].plot(t, consSrc.M, label='Src')
    ax[0, 0].plot(t, consSrcGrv.M, label='Src-Grav')
    ax[0, 0].plot(t, consSrcSnk.M, label='Src-Snk')
    ax[1, 0].plot(t, consTot.M - consTot.M[0], label='Total')
    ax[1, 0].plot(t, consFlux.M, label='Flux')

    ax[0, 1].plot(t, consTot.J)
    ax[0, 1].plot(t, consGas.J)
    ax[0, 1].plot(t, consPl.J)
    ax[0, 1].plot(t, consFlux.J)
    ax[1, 1].plot(t, consTot.J - consTot.J[0])
    ax[1, 1].plot(t, consFlux.J)

    ax[0, 2].plot(t, consTot.E)
    ax[0, 2].plot(t, consGas.K + consGas.Eint)
    ax[0, 2].plot(t, consPl.E)
    ax[0, 2].plot(t, consFlux.E)
    ax[0, 2].plot(t, consGas.U)
    ax[1, 2].plot(t, consTot.E - consTot.E[0])
    ax[1, 2].plot(t, consFlux.E)

    ax[2, 0].plot(t, consTot.x)
    ax[2, 0].plot(t, consGas.x)
    ax[2, 0].plot(t, consPl.x)
    ax[3, 0].plot(t, consTot.x - consTot.x[0])

    ax[2, 1].plot(t, consTot.y)
    ax[2, 1].plot(t, consGas.y)
    ax[2, 1].plot(t, consPl.y)
    ax[3, 1].plot(t, consTot.y - consTot.y[0])

    ax[2, 2].plot(t, consTot.z)
    ax[2, 2].plot(t, consGas.z)
    ax[2, 2].plot(t, consPl.z)
    ax[3, 2].plot(t, consTot.z - consTot.z[0])

    ax[4, 0].plot(t, consTot.Px)
    ax[4, 0].plot(t, consGas.Px)
    ax[4, 0].plot(t, consPl.Px)
    ax[5, 0].plot(t, consTot.Px - consTot.Px[0])

    ax[4, 1].plot(t, consTot.Py)
    ax[4, 1].plot(t, consGas.Py)
    ax[4, 1].plot(t, consPl.Py)
    ax[5, 1].plot(t, consTot.Py - consTot.Py[0])

    ax[4, 2].plot(t, consTot.Pz)
    ax[4, 2].plot(t, consGas.Pz)
    ax[4, 2].plot(t, consPl.Pz)
    ax[5, 2].plot(t, consTot.Pz - consTot.Pz[0])

    ax[0, 0].set(xlabel='t', ylabel='Mass')
    ax[1, 0].set(xlabel='t', ylabel='Mass Change')
    
    ax[0, 1].set(xlabel='t', ylabel='Angular Momentum')
    ax[1, 1].set(xlabel='t', ylabel='Angular Momentum Change')
    
    ax[0, 2].set(xlabel='t', ylabel='Energy')
    ax[1, 2].set(xlabel='t', ylabel='Energy Change')
    
    ax[2, 0].set(xlabel='t', ylabel=r'$X_{COM}$')
    ax[3, 0].set(xlabel='t', ylabel=r'$X_{COM}$ Change')
    
    ax[2, 1].set(xlabel='t', ylabel=r'$Y_{COM}$')
    ax[3, 1].set(xlabel='t', ylabel=r'$Y_{COM}$ Change')
    
    ax[2, 2].set(xlabel='t', ylabel=r'$Z_{COM}$')
    ax[3, 2].set(xlabel='t', ylabel=r'$Z_{COM}$ Change')
    
    ax[4, 0].set(xlabel='t', ylabel=r'$P_X$')
    ax[5, 0].set(xlabel='t', ylabel=r'$P_X$ Change')
    
    ax[4, 1].set(xlabel='t', ylabel=r'$P_Y$')
    ax[5, 1].set(xlabel='t', ylabel=r'$P_Y$ Change')
    
    ax[4, 2].set(xlabel='t', ylabel=r'$P_Z$')
    ax[5, 2].set(xlabel='t', ylabel=r'$P_Z$ Change')

    ax[0, 0].legend()
    ax[1, 0].legend()

    fig.tight_layout()

    figname = "conservation_check_total.png"
    print("Saving", figname)
    fig.savefig(figname)
    plt.close(fig)


    Kpl = ((report.KIN_PR**2 + (report.KIN_LL/report.KIN_R)**2)
           / (2*report.KIN_M)).sum(axis=0)

    dt_rep = np.diff(report.t)
    dt = np.diff(t)
    tc = 0.5*(t[1:] + t[:-1])
    tc_rep = 0.5*(report.t[1:] + report.t[:-1])

    dKpl = np.diff(Kpl, prepend=0.0)
    dKpl_dt = dKpl[1:] / dt_rep

    dKpl_aux = report.GRV_K.sum(axis=0) + report.SNK_K.sum(axis=0)
    dKpl_dt_aux = dKpl_aux[1:] / dt_rep

    pl_PX = np.cos(report.KIN_PHI) * report.KIN_PR\
            - np.sin(report.KIN_PHI) * report.KIN_LL / report.KIN_R
    pl_PY = np.sin(report.KIN_PHI) * report.KIN_PR\
            + np.cos(report.KIN_PHI) * report.KIN_LL / report.KIN_R
    pl_VX = pl_PX / report.KIN_M
    pl_VY = pl_PY / report.KIN_M
    pl_VZ = report.KIN_PZ / report.KIN_M

    pl_VX[1:] = 0.5*(pl_VX[1:] + pl_VX[:-1])
    pl_VY[1:] = 0.5*(pl_VY[1:] + pl_VY[:-1])
    pl_VZ[1:] = 0.5*(pl_VZ[1:] + pl_VZ[:-1])
    
    dKpl_aux2 = (pl_VX*(report.GRV_PX+report.SNK_PX)
                 + pl_VY*(report.GRV_PY+report.SNK_PY)
                 + pl_VZ*(report.GRV_PZ+report.SNK_PZ)).sum(axis=0)
    dKpl_dt_aux2 = dKpl_aux2[1:] / dt_rep

    dKpl_chk = np.diff(consPl.K, prepend=0.0)
    dKpl_dt_chk = dKpl_chk[1:] / dt

    Uf = report.Uf.sum(axis=0)
    dUf = np.diff(Uf, prepend=0.0)
    dUf_dt = dUf[1:] / dt_rep

    dUf_chk = np.diff(consGas.U, prepend=0.0)
    dUf_dt_chk = dUf_chk[1:] / dt

    dUf_aux = report.GRV_U.sum(axis=0) + report.SNK_U.sum(axis=0)
    dUf_dt_aux = dUf_aux[1:] / dt_rep


    fig, ax = plt.subplots(3, 5, figsize=(15, 8))

    ax[0, 0].plot(report.t, report.cons[0], '-')
    ax[0, 0].plot(t, consGas.M, '.')
    ax[0, 1].plot(report.t, report.cons[3], '-')
    ax[0, 1].plot(t, consGas.J, '.')
    ax[0, 2].plot(report.t, report.cons[1], '-')
    ax[0, 2].plot(t, consGas.K+consGas.Eint, '.')
    
    ax[0, 3].plot(report.t, Kpl, '-', label='report')
    ax[0, 3].plot(t, consPl.K, '.', label='checkpoints')
    ax[0, 4].plot(report.t, Uf, '-', label='report')
    ax[0, 4].plot(t, consGas.U, '.', label='checkpoints')

    ax[1, 3].plot(tc_rep, dKpl_dt, '-')
    ax[1, 3].plot(tc, dKpl_dt_chk, '.')
    ax[1, 3].plot(tc_rep, dKpl_dt_aux, '-')
    ax[1, 3].plot(tc_rep, dKpl_dt_aux2, '-')

    ax[1, 4].plot(tc_rep, dUf_dt, '-')
    ax[1, 4].plot(tc, dUf_dt_chk, '.')
    ax[1, 4].plot(tc_rep, dUf_dt_aux, '-')

    ax[0, 0].set(xlabel=r'$t$', ylabel=r'$M_f$')
    ax[0, 1].set(xlabel=r'$t$', ylabel=r'$J_f$')
    ax[0, 2].set(xlabel=r'$t$', ylabel=r'$E_f$')
    ax[0, 3].set(xlabel=r'$t$', ylabel=r'$K_p$')
    ax[0, 4].set(xlabel=r'$t$', ylabel=r'$U_g$')

    fig.tight_layout()

    figname = "planet_evol_check.png"
    print("Saving", figname)
    fig.savefig(figname)
    plt.close(fig)



def runReport(report):

    fig, ax = plt.subplots(3, 1)

    colors = ['C{0:d}'.format(i) for i in range(10)]

    j_grv = report.GRV_JZ.cumsum(axis=1)
    l_grv = report.GRV_LZ.cumsum(axis=1)
    l_kin = report.KIN_LL - report.KIN_LL[:, 0:1]

    for i in range(report.Npl):
        ax[0].plot(report.t, j_grv[i], color=colors[i], ls='--', lw=1)
        ax[0].plot(report.t, l_grv[i], color=colors[i], ls=':', lw=2)
        ax[0].plot(report.t, l_kin[i], color=colors[i], ls='-', lw=1)

    ax[0].plot(report.t, j_grv.sum(axis=0), color='k', ls='--', lw=1)
    ax[0].plot(report.t, l_grv.sum(axis=0), color='k', ls=':', lw=2)
    ax[0].plot(report.t, l_kin.sum(axis=0), color='k', ls='-', lw=1)

    for i in range(report.Npl):
        ax[1].plot(report.t, j_grv[i]-l_grv[i], color=colors[i], ls='-', lw=1)

    ax[1].plot(report.t, (j_grv-l_grv).sum(axis=0), color='k', ls='-', lw=1)
    
    ax[2].plot(report.t, (j_grv-l_kin).sum(axis=0), color='k', ls='--', lw=1)



    figname = "conservation_report_check.png"
    print("Saving", figname)
    fig.savefig(figname)
    plt.close(fig)


def reintegrate_timeseries(t1, df1, t2):

    N = len(t2)
    df2 = np.zeros(N)
    for j in range(1, N):
        df2[j] = integrate_timeseries(t1, df1, t2[j-1], t2[j])

    return df2


def integrate_timeseries(t, df, t1, t2):

    if t1 < t[0]:
        t1 = t[0]

    if t2 > t[-1]:
        t2 = t[-1]

    N = len(t)
    dt = t[1:] - t[:-1]

    j1 = np.searchsorted(t, t1)
    if j1 < 1:
        j1 = 1

    j2 = np.searchsorted(t, t2)
    if j2 > N-1:
        j2 = N-1

    if j1 > j2:
        return 0.0
    
    if j1 == j2:
        return df[j1] * (t2-t1) / dt[j1-1]

    part1 = df[j1] * (t[j1] - t1) / dt[j1-1]
    part2 = df[j2] * (t2 - t[j2-1]) / dt[j2-1]

    total = part1 + part2
    if j2 > j1 + 1:
        total += df[j1+1:j2].sum()

    return total


def makePlotStack(ax, name, t, F, DF_cmps=None, DF_plts=None):

    dt = t[1:] - t[:-1]
    tc = 0.5*(t[1:] + t[:-1])

    f, label = F

    fdot = np.diff(f) / dt

    ax[0].plot(t, f, label=label, color='k', lw=2)
    ax[1].plot(tc, fdot, color='k', lw=2)

    if DF_cmps is not None:
        for df_cmp, label_cmp in DF_cmps:
            f_cmp = np.concatenate(([f[0]], f[0] + np.cumsum(df_cmp)))
            fdot_cmp = df_cmp / dt

            ax[0].plot(t, f_cmp, label=label_cmp)
            ax[1].plot(tc, fdot_cmp)
            ax[2].plot(t, f_cmp - f)
            ax[3].plot(tc, fdot_cmp - fdot)
    
    if DF_plts is not None:
        for df_plt, label_plt in DF_plts:
            f_plt = np.concatenate(([f[0]], f[0] + np.cumsum(df_plt)))
            fdot_plt = df_plt / dt

            ax[0].plot(t, f_plt, label=label_plt)
            ax[1].plot(tc, fdot_plt)

    ax[0].set(xlabel=r'$t$', ylabel=r"${0:s}$".format(name))
    ax[1].set(xlabel=r'$t$', ylabel=r"$\dot{{{0:s}}}$".format(name))
    ax[2].set(xlabel=r'$t$', ylabel=r"$\delta {0:s}$".format(name))
    ax[3].set(xlabel=r'$t$', ylabel=r"$\delta \dot{{{0:s}}}$".format(name))

    ax[0].legend(fontsize=8)


def makePlanetFigures(report):

    t = report.t

    ax_w = 3
    ax_h = 8/3

    for p in range(report.Npl):

        fig, ax = plt.subplots(4, 10, figsize=(10*ax_w, 4*ax_h))
        figname = "conservation_check_planet_{0:d}.png".format(p)
        print("Plotting planet", p)

        M = report.KIN_M[p]
        dM_snk = report.SNK_M[p][1:]

        print(t.shape)
        print(report.KIN_M.shape)
        print(M.shape)
        print(report.SNK_M.shape)
        print(dM_snk.shape)

        r = report.KIN_R[p]
        phi = report.KIN_PHI[p]
        cos = np.cos(phi)
        sin = np.sin(phi)

        Pr = report.KIN_PR[p]
        L = report.KIN_LL[p]
        Pz = report.KIN_PZ[p]
        S = report.KIN_SZ[p]

        J = L + S

        Px = cos * Pr - sin * L/r
        Py = sin * Pr + cos * L/r

        dPx_grv = report.GRV_PX[p][1:]
        dPy_grv = report.GRV_PY[p][1:]
        dPz_grv = report.GRV_PZ[p][1:]
        dPx_snk = report.SNK_PX[p][1:]
        dPy_snk = report.SNK_PY[p][1:]
        dPz_snk = report.SNK_PZ[p][1:]
        dPx_ext = report.EXT_PX[p][1:]
        dPy_ext = report.EXT_PY[p][1:]
        dPz_ext = report.EXT_PZ[p][1:]
        
        dPx_aux = dPx_grv + dPx_snk + dPx_ext
        dPy_aux = dPy_grv + dPy_snk + dPy_ext
        dPz_aux = dPz_grv + dPz_snk + dPz_ext

        dJ_grv = report.GRV_JZ[p][1:]
        dJ_snk = report.SNK_JZ[p][1:]
        dS_snk = report.SNK_SZ[p][1:]
        dL_snk = report.SNK_LZ[p][1:]
        dJ_ext = report.EXT_JZ[p][1:]
        dJ_snk2 = dL_snk + dS_snk
        dL_snk2 = dJ_snk - dS_snk
        dS_snk2 = dJ_snk - dL_snk

        dJ_aux = dJ_grv + dJ_snk + dJ_ext
        dJ_aux2 = dJ_grv + dJ_snk2 + dJ_ext
        dL_aux = dJ_grv + dL_snk + dJ_ext
        dL_aux2 = dJ_grv + dL_snk2 + dJ_ext
        dS_aux = dS_snk
        dS_aux2 = dS_snk2

        K = (Pr**2 + (L/r)**2 + Pz**2) / (2*M)
        Q = report.KIN_EINT[p]
        Uf = report.Uf[p]

        dK_grv = report.GRV_K[p][1:]
        dK_snk = report.SNK_K[p][1:]
        dK_ext = report.EXT_K[p][1:]
        dK_aux = dK_grv + dK_snk + dK_ext

        dQ_snk = report.SNK_EINT[p][1:]

        dUf_grv = report.GRV_U[p][1:]
        dUf_snk = report.SNK_UGAS[p][1:]
        
        dUf_aux = dUf_grv + dUf_snk

        makePlotStack(ax[:, 0], "M", t, (M, "Actual"),
                        [(dM_snk, "Tracked")])
        makePlotStack(ax[:, 1], "P_x", t, (Px, "Actual"),
                        [(dPx_aux, "Tracked")],
                        [(dPx_grv, "grv"), (dPx_snk, "s"), (dPx_ext, "ext")])
        makePlotStack(ax[:, 2], "P_y", t, (Py, "Actual"),
                        [(dPy_aux, "Tracked")],
                        [(dPy_grv, "grv"), (dPy_snk, "snk"), (dPy_ext, "ext")])
        makePlotStack(ax[:, 3], "P_z", t, (Pz, "Actual"),
                        [(dPz_aux, "Tracked")],
                        [(dPz_grv, "grv"), (dPz_snk, "snk"), (dPz_ext, "ext")])
        makePlotStack(ax[:, 4], "J", t, (J, "Actual"),
                        [(dJ_aux, "Tracked"), (dJ_aux2, "Tracked L+S")],
                        [(dJ_grv, "grv"), (dJ_ext, "ext"),
                            (dL_snk, "snk-L"), (dS_snk, "snk-S")])
        makePlotStack(ax[:, 5], "L", t, (L, "Actual"),
                        [(dL_aux, "Tracked"), (dL_aux2, "Tracked J-S")])
        makePlotStack(ax[:, 6], "S", t, (S, "Actual"),
                        [(dS_aux, "Tracked"), (dS_aux2, "Tracked J-L")])
        makePlotStack(ax[:, 7], "K", t, (K, "Actual"),
                        [(dK_aux, "Tracked")],
                        [(dK_grv, "grv"), (dK_snk, "snk"), (dK_ext, "ext")])
        makePlotStack(ax[:, 8], "U_f", t, (Uf, "Actual"),
                        [(dUf_aux, "Tracked")],
                        [(dUf_grv, "grv"), (dUf_snk, "snk")])
        makePlotStack(ax[:, 9], "Q", t, (Q, "Actual"),
                        [(dQ_snk, "Tracked")])


        fig.tight_layout()
        print("Saving", figname)
        fig.savefig(figname, dpi=100)
        plt.close(fig)




if __name__ == "__main__":

    reportFilename = Path(sys.argv[1])
    checkpointFilenames = [Path(x) for x in sys.argv[2:]]
    
    report = util.DiscoReportTestLive(reportFilename)

    runCheckpoints(checkpointFilenames, report)

    runReport(report)

    makePlanetFigures(report)
