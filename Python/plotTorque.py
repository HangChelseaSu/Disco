from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.geom as geom


def getCheckpointStuff(filename):

    pars = util.loadPars(filename)
    opts = util.loadOpts(filename)
    t, rf, zf, f, fv, dt, rjph, zkph = util.loadFluxR(filename)
    t, rs, zs, s, sg, sv, ss, sc, sd, dt, rjph, zkph = util.loadSource(filename)
    t, r, phi, z, prim, dat = util.loadCheckpoint(filename)
    _, _, _, diag, _, _ = util.loadDiagRZ(filename)
    
    Nr = len(rjph)-3
    dr = rjph[1:] - rjph[:-1]

    JdotG = sg[0, :Nr, 3].sum()
    MdotS = ss[0, :Nr, 0].sum()
    JdotS = ss[0, :Nr, 3].sum()

    dV_ann = 2*np.pi * (0.5*(rjph[1:]+rjph[:-1])) * dr

    dT_diag = diag[0, :, 4] * dV_ann

    JdotG_diag = dT_diag[:Nr].sum()


    Mdot_flux = f[0, Nr-1, 0]
    Jdot_flux = f[0, Nr-1, 3] + fv[0, Nr-1, 3]
    Rmax = rjph[Nr]

    rho = prim[:, 0]
    om = prim[:, 3]
    dV = geom.getDV(dat, opts, pars)

    M = (rho*dV)[r < Rmax].sum()
    J = (rho*r*r*om*dV)[r < Rmax].sum()

    name = (filename.stem).split("_")[-1]

    fig, ax = plt.subplots(2, 2)


    ax[0, 0].plot(rs[0, :-2], sg[0, :Nr, 3] / dr[:Nr], label='Grav Source')
    ax[0, 0].plot(rs[0, :-2], dT_diag[:Nr] / dr[:Nr], label='Diagnostics')
    ax[1, 0].plot(rs[0, :-2], (sg[0, :Nr, 3] - dT_diag[:Nr]) / dr[:Nr],
                label='Difference')

    ax[0, 1].axhline(JdotG, ls='--', alpha=0.5, color='C0', label='Grav Source')
    ax[0, 1].axhline(JdotG_diag, ls='--', alpha=0.5, color='C1', label='Diagnostics')
    ax[1, 1].axhline(JdotG-JdotG_diag, ls='--', alpha=0.5, color='C0', label='Difference')

    ax[0, 0].legend()
    ax[1, 0].legend()
    ax[0, 1].legend()
    ax[1, 1].legend()

    ax[0, 0].set(xlabel=r'$R$', ylabel=r'$dT_g/dR$')
    ax[1, 0].set(xlabel=r'$R$', ylabel=r'$\delta\ dT_g/dR$')
    ax[0, 1].set(xlabel=r'$R$', ylabel=r'$T_g$')
    ax[1, 1].set(xlabel=r'$R$', ylabel=r'$\delta T_g$')

    fig.tight_layout()
    figname = "torque_prof_{0:s}.png".format(name)
    print("Saving", figname)
    fig.savefig(figname)
    plt.close(fig)



    return t, MdotS, JdotG, JdotS, Mdot_flux, Jdot_flux, M, J, JdotG_diag
    


if __name__ == "__main__":

    reportname = Path(sys.argv[1])

    chkname = [Path(f) for f in sys.argv[2:]]

    t, Mg, Jg, dLg1, dLg2, dM1, dM2, dLa1, dLa2, dLs1, dLs2, dEg1, dEg2 = np.loadtxt(reportname,
            unpack=True, usecols=[0, 1, 4, 7, 8, 9, 10, 11, 12, 13, 14, 17, 18])

    tc = []
    dMgas = []
    dJggas = []
    dJagas = []
    dMf = []
    dJf = []
    Mc = []
    Jc = []
    dJggas_diag = []

    for name in chkname:
        ti, dMi, dJgi, dJai, dMfi, dJfi, Mi, Ji, dJgi_diag\
                = getCheckpointStuff(name)

        tc.append(ti)
        dMgas.append(dMi)
        dJggas.append(dJgi)
        dJagas.append(dJai)
        dMf.append(dMfi)
        dJf.append(dJfi)
        Mc.append(Mi)
        Jc.append(Ji)
        dJggas_diag.append(dJgi_diag)

    tc = np.array(tc)
    dMgas = np.array(dMgas)
    dJggas = np.array(dJggas)
    dJagas = np.array(dJagas)
    dJggas_diag = np.array(dJggas_diag)

    dMf = np.array(dMf)
    dJf = np.array(dJf)
    Mc = np.array(Mc)
    Jc = np.array(Jc)
    dMf[0] = 0.0
    dJf[0] = 0.0
    dMgas[0] = 0.0
    dJggas[0] = 0.0
    dJagas[0] = 0.0
    dJggas_diag[0] = 0.0

    tc1 = 0.5*(tc[1:] + tc[:-1]) / (2*np.pi)
    dtc = tc[1:] - tc[:-1]

    Mdot_gas = dMgas[1:]
    Jdotg_gas = dJggas[1:]
    Jdota_gas = dJagas[1:]
    Jdot_gas = Jdotg_gas + Jdota_gas
    Jdotg_gas_diag = dJggas_diag[1:]

    t1 = 0.5*(t[1:]+t[:-1]) / (2*np.pi)
    dt = t[1:] - t[:-1]

    Mdot1 = dM1[1:] / dt
    Mdot2 = dM2[1:] / dt

    Jdotg1 = (dLg1)[1:] / dt
    Jdotg2 = (dLg2)[1:] / dt
    Jdota1 = (dLa1 + dLs1)[1:] / dt
    Jdota2 = (dLa2 + dLs2)[1:] / dt
    Jdot1 = (dLg1 + dLa1 + dLs1)[1:] / dt
    Jdot2 = (dLg2 + dLa2 + dLs2)[1:] / dt

    Edotg1 = dEg1[1:] / dt
    Edotg2 = dEg2[1:] / dt

    fig, ax = plt.subplots(3, 3, figsize=(15, 10))
    print(dMf)
    print(dMgas)
    ax[0, 0].plot(t/(2*np.pi), Mg, label='report.dat')
    ax[0, 0].plot(tc[1:]/(2*np.pi),
                  Mg[0] + np.cumsum((-dMf[1:] + dMgas[1:])*dtc), '+',
                  label='FluxAvg + SourceAvg')
    ax[0, 0].plot(tc/(2*np.pi), Mc, '.', label='checkpoints')
    ax[1, 0].plot(t/(2*np.pi), Jg)
    ax[1, 0].plot(tc[1:]/(2*np.pi),
            Jg[0] + np.cumsum((-dJf + dJagas + dJggas)[1:] * dtc), '+')
    ax[1, 0].plot(tc/(2*np.pi), Jc, '.')


    ax[0, 1].plot(t/(2*np.pi), np.cumsum(dM1), label='report - P1')
    ax[0, 1].plot(t/(2*np.pi), np.cumsum(dM2), label='report - P2')
    ax[0, 1].plot(t/(2*np.pi), np.cumsum(dM1 + dM2), label='report - Tot')
    ax[0, 1].plot(tc[1:]/(2*np.pi), -np.cumsum(Mdot_gas * dtc), '.',
                  label='SourceAvg')

    ax[1, 1].plot(t/(2*np.pi), np.cumsum(dLg1))
    ax[1, 1].plot(t/(2*np.pi), np.cumsum(dLg2))
    ax[1, 1].plot(t/(2*np.pi), np.cumsum(dLg1 + dLg2))
    ax[1, 1].plot(tc[1:]/(2*np.pi), -np.cumsum(Jdotg_gas * dtc), '.')
    ax[1, 1].plot(tc[1:]/(2*np.pi), -np.cumsum(Jdotg_gas_diag * dtc), 'x')

    ax[2, 1].plot(t/(2*np.pi), np.cumsum(dLa1 + dLs1))
    ax[2, 1].plot(t/(2*np.pi), np.cumsum(dLa2 + dLs2))
    ax[2, 1].plot(t/(2*np.pi), np.cumsum(dLa1 + dLa2 + dLs1 + dLs2))
    ax[2, 1].plot(tc[1:]/(2*np.pi), -np.cumsum(Jdota_gas * dtc), '.')

    ax[0, 2].plot(t1, Mdot1)
    ax[0, 2].plot(t1, Mdot2)
    ax[0, 2].plot(t1, Mdot1 + Mdot2)
    ax[0, 2].plot(tc1, -Mdot_gas, '.')
    
    ax[1, 2].plot(t1, Jdotg1)
    ax[1, 2].plot(t1, Jdotg2)
    ax[1, 2].plot(t1, Jdotg1 + Jdotg2)
    ax[1, 2].plot(tc1, -Jdotg_gas, '.')
    ax[1, 2].plot(tc1, -Jdotg_gas_diag, 'x')
    
    ax[2, 2].plot(t1, Jdota1)
    ax[2, 2].plot(t1, Jdota2)
    ax[2, 2].plot(t1, Jdota1 + Jdota2)
    ax[2, 2].plot(tc1, -Jdota_gas, '.')

    ax[0, 0].legend()
    ax[0, 1].legend()


    ax[0, 0].set(xlabel=r'$t$ (orbits)', ylabel=r'$M_{\rm gas}$')
    ax[1, 0].set(xlabel=r'$t$ (orbits)', ylabel=r'$J_{\rm gas}$')
    ax[0, 1].set(xlabel=r'$t$ (orbits)', ylabel=r'$\delta M_{\rm bin}$')
    ax[1, 1].set(xlabel=r'$t$ (orbits)', ylabel=r'$\delta J_{\rm bin, grav}$')
    ax[2, 1].set(xlabel=r'$t$ (orbits)', ylabel=r'$\delta J_{\rm bin, acc}$')
    ax[0, 2].set(xlabel=r'$t$ (orbits)', ylabel=r'$\dot{M}$')
    ax[1, 2].set(xlabel=r'$t$ (orbits)', ylabel=r'$\dot{J}_{\rm grav}$')
    ax[2, 2].set(xlabel=r'$t$ (orbits)', ylabel=r'$\dot{J}_{\rm acc}$')

    fig.tight_layout()

    fig.savefig("sum.pdf")

    #plt.show()
