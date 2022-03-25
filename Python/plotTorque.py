from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
import discopy.util as util


def getCheckpointStuff(filename):

    t, rs, zs, s, sg, sv, ss, sc, sd, dt, rjph, zkph = util.loadSource(filename)

    JdotG = sg[0, :, 3].sum()
    MdotS = ss[0, :, 0].sum()
    JdotS = ss[0, :, 3].sum()

    return t, MdotS, JdotG, JdotS
    


if __name__ == "__main__":

    reportname = Path(sys.argv[1])

    chkname = [Path(f) for f in sys.argv[2:]]

    t, Mg, Jg, dM1, dM2, dLg1, dLg2, dLa1, dLa2 = np.loadtxt(reportname,
            unpack=True, usecols=[0, 1, 4, 9, 10, 7, 8, 11, 12])

    tc = []
    dMgas = []
    dJggas = []
    dJagas = []

    for name in chkname:
        ti, dMi, dJgi, dJai = getCheckpointStuff(name)

        tc.append(ti)
        dMgas.append(dMi)
        dJggas.append(dJgi)
        dJagas.append(dJai)

    tc = np.array(tc)
    dMgas = np.array(dMgas)
    dJggas = np.array(dJggas)
    dJagas = np.array(dJagas)

    tc1 = 0.5*(tc[1:] + tc[:-1]) / (2*np.pi)
    dtc = tc[1:] - tc[:-1]

    Mdot_gas = dMgas[1:] / dtc
    Jdotg_gas = dJggas[1:] / dtc
    Jdota_gas = dJagas[1:] / dtc
    Jdot_gas = Jdotg_gas + Jdota_gas

    print(Mdot_gas)
    print(Jdot_gas)

    t1 = t[1:] / (2*np.pi)

    Mdot1 = dM1[1:] / (t[1:] - t[:-1])
    Mdot2 = dM2[1:] / (t[1:] - t[:-1])

    Jdot1 = (dLg1)[1:] / (t[1:] - t[:-1])
    Jdot2 = (dLg2)[1:] / (t[1:] - t[:-1])
    Jdot1 = (dLg1 + dLa1)[1:] / (t[1:] - t[:-1])
    Jdot2 = (dLg2 + dLa2)[1:] / (t[1:] - t[:-1])

    fig, ax = plt.subplots(2, 1, figsize=(12, 9))

    ax[0].plot(t1, Mdot1)
    ax[0].plot(t1, Mdot2)
    ax[0].plot(t1, Mdot1 + Mdot2)
    ax[0].plot(tc1, -Mdot_gas)
    ax[1].plot(t1, Jdot1)
    ax[1].plot(t1, Jdot2)
    ax[1].plot(t1, Jdot1 + Jdot2)
    ax[1].plot(tc1, -Jdot_gas)

    ax[0].set(xlabel=r'$t$ (orbits)', ylabel=r'$\dot{M}$')
    ax[1].set(xlabel=r'$t$ (orbits)', ylabel=r'$\dot{J}$')

    fig.savefig("sum.png")
