import sys
import math
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
import discopy.util as util
import discopy.geom as geom
import discopy.plot as plot


def analyzeSingle(filename):

    opts = util.loadOpts(filename)
    pars = util.loadPars(filename)
    print("Loading " + filename)
    t, x1, x2, x3, prim, dat = util.loadCheckpoint(filename)

    rjph = dat[0]
    piph = dat[3]

    nx = pars['Num_R']
    gam = pars['Adiabatic_Index']

    rho = prim[:, 0]
    P = prim[:, 1]
    v1 = prim[:, 2]
    v2 = prim[:, 3]
    v3 = prim[:, 4]
    q = prim[:, 5]

    x, y, z = geom.getXYZ(x1, x2, x3, opts, pars)
    vx, vy, vz = geom.getVXYZ(x1, x2, x3, v1, v2, v3, opts)

    X = np.linspace(-rjph[-1], rjph[-1], 500)
    Y = np.linspace(-rjph[-1], rjph[-1], 495)

    XX, YY = np.meshgrid(X, Y)

    vxXY = np.empty(XX.shape)
    vyXY = np.empty(XX.shape)

    print("Interpolating") 
    vxXY.flat[:] = scipy.interpolate.griddata((x, y), vx, (XX.flat, YY.flat))
    vyXY.flat[:] = scipy.interpolate.griddata((x, y), vy, (XX.flat, YY.flat))

    dy = 0.5
    rout = 0.5*(rjph[-3]+rjph[-2])
    yyy = np.linspace(-dy, dy, 10)
    xxx = -0.5 * np.ones(yyy.shape)
    start_points = np.array([xxx, yyy]).T

    print("Plotting") 
    fig, ax = plt.subplots(1, 1)
    plot.plotZSlice(fig, ax, rjph, piph, x1, rho, x3.mean(), r"$\rho$",
                    pars, opts)
    ax.streamplot(X, Y, vxXY, vyXY, start_points=start_points, density=10,
                  integration_direction='both', minlength=1.0)
    ax.set_aspect('equal')
    print("Saving") 
    fig.savefig('stream.png')
    plt.close(fig)

    print("Calculating div")

    divV = geom.calculateDivV(x1, x2, x3, v1, v2, v3, dat, opts, pars)

    interior = (x1 > pars['R_Min']) & (x1 < pars['R_Max'])
    minDiv = divV[interior].min()
    maxDiv = divV[interior].max()

    print("Plotting") 
    fig, ax = plt.subplots(1, 1)
    plot.plotZSlice(fig, ax, rjph, piph, x1, divV, x3.mean(),
                    r"$\nabla\cdot v$", pars, opts, vmin=minDiv, vmax=maxDiv)
    ax.streamplot(X, Y, vxXY, vyXY, start_points=start_points, density=10,
                  integration_direction='both', minlength=1.0)
    ax.set_aspect('equal')
    print("Saving") 
    fig.savefig('stream_divV.pdf')
    plt.close(fig)

    return t, nx

def analyze(filenames):

    N = len(filenames)

    t = np.empty(N)
    nx = np.empty(N)

    for i, f in enumerate(filenames):
        dat = analyzeSingle(f)
        t[i] = dat[0]
        nx[i] = dat[1]


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need a checkpoint dude!")
        sys.exit()

    filenames = sys.argv[1:]
    analyze(filenames)
