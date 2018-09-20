import sys
import math
import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import discoGeom as dg

def loadPars(filename):

    f = h5.File(filename, "r")

    pars = {}
    for key in f['Pars']:
        pars[key] = f['Pars'][key][0]

    f.close()

    return pars

def loadOpts(filename):

    f = h5.File(filename, "r")

    opts = {}
    for key in f['Opts']:
        opts[key] = f['Opts'][key][0]

    f.close()

    return opts

def loadCheckpointPrims(filename):
   
    f = h5.File(filename, "r")
    NUM_C = f['Opts']['NUM_C'][0]
    NUM_N = f['Opts']['NUM_N'][0]
    NUM_Q = NUM_N + NUM_C
    prim = f['Data']['Cells'][:,:NUM_Q][...]
    f.close()
    return prim

def loadCheckpoint(filename):

    f = h5.File(filename, "r")

    NUM_C = f['Opts']['NUM_C'][0]
    NUM_N = f['Opts']['NUM_N'][0]
    NUM_Q = NUM_N + NUM_C

    piph = f['Data']['Cells'][:,-1][...]
    prim = f['Data']['Cells'][:,:NUM_Q][...]
    Phi = f['Data']['Cells'][:,NUM_Q:-1][...]
    index = f['Grid']['Index'][...]
    idPhi0 = f['Grid']['Id_phi0'][...]
    nphi = f['Grid']['Np'][...]
    index = f['Grid']['Index'][...]
    t = f['Grid']['T'][0]
    riph = f['Grid']['r_jph'][...]
    ziph = f['Grid']['z_kph'][...]
    planetDat = f['Data']['Planets'][...]
    pmax = f['Pars']['Phi_Max']

    primPhi0 = np.zeros((index.shape[0], index.shape[1], prim.shape[1]))
    for k in range(index.shape[0]):
        for j in range(index.shape[1]):
            primPhi0[k,j,:] = prim[idPhi0[k,j],:]

    dat = (riph, ziph, primPhi0, piph, planetDat, Phi, idPhi0, nphi, index)
    opts = loadOpts(filename)
    pars = loadPars(filename)

    r, phi, z = dg.getCellCentroids(dat, opts, pars)

    return t, r, phi, z, prim, dat

def loadDiagRZ(filename):

    f = h5.File(filename, "r")

    t = f['Grid']['T'][0]
    rjph = f['Grid']['r_jph'][...]
    zkph = f['Grid']['z_kph'][...]
    diag = f['Data']['Diagnostics'][...]

    f.close()

    Nr = rjph.shape[0]-1
    Nz = zkph.shape[0]-1
    Nq = diag.shape[-1]
    #diag = np.resize(diag, (Nz,Nr,Nq))

    R = 0.5*(rjph[1:]+rjph[:-1])
    Z = 0.5*(zkph[1:]+zkph[:-1])

    r = np.empty((Nz,Nr))
    z = np.empty((Nz,Nr))

    r[:,:] = R[None,:]
    z[:,:] = Z[:,None]

    return t, r, z, diag, rjph, zkph

def plotAx(ax, x, y, xscale, yscale, xlabel, ylabel, *args, **kwargs):
    ax.plot(x, y, *args, **kwargs)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    if (y>0).any():
        ax.set_yscale(yscale)
    else:
        ax.set_yscale("linear")

def getVarNames(filename):
    f = h5.File(filename, "r")
    hydro = str(f['Opts']['HYDRO'][0])
    num_c = int(f['Opts']['NUM_C'][0])
    num_n = int(f['Opts']['NUM_N'][0])
    num_q = num_c + num_n

    names = None
    texnames = None

    if hydro == 'euler' and num_c <= 5:
        names = ['rho', 'P', 'vr', 'om', 'vz'][:num_c]
        texnames = [r'$\rho$', r'$P$', r'$v_r$', r'$\Omega$', r'$v_z$'][:num_c]
        for q in range(num_n):
            names.append('q{0:d}'.format(q))
            texnames.append(r'$q_{0:d}$'.format(q))
    
    elif hydro == 'euler_cart' and num_c <= 5:
        names = ['rho', 'P', 'vx', 'vy', 'vz'][:num_c]
        texnames = [r'$\rho$', r'$P$', r'$v_x$', r'$v_y$', r'$v_z$'][:num_c]
        for q in range(num_n):
            names.append('q{0:d}'.format(q))
            texnames.append(r'$q_{0:d}$'.format(q))
    
    elif hydro == 'euler_sph' and num_c <= 5:
        names = ['rho', 'P', 'vr', 'om', 'vt'][:num_c]
        texnames = [r'$\rho$', r'$P$', r'$v_r$', r'$\Omega$', r'$v^\theta$'][:num_c]
        for q in range(num_n):
            names.append('q{0:d}'.format(q))
            texnames.append(r'$q_{0:d}$'.format(q))
    
    elif hydro == 'mhd' and num_c <= 8:
        names = ['rho', 'P', 'vr', 'om', 'vz', 'Br', 'Bp', 'Bz'][:num_c]
        texnames = [r'$\rho$', r'$P$', r'$v_r$', r'$\Omega$', r'$v_z$', 
                        r'$B_r$', r'$B_\phi$', r'$B_z$'][:num_c]
        for q in range(num_n):
            names.append('q{0:d}'.format(q))
            texnames.append(r'$q_{0:d}$'.format(q))

    elif hydro == 'mhd_cart' and num_c <= 8:
        names = ['rho', 'P', 'vx', 'vy', 'vz', 'Bx', 'By', 'Bz'][:num_c]
        texnames = [r'$\rho$', r'$P$', r'$v_x$', r'$v_y$', r'$v_z$', 
                        r'$B_x$', r'$B_y$', r'$B_z$'][:num_c]
        for q in range(num_n):
            names.append('q{0:d}'.format(q))
            texnames.append(r'$q_{0:d}$'.format(q))

    elif hydro == 'greuler' and num_c <= 5:
        names = ['rho', 'P', 'u_r', 'u_p', 'u_z'][:num_c]
        texnames = [r'$\rho$', r'$P$', r'$u_r$', r'$u_\phi$', r'$u_z$'][:num_c]
        for q in range(num_n):
            names.append('q{0:d}'.format(q))
            texnames.append(r'$q_{0:d}$'.format(q))

    elif (hydro == 'grmhd' or hydro == 'grmhd2') and num_c <= 8:
        names = ['rho', 'P', 'u_r', 'u_p', 'u_z', 'Br', 'Bp', 'Bz'][:num_c]
        texnames = [r'$\rho$', r'$P$', r'$u_r$', r'$u_\phi$', r'$u_z$', 
                        r'$B^r$', r'$B^\phi$', r'$B^z$'][:num_c]
        for q in range(num_n):
            names.append('q{0:d}'.format(q))
            texnames.append(r'$q_{0:d}$'.format(q))

    if names == None or texnames == None:
        names = ["p{0:d}".format(q) for q in range(num_c)]
        texnames = [r"$p_{0:d}$".format(q) for q in range(num_c)]
        for q in range(num_n):
            names.append('q{0:d}'.format(q))
            texnames.append(r'$q_{0:d}$'.format(q))

    return names, texnames, num_c, num_n

