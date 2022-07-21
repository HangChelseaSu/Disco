import h5py as h5
import numpy as np
from . import geom as dg


def loadGitVersion(filename):

    f = h5.File(filename, "r")

    gv = f['GIT_VERSION'][0]

    f.close()

    return gv


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
        try:
            opts[key] = int(f['Opts'][key][0])
        except ValueError:
            opts[key] = str(f['Opts'][key][0], encoding='utf-8')

    if 'GEOMETRY' not in opts:
        opts['GEOMETRY'] = 'cylindrical'

    f.close()

    return opts


def loadCheckpointPrims(filename):

    f = h5.File(filename, "r")
    NUM_C = f['Opts']['NUM_C'][0]
    NUM_N = f['Opts']['NUM_N'][0]
    NUM_Q = NUM_N + NUM_C
    prim = f['Data']['Cells'][:, :NUM_Q][...]
    f.close()
    return prim


def loadCheckpoint(filename):

    f = h5.File(filename, "r")

    NUM_C = f['Opts']['NUM_C'][0]
    NUM_N = f['Opts']['NUM_N'][0]
    NUM_Q = NUM_N + NUM_C

    piph = f['Data']['Cells'][:, -1][...]
    prim = f['Data']['Cells'][:, :NUM_Q][...]
    Phi = f['Data']['Cells'][:, NUM_Q:-1][...]
    index = f['Grid']['Index'][...]
    idPhi0 = f['Grid']['Id_phi0'][...]
    nphi = f['Grid']['Np'][...]
    index = f['Grid']['Index'][...]
    t = f['Grid']['T'][0]
    riph = f['Grid']['r_jph'][...]
    ziph = f['Grid']['z_kph'][...]
    planetDat = f['Data']['Planets'][...]

    primPhi0 = np.zeros((index.shape[0], index.shape[1], prim.shape[1]))
    for k in range(index.shape[0]):
        for j in range(index.shape[1]):
            primPhi0[k, j, :] = prim[idPhi0[k, j], :]

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
    # diag = np.resize(diag, (Nz,Nr,Nq))

    R = 0.5*(rjph[1:]+rjph[:-1])
    Z = 0.5*(zkph[1:]+zkph[:-1])

    r = np.empty((Nz, Nr))
    z = np.empty((Nz, Nr))

    r[:, :] = R[None, :]
    z[:, :] = Z[:, None]

    return t, r, z, diag, rjph, zkph


def loadFluxR(filename):

    f = h5.File(filename, "r")

    t = f['Grid']['T'][0]
    rjph = f['Grid']['r_jph'][...]
    zkph = f['Grid']['z_kph'][...]
    fluxHydro = f['Data']['FluxHydroAvgR'][...]
    fluxVisc = f['Data']['FluxViscAvgR'][...]
    dt = f['Data']['Diagnostics_DT'][0]

    f.close()

    Nr = rjph.shape[0]-1
    Nz = zkph.shape[0]-1

    if Nr == 1:
        print("No radial faces for fluxes")
        return None

    opts = loadOpts(filename)

    Z = dg.getCentroid(zkph[:-1], zkph[1:], 2, opts)

    r = np.empty((Nz, Nr-1))
    z = np.empty((Nz, Nr-1))

    r[:, :] = rjph[None, 1:-1]
    z[:, :] = Z[:, None]

    return t, r, z, fluxHydro, fluxVisc, dt, rjph, zkph


def loadFluxZ(filename):

    f = h5.File(filename, "r")

    t = f['Grid']['T'][0]
    rjph = f['Grid']['r_jph'][...]
    zkph = f['Grid']['z_kph'][...]
    fluxHydro = f['Data']['FluxHydroAvgZ'][...]
    fluxVisc = f['Data']['FluxViscAvgZ'][...]

    f.close()

    Nr = rjph.shape[0]-1
    Nz = zkph.shape[0]-1

    if Nz == 1:
        print("No z-faces for fluxes")
        return None

    opts = loadOpts(filename)
    dt = 1.0

    R = dg.getCentroid(rjph[:-1], rjph[1:], 1, opts)

    r = np.empty((Nz-1, Nr))
    z = np.empty((Nz-1, Nr))

    r[:, :] = R[None, :]
    z[:, :] = zkph[1:-1, None]

    return t, r, z, fluxHydro, fluxVisc, dt, rjph, zkph


def loadSource(filename):

    f = h5.File(filename, "r")

    t = f['Grid']['T'][0]
    rjph = f['Grid']['r_jph'][...]
    zkph = f['Grid']['z_kph'][...]
    srcHydro = f['Data']['SourceHydroAvg'][...]
    srcGrav = f['Data']['SourceGravAvg'][...]
    srcVisc = f['Data']['SourceViscAvg'][...]
    srcSink = f['Data']['SourceSinkAvg'][...]
    srcCool = f['Data']['SourceCoolAvg'][...]
    srcDamp = f['Data']['SourceDampAvg'][...]

    f.close()

    Nr = rjph.shape[0]-1
    Nz = zkph.shape[0]-1

    opts = loadOpts(filename)
    dt = 1.0

    R = dg.getCentroid(rjph[:-1], rjph[1:], 1, opts)
    Z = dg.getCentroid(zkph[:-1], zkph[1:], 2, opts)

    r = np.empty((Nz, Nr))
    z = np.empty((Nz, Nr))

    r[:, :] = R[None, :]
    z[:, :] = Z[:, None]

    return t, r, z, srcHydro, srcGrav, srcVisc, srcSink, srcCool, srcDamp,\
            dt, rjph, zkph


def plotAx(ax, x, y, xscale, yscale, xlabel, ylabel, *args, **kwargs):
    ax.plot(x, y, *args, **kwargs)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    if (y > 0).any():
        ax.set_yscale(yscale)
    else:
        ax.set_yscale("linear")


def getVarNames(filename):
    f = h5.File(filename, "r")
    hydro = str(f['Opts']['HYDRO'][0], encoding='utf-8')
    num_c = int(f['Opts']['NUM_C'][0])
    num_n = int(f['Opts']['NUM_N'][0])
    f.close()

    names = None
    texnames = None

    if hydro == 'euler' and num_c <= 5:
        names = ['rho', 'P', 'vr', 'om', 'vz'][:num_c]
        texnames = [r'$\rho$', r'$P$', r'$v_r$', r'$\Omega$', r'$v_z$'
                    ][:num_c]
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
        texnames = [r'$\rho$', r'$P$', r'$v_r$', r'$\Omega$', r'$v^\theta$'
                    ][:num_c]
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

    elif hydro == 'mhd_sph' and num_c <= 8:
        names = ['rho', 'P', 'vr', 'om', 'vt', 'Br', 'Bp', 'Bt'][:num_c]
        texnames = [r'$\rho$', r'$P$', r'$v^r$', r'$\Omega$', r'$v^\theta$',
                    r'$B^r$', r'$B^\phi$', r'$B^\theta$'][:num_c]
        for q in range(num_n):
            names.append('q{0:d}'.format(q))
            texnames.append(r'$q_{0:d}$'.format(q))

    elif hydro == 'greuler' and num_c <= 5:
        names = ['rho', 'P', 'ur', 'up', 'uz'][:num_c]
        texnames = [r'$\rho$', r'$P$', r'$u_r$', r'$u_\phi$', r'$u_z$'][:num_c]
        for q in range(num_n):
            names.append('q{0:d}'.format(q))
            texnames.append(r'$q_{0:d}$'.format(q))

    elif (hydro == 'grmhd' or hydro == 'grmhd2') and num_c <= 8:
        names = ['rho', 'P', 'ur', 'up', 'uz', 'Br', 'Bp', 'Bz'][:num_c]
        texnames = [r'$\rho$', r'$P$', r'$u_r$', r'$u_\phi$', r'$u_z$',
                    r'$B^r$', r'$B^\phi$', r'$B^z$'][:num_c]
        for q in range(num_n):
            names.append('q{0:d}'.format(q))
            texnames.append(r'$q_{0:d}$'.format(q))

    if names is None or texnames is None:
        names = ["p{0:d}".format(q) for q in range(num_c)]
        texnames = [r"$p_{0:d}$".format(q) for q in range(num_c)]
        for q in range(num_n):
            names.append('q{0:d}'.format(q))
            texnames.append(r'$q_{0:d}$'.format(q))

    return names, texnames, num_c, num_n


class DiscoReport:

    def __init__(self, filename=None):

        if filename is None:
            return

        hdr = self._readHeader(filename)

        self.NUM_C = hdr['NUM_C']
        self.NUM_N = hdr['NUM_N']
        self.NUM_Q = self.NUM_C + self.NUM_N
        self.Npl = hdr['Npl']
        self.NUM_PL_KIN = hdr['NUM_PL_KIN']
        self.NUM_PL_AUX = hdr['NUM_PL_AUX']
        self.N_shared = hdr['N_shared_reports']
        self.N_dist_aux = hdr['N_distributed_aux_reports']
        self.N_dist_int = hdr['N_distributed_integral_reports']

        self.t = np.loadtxt(filename, unpack=True, usecols=0, comments='#')
        self.N = self.t.shape[0]

        self.cons = np.loadtxt(filename, comments='#',
                               usecols=range(1, self.NUM_Q+1)).T

        print(self.N_shared, self.N_dist_aux, self.N_dist_int)

        a_shared = 1+self.NUM_Q
        a_dist_aux = a_shared + self.N_shared
        a_dist_int = a_dist_aux + self.N_dist_aux

        print(a_shared, a_dist_aux, a_dist_int, a_dist_int+self.N_dist_int)

        self.shared = None
        if self.N_shared > 0:
            self.shared = np.loadtxt(filename, comments='#',
                    usecols=range(a_shared, a_shared+self.N_shared)).T
        
        self.dist_aux = None
        if self.N_dist_aux > 0:
            self.dist_aux = np.loadtxt(filename, comments='#',
                    usecols=range(a_dist_aux, a_dist_aux+self.N_dist_aux)).T
        
        self.dist_int = None
        if self.N_dist_int > 0:
            self.dist_int = np.loadtxt(filename, comments='#',
                    usecols=range(a_dist_int, a_dist_int+self.N_dist_int)).T


    @property
    def dt(self):
        tdiff = np.zeros(self.t.shape)
        tdiff[1:] = self.t[1:] - self.t[:-1]
        return tdiff

    def _readHeader(self, filename):

        hdr = {}

        with open(filename, 'r') as f:
            for line in f:
                if line[0] != '#':
                    break
                words = line.split()
                key = words[1]
                try:
                    val = int(words[2])
                except ValueError:
                    try:
                        val = float(words[2])
                    except ValueError:
                        val = words[2]
                hdr[key] = val

        return hdr


class DiscoReportTestLive(DiscoReport):

    def __init__(self, filename):
        super().__init__(filename)
        self.stride = self.NUM_PL_KIN + self.NUM_PL_AUX

    @property
    def Uf(self):
        return self.dist_int[:, :]

    @property
    def KIN_M(self):
        return self.shared[0::self.stride, :]
    @property
    def KIN_R(self):
        return self.shared[1::self.stride, :]
    @property
    def KIN_PHI(self):
        return self.shared[2::self.stride, :]
    @property
    def KIN_Z(self):
        return self.shared[3::self.stride, :]
    @property
    def KIN_PR(self):
        return self.shared[4::self.stride, :]
    @property
    def KIN_LL(self):
        return self.shared[5::self.stride, :]
    @property
    def KIN_PZ(self):
        return self.shared[6::self.stride, :]
    @property
    def KIN_SZ(self):
        return self.shared[7::self.stride, :]
    @property
    def KIN_EINT(self):
        return self.shared[8::self.stride, :]

    @property
    def SNK_M(self):
        return self.shared[0+self.NUM_PL_KIN::self.stride, :]
    @property
    def GRV_PX(self):
        return self.shared[1+self.NUM_PL_KIN::self.stride, :]
    @property
    def GRV_PY(self):
        return self.shared[2+self.NUM_PL_KIN::self.stride, :]
    @property
    def GRV_PZ(self):
        return self.shared[3+self.NUM_PL_KIN::self.stride, :]
    @property
    def GRV_JZ(self):
        return self.shared[4+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_PX(self):
        return self.shared[5+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_PY(self):
        return self.shared[6+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_PZ(self):
        return self.shared[7+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_JZ(self):
        return self.shared[8+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_SZ(self):
        return self.shared[9+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_MX(self):
        return self.shared[10+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_MY(self):
        return self.shared[11+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_MZ(self):
        return self.shared[12+self.NUM_PL_KIN::self.stride, :]
    @property
    def GRV_EGAS(self):
        return self.shared[13+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_EGAS(self):
        return self.shared[14+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_UGAS(self):
        return self.shared[15+self.NUM_PL_KIN::self.stride, :]
    @property
    def GRV_LZ(self):
        return self.shared[16+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_LZ(self):
        return self.shared[17+self.NUM_PL_KIN::self.stride, :]
    @property
    def GRV_K(self):
        return self.shared[18+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_K(self):
        return self.shared[19+self.NUM_PL_KIN::self.stride, :]
    @property
    def GRV_U(self):
        return self.shared[20+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_U(self):
        return self.shared[21+self.NUM_PL_KIN::self.stride, :]
    @property
    def SNK_EINT(self):
        return self.shared[22+self.NUM_PL_KIN::self.stride, :]
    @property
    def EXT_PX(self):
        return self.shared[23+self.NUM_PL_KIN::self.stride, :]
    @property
    def EXT_PY(self):
        return self.shared[24+self.NUM_PL_KIN::self.stride, :]
    @property
    def EXT_PZ(self):
        return self.shared[25+self.NUM_PL_KIN::self.stride, :]
    @property
    def EXT_JZ(self):
        return self.shared[26+self.NUM_PL_KIN::self.stride, :]
    @property
    def EXT_K(self):
        return self.shared[27+self.NUM_PL_KIN::self.stride, :]
    @property
    def EXT_U(self):
        return self.shared[28+self.NUM_PL_KIN::self.stride, :]


def interpolate_timeseries(t1, x1, t2):

    sh1 = x1.shape

    N1 = len(t1)
    N2 = len(t2)

    if len(sh1) == 1:
        sh2 = (N2,) 
    else:
        sh2 = (N2, sh1[1:])

    x2 = np.empty(sh2)

    idx = np.searchsorted(t1, t2)

    for i in range(N2):
        if idx[i] == 0:
            x2[i] = x1[0]
        elif idx[i] >= N1:
            x2[i] = x1[N1-1]
        else:
            idxa = idx[i] - 1
            idxb = idx[i]
            wa = (t1[idxb] - t2[i]) / (t1[idxb] - t1[idxa])
            wb = (t2[i] - t1[idxa]) / (t1[idxb] - t1[idxa])
            x2[i, ...] = wa * x1[idxa, ...] + wb * x1[idxb, ...]

    return x2
