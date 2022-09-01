### Remaps, shifting to the co-rotating frame ###

import sys
import numpy as np
import discopy as dp

def remap(d1, d2):

    z = d2.Z
    r = d2.R
    nq = d1.numQ
    pi2 = 2.0*np.pi

    prim = []
    for k in range(d2.numZ):
        sheet = []
        for j in range(d2.numR):
            sheet.append(np.empty((d2.numPhi[k,j], nq),np.float))
        prim.append(sheet)

    for k in range(d2.numZ):
        for j in range(d2.numR):
            print("   Processing k = {0:04d}, j = {1:04d}".format(k, j))
            phi = d2.Phi(k, j)
            for i in range(d2.numPhi[k,j]):
                prim1 = d1.getPrimAt(r[j], phi[i]+(d1.t%pi2), z[k])
                prim[k][j][i,:] = prim1[:]

    d2.t = d1.t
    d2.prim = prim
    d2._opts = d1._opts.copy()
    d2.planets = d1.planets.copy()



if __name__ == "__main__":

    if len(sys.argv) < 3:
        print("usage: $ python remapCheckpoint.py <checkpoint.h5> <in.par> <filename>")
        print("Remaps checkpoint to grid given by in.par. Saves to file <filename>")
        sys.exit()

    chkin = sys.argv[1]
    parfile = sys.argv[2]
    chkout = sys.argv[3]

    d1 = dp.DiscoDomain(chkin)
    d2 = dp.DiscoDomain(parFile=parfile)
    remap(d1, d2)
    d2.writeCheckpoint(chkout)
