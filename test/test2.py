import dockml, mdanaly
import numpy as np
import os, sys, math
from mpi4py import MPI

def icriticalMap(icritical, dat) :

    map = np.sum( np.greater(dat, icritical), axis=0 ) / dat.shape[0]

    return map

def diagonalZeroed(cmap, xsize, outfile) :
    from mdanaly import matrix

    map = np.reshape(cmap, (xsize, cmap.shape[0]/xsize ))

    # matrix 2 xyz
    xyz = matrix.MatrixHandle().matrix2xyz(map)
    mtx = matrix.MatrixHandle().neiborhood2zero(xyz, neiborsize=4, outtype='mtx')

    np.savetxt(outfile, mtx, fmt="%3.1f", delimiter=" ")
    return 1

def generateCmap() :

    dat = np.loadtxt("all_cmapnbyn.csv", delimiter=",")
    Icritical = np.arange(0, 15, 1.0)

    for ic in Icritical :
        imap = icriticalMap(ic, dat)

        diagonalZeroed(imap, int(np.sqrt(dat.shape[1])), "cmap_I_%f"%ic)

if __name__ == "__main__" :

    ligpdb = 'M1.pdbqt'
    features_out = []

    '''with open(ligpdb) as lines :
        clines = [ x for x in lines if "ATOM" in x ]
        ligcoords = dockml.coordinatesPDB().getAtomCrdFromLines(clines) '''

    gf = dockml.GridBasedFeature(ligpdb, ligpdb, 0.8)
    ligcoords, ndx = gf.ligCoords(" ", )

    if not os.path.exists("GRID") :
        pocket = gf.createCubicPocket(ligcoords, extension=4)
        grids = gf.generateGrids(pocket)
    else :
        grids = np.loadtxt("GRID")

    gmap = gf.gridAtomMap(grids, ndx, ligcoords, count_cutoff=0.2)
    features = gf.gridBinProperty(gmap, )

    all_features = []

    for key in features.keys() :
        all_features += features[key]

    #print(all_features)
    #print(np.sum(np.array(features.values()), axis=0))

    print("%s %20s completed " % (0, ligpdb))
