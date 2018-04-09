import dockml, mdanaly
import numpy as np
import os, sys, math
from mpi4py import MPI

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
