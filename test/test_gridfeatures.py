import dockml, mdanaly
import numpy as np
import os, sys, math
from mpi4py import MPI
from datetime import datetime

if __name__ == "__main__" :

    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.rank

    oldtime = datetime.now()

    inp = sys.argv[1]
    out = sys.argv[2]

    with open(inp) as lines :
        filesList = [ x.split()[0] for x in lines if ("#" not in x and len(x.split())) ]

    if rank == 0:
        load4each = int(math.ceil(float(len(filesList)) / float(size)))
        rankfilesList = []

        for i in range(size - 1):
            rankfilesList.append(filesList[i * load4each: load4each * (i + 1)])
        rankfilesList.append(filesList[(size - 1) * load4each:])
        print(rankfilesList)
    else:
        rankfilesList = None

    ## scatter data to sub-processers/threads
    rankfilesList = comm.scatter(rankfilesList, root=0)

    features_out = []

    #print(rank, rankfilesList)

    gf = dockml.GridBasedFeature(rankfilesList[0], rankfilesList[0], 1.2)
    ligcoords = gf.ligCoords(" ")

    if not os.path.exists("GRID"):
        print(list(ligcoords.values()))
        pocket = gf.createCubicPocket(coordinates=list(ligcoords.values()), extension=0)
        grids  = gf.generateGrids(pocket)
    else:
        grids = np.loadtxt("GRID")

    features_out = []
    for ligpdb in rankfilesList :
        print("Rank %d Initialize !"%rank)
        #ligpdb = 'M1.pdbqt'

        '''with open(ligpdb) as lines :
            clines = [ x for x in lines if "ATOM" in x ]
            ligcoords = dockml.coordinatesPDB().getAtomCrdFromLines(clines) '''

        gf = dockml.GridBasedFeature(ligpdb, ligpdb, 1.0)
        ligcoords = gf.ligCoords(" ")

        print("Rank %d GRID Loading completed!" % rank)

        gmap = gf.gridAtomMap(grids, ligcoords)
        properties = gf.atomProperties(gf.ligndx)
        features = gf.gridBinProperty(gmap, properties, count_cutoff=0.2)

        print("Rank %d Gather data correctly! "%rank)

        print(features)

        feature_vector = np.reshape(np.array(list(features.values())), (1, -1))[0]
        print(len(feature_vector))

        features_out.append(feature_vector)

        print(np.array(features_out).shape)
        #print(np.sum(np.array(features.values()), axis=0))

        print("Rank %d %20s completed " % (rank, ligpdb))
        delta = datetime.now() - oldtime
        print("Rank %d Progress: Time Usage for 1 frame %d seconds" % (rank, delta.total_seconds()))
        oldtime = datetime.now()

    if rank == 0 :
        #print(features_out)
        overallFeatures = comm.gather(features_out, root=0)[0]

        fsize = int(len(overallFeatures) / 5)

        HEAD = []

        for i in range(fsize) :
            fs = [ 'mass', 'vdwe', 'vdws', 'charge', 'negativity']
            HEAD += [ x + "_" + str(i) for x in fs ]
        HEAD = ",".join(HEAD)

        print(np.asarray(overallFeatures).shape)
        np.savetxt(out, overallFeatures, delimiter=',', fmt="%12.6f", header=HEAD)

        print("Completed!")
    MPI.Finalize()
    sys.exit(1)
