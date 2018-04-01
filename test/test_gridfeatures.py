import dockml, mdanaly
import numpy as np


ligpdb = 'M1.pdbqt'

'''with open(ligpdb) as lines :
    clines = [ x for x in lines if "ATOM" in x ]
    ligcoords = dockml.coordinatesPDB().getAtomCrdFromLines(clines) '''



gf = dockml.GridBasedFeature(ligpdb, ligpdb, )

ligcoords, ndx = gf.ligCoords(" ", )
pocket = gf.createCubicPocket(ligcoords, extension=4)

print(pocket)

grids = gf.generateGrids(pocket, gridsize=0.8)

print(grids.shape)

gmap = gf.gridAtomMap(grids, ndx, ligcoords, count_cutoff=0.2)

features = gf.gridBinProperty(gmap, )
print(features)
#print(np.sum(np.array(features.values()), axis=0))