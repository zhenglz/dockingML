# -*- coding: utf-8 -*-

import sys
import dockml
import numpy as np
from scipy.spatial import ConvexHull

class lipidThickness :

    def __init__(self, pdb, lipRes=['DOPC'], headatoms=["P"]):

        self.pdb = pdb
        self.lip = lipRes
        self.head= headatoms

    def getZvalues(self):
        '''
        get all z coordinates of selected residues and their head group atoms
        :return:
        '''

        with open(self.pdb)as lines :
            lines = [ x for x in lines if "ATOM" in x ]
            plines = [ x for x in lines if ( (x[17:20] in self.lip) and (x.split()[2] in self.head) ) ]

            coord = dockml.coordinatesPDB()
            zvalues = np.asarray(coord.getAtomCrdFromLines(plines))[:, 2]

        return zvalues

    def deltaZcoord(self, zvalues, numbins=20):
        '''
        find the up leaflet and low leaflet average Z values
        :param zvalues:
        :param numbins:
        :return:
        '''

        hist, edge = np.histogram(zvalues, bins=int(numbins))

        middle = edge[int(numbins/2)]

        upleaflet = [ x for x in zvalues if x > middle ]
        lowleaflet= [ x for x in zvalues if x < middle ]

        up_aver = np.mean(upleaflet)
        low_aver = np.mean(lowleaflet)

        deltaZ = up_aver - low_aver

        return deltaZ, middle, up_aver, low_aver

    def lipidsNum(self, zvalues, middle):

        upleaflet = [x for x in zvalues if x > middle]
        lowleaflet = [x for x in zvalues if x < middle]

        return len(upleaflet), len(lowleaflet)

class areaPerLipid :

    def __init__(self, pdb):
        self.pdb = pdb

    def proteinArea(self, proCrds):
        '''
        calculate protein area given a set of 2d points (xy data only)
        convex hull algorithm
        http://scipy-cookbook.readthedocs.io/items/Finding_Convex_Hull.html
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html
        :param proCrds: a list of xyz coordinates, 3*N array
        :return: area, float
        '''

        xy_coords = np.asarray(proCrds)[:, :2]

        hull = ConvexHull(xy_coords)
        area = hull.area

        return area

    def selectProteinAtomsCrds(self, zrange=[0, 1.0]):

        pdbio = dockml.parsePDB()
        atominfor = pdbio.atomInformation(self.pdb)

        ndx = atominfor.keys()

        # a list of atom index, strings
        selected_ndx = [ x for x in ndx if atominfor[x][1] == "Protein" ]

        cpdb = dockml.coordinatesPDB()
        crds = np.asarray(cpdb.getAtomCrdByNdx(self.pdb, selected_ndx))

        selected_crds = [ x for x in crds if ((x[2] > zrange[0]) and (x[2] < zrange[1]))]

        return np.asarray(selected_crds)

    def totalArea(self, vectors):

        return vectors[0] * vectors[1]

def main() :

    N=5
    layer_step = 0.5

    # extract all files
    files = [ "S_%d.pdb"%x for x in range(1, N) ]

    for f in files :
        print("Frame %s " % (f) )
        lip = lipidThickness(f, ["POP"], ["P"])
        zv = lip.getZvalues()

        thick, middle, up, low = lip.deltaZcoord(zv, 20)

        up_num_lip, low_num_lip = lip.lipidsNum(zv, middle)

        up_zrange = [ up - layer_step, up + layer_step]
        low_zrange= [ low - layer_step, low+ layer_step]

        apl = areaPerLipid(f)
        up_proarea = apl.proteinArea(apl.selectProteinAtomsCrds(up_zrange))
        low_proarea = apl.proteinArea(apl.selectProteinAtomsCrds(low_zrange))

        total_area = apl.totalArea(vectors=[91.1,91.1,97.9])

        up_alp = (total_area - up_proarea )  / float(up_num_lip)
        low_alp= (total_area - low_proarea) / float(low_num_lip)

        alp_aver = (total_area * 2 - up_proarea - low_proarea) / float(up_num_lip+low_num_lip)

        print("%5d %5d %8.3f %8.3f "%(up_num_lip, low_num_lip, up_proarea, low_proarea))
        print("%8.3f %8.3f %8.3f %8.3f"%(thick, up_alp, low_alp, alp_aver))