# -*- coding: utf-8 -*-

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

        return deltaZ, middle

class areaPerLip :

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

    def selectProteinAtomsCrds(self, pdb, zrange=[0, 1.0]):

        pdbio = dockml.parsePDB()
        atominfor = pdbio.atomInformation(pdb)

        ndx = atominfor.keys()

        # a list of atom index, strings
        selected_ndx = [ x for x in ndx if atominfor[x][1] == "Protein" ]

        cpdb = dockml.coordinatesPDB()
        crds = np.asarray(cpdb.getAtomCrdByNdx(pdb, selected_ndx))

        selected_crds = [ x for x in crds if ((x[2] > zrange[0]) and (x[2] < zrange[1]))]

        return np.asarray(selected_crds)
