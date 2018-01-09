# -*- coding: utf-8 -*-

import sys
import dockml
import numpy as np
from scipy.spatial import ConvexHull
import argparse
from argparse import RawTextHelpFormatter

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
        '''
        number of lipids in upper and lower leaflet
        :param zvalues:
        :param middle:
        :return:
        '''

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
        '''
        get protein atoms' coordinates if the protein atoms z values in zrange
        :param zrange:
        :return:
        '''

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
        '''
        get PBC vectors, and return the total area of the bilayer lipids
        :param vectors:
        :return:
        '''

        return vectors[0] * vectors[1]

def arguments() :

    d = '''
    Descriptions of the lipid Thickness and Lipid Per Area calculation
    
    '''
    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)

    parser.add_argument('-pdbf', type=str, default="S_%d.pdb",
                        help="Format of pdb files. For example, \n"
                             "you have S_1.pdb to S_100.pdb, you could\n"
                             "set this argument as \'S_%d.pdb\' \n"
                             "Default is S_%d.pdb")
    parser.add_argument('-out', type=str, default="thickness_APL.dat",
                        help="Output file name. Default is thickness_APL.dat \n")
    parser.add_argument('-num', type=int, default=1,
                        help="Number of files for analysis. Default is 1 .\n")
    parser.add_argument('-layer', type=float, default=2.5,
                        help="To select a layer of protein atoms, you set a\n"
                             "thickness for the layer. Default is 2.5 Angstrom.\n")
    parser.add_argument('-res', type=str, nargs="+", default=['POC'],
                        help="Lipid residues for analysis. \n"
                             "Default is [ POC ]. \n")
    parser.add_argument('-head', type=str, nargs='+', default=['P'],
                        help="Head atoms for lipid thickness and APL calculations.\n"
                             "Default is [ P ]. \n")
    parser.add_argument('-pbc', type=float, nargs="+", default=[0.0, 0.0, 0.0 ],
                        help="PBC vectors (3 elements) for analysis. Unit is Angstrom. \n")
    parser.add_argument('-grid', type=int, default=20,
                        help="Grid number for P z coordinates analysis. \n"
                             "Default is 20. \n")

    args = parser.parse_args()

    return args

def main() :

    args = arguments()

    N = args.num
    layer_step = args.layer
    ntype = args.pdbf

    # extract all files
    files = [ ntype.split("_")[0]+"_"+str(x)+".pdb" for x in range(1, N+1) ]

    tofile = open(args.out, 'w')
    tofile.write("#Thickness Average_APL Upper_APL Lower_APL \n")

    for f in files :
        print("Frame %s " % (f) )
        lip = lipidThickness(f, args.res, args.head)
        zv = lip.getZvalues()

        thick, middle, up, low = lip.deltaZcoord(zv, args.grid)

        up_num_lip, low_num_lip = lip.lipidsNum(zv, middle)

        up_zrange = [ up - layer_step, up + layer_step]
        low_zrange= [ low - layer_step, low+ layer_step]

        apl = areaPerLipid(f)
        up_proarea = apl.proteinArea(apl.selectProteinAtomsCrds(up_zrange))
        low_proarea = apl.proteinArea(apl.selectProteinAtomsCrds(low_zrange))

        total_area = apl.totalArea(vectors=args.pbc)

        up_alp = (total_area - up_proarea )  / float(up_num_lip)
        low_alp= (total_area - low_proarea) / float(low_num_lip)

        alp_aver = (total_area * 2 - up_proarea - low_proarea) / float(up_num_lip+low_num_lip)

        print("%5d %5d %8.3f %8.3f "%(up_num_lip, low_num_lip, up_proarea, low_proarea))
        print("%8.3f %8.3f %8.3f %8.3f"%(thick, up_alp, low_alp, alp_aver))

        tofile.write("%6.3f %6.3f %6.3f %6.3f \n"%(thick, alp_aver, up_alp, low_alp))

    tofile.close()

    print("Calculations of Thickness and APL completed!")