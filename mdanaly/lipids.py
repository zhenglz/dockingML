# -*- coding: utf-8 -*-

import os, sys
from dockml.pdbIO import coordinatesPDB, handlePBC, parsePDB
import numpy as np
from scipy.spatial import ConvexHull
import argparse
from argparse import RawTextHelpFormatter


class LipidThickness(object):
    """Lipid Thickness calculation and lipid surface calculation.

    Parameters
    ----------
    pdb : str,
        The input pdb file name
    lipRes : list
        The list of lipid residue names in the pdb file
    headatoms : list
        The list of head atom name list in the pdb file

    Attributes
    ----------
    pdb
    lip
    head

    Notes
    -----
    The surface area of the protein around the upper or lower
    leaflets, could be estimated with convex hull algorithm.

    References
    ----------
    Detail method could be found here:
    http://scipy-cookbook.readthedocs.io/items/Finding_Convex_Hull.html
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html

    """

    def __init__(self, pdb, lipRes=['DOPC'], headatoms=["P"]):

        self.pdb = pdb
        self.lip = lipRes
        self.head = headatoms

    def getZvalues(self):
        """
        Get all z coordinates of selected residues and their head group atoms.

        Returns
        -------
        zvalues: numpy ndarray, shape=[N, 1]
            the z-coordinates of lines, N is number of atoms
        """

        if not os.path.exists(self.pdb):
            print("PDB file %s not exist. Exit now! " % self.pdb)
            sys.exit(0)

        with open(self.pdb)as lines:
            lines = [x for x in lines if "ATOM" in x]
            plines = [x for x in lines if ((x[17:20] in self.lip) and (x.split()[2] in self.head))]

            coord = coordinatesPDB()
            zvalues = np.asarray(coord.getAtomCrdFromLines(plines))[:, 2]

        return zvalues

    def deltaZcoord(self, zvalues, numbins=20):
        """
        Find the up leaflet and low leaflet average Z values.

        Parameters
        ----------
        zvalues: numpy ndarray,
            the z-coordinates of a group of atoms
        numbins: int,
            number of bins for histogram

        Returns
        -------
        deltaZ: float,
            the average distance between upper leaflet and lower leaflet
        middle: float,
            the middle layer z-coordinates
        up_aver: float
            the average z-coordinate of the upper leaflet
        low_aver: float
            the average z-coordinate of the lower leaflet

        """

        hist, edge = np.histogram(zvalues, bins=int(numbins))

        middle = edge[int(numbins/2)]

        upleaflet = [x for x in zvalues if x > middle]
        lowleaflet= [x for x in zvalues if x < middle]

        up_aver = np.mean(upleaflet)
        low_aver = np.mean(lowleaflet)

        deltaZ = up_aver - low_aver

        return deltaZ, middle, up_aver, low_aver

    def lipidsNum(self, zvalues, middle):
        """Number of lipids in upper and lower leaflet

        Parameters
        ----------
        zvalues : np.ndarray, or list, shape = [ N, ]
            the Z-coordinates
        middle : np.ndarray, or list, shape = [ N, ]
            the middle coordinates

        Returns
        -------
        num_upleaflet : int
        num_lowleaflet : int
        """

        upleaflet = [x for x in zvalues if x > middle]
        lowleaflet = [x for x in zvalues if x < middle]

        return len(upleaflet), len(lowleaflet)


class AreaPerLipid(object):

    def __init__(self, pdb):
        self.pdb = pdb

    def proteinArea(self, proCrds, vectors, restorePBC=False):
        """Calculate protein area given a set of 2d points (xy data only)
        with convex hull algorithm.

        Detail method could be found here:
        http://scipy-cookbook.readthedocs.io/items/Finding_Convex_Hull.html
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html

        Parameters
        ----------
        proCrds : np.ndarray
            a list of xyz coordinates, 3*N array
        vectors : list, or np.array, shape = [ 3, ]
            pbc xyz vector
        restorePBC : bool
            whether restore PBC box

        Returns
        -------
        area: float
            the area of the protein area, unit square angstrom

        """

        if restorePBC:
            pbc = handlePBC()
            pbcVec = [
                [0.0, vectors[0]],
                [0.0, vectors[1]],
                [0.0, vectors[2]],
            ]
            coords = pbc.crdPBCRestore(proCrds, pbcVec)
        else :
            coords = proCrds

        xy_coords = np.asarray(coords)[:, :2]

        if xy_coords.shape[0]:
            vdw_dummy = self.atomVdWBoundary(xy_coords)

            # convex hull algorithm to determine the area occupied by atoms
            hull = ConvexHull(np.concatenate((xy_coords, vdw_dummy)))
            # TODO: check the area calculation result
            area = hull.area

            return area
        else:
            return 0.0

    def selectProteinAtomsCrds(self, zrange=[0, 1.0]):
        """get protein atoms' coordinates if the protein
        atoms z coordinates in zrange

        Parameters
        ----------
        zrange : iterable, length = 2
            up and low boundaries

        Returns
        -------
        selected_atoms : np.ndarray, shape = [ N, 3]
            the xyz coordinates of selected atoms
            N is number of atoms selected
        """

        pdbio = parsePDB()
        atominfor = pdbio.atomInformation(self.pdb)

        ndx = atominfor.keys()

        # a list of atom index, strings
        selected_ndx = [x for x in ndx
                        if atominfor[x][1] == "Protein"]

        cpdb = coordinatesPDB()
        crds = np.asarray(cpdb.getAtomCrdByNdx(self.pdb, selected_ndx))

        selected_crds = [x for x in crds
                         if ((x[2] > zrange[0]) and (x[2] < zrange[1]))]

        return np.asarray(selected_crds)

    def atomVdWBoundary(self, points, vcutoff=1.4):
        """given a list of points, for each of points,
        find it's 4 neihboring dummy points
        to represent its vdw boundary
                   dummy2
                     |
          dummy1-- atom --  dummy3
                     |
                   dummy4

        Parameters
        ----------
        points : np.ndarray, shape = [ N, 3]
            the coordinates of a list of atoms
            N is the number of atoms
        vcutoff : float
            the van del waals distance cutoff

        Returns
        -------
        dummy_atoms : np.ndarray, shape = [ 4*N, 2]
            the coordinates of 4*N dummy atoms
            N is the number of real atoms
        """

        dummy = []
        for p in points:
            d_p = map(lambda x, y: [x + vcutoff*p[0], y+vcutoff*p[1]],
                      [-1.0, 1.0, -1.0, 1.0], [-1, 1.0, 1.0, -1.0]
                      )
            dummy += list(d_p)

        return np.asarray(dummy)

    def totalArea(self, vectors):
        """get PBC vectors, and return the total area of
        the bilayer lipids total area = x_length * y_length

        Parameters
        ----------
        vectors : iterable, length = 2
            the vectors of x_length and y_length

        Returns
        -------
        area : float
            the total areal of lipids

        """

        return vectors[0] * vectors[1]


def arguments():

    d = '''
    Descriptions of the lipid Thickness and Lipid Per Area calculation
    
    '''
    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-pdbf', type=str, default="S_%d.pdb",
                        help="Format of pdb files. For example, \n"
                             "you have S_1.pdb to S_100.pdb, you could\n"
                             "set this argument as S_d.pdb \n"
                             "Default is S_d.pdb")
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
    parser.add_argument('-pbc', type=float, nargs="+", default=[ ],
                        help="PBC vectors (3 elements) for analysis. Unit is Angstrom. \n")
    parser.add_argument('-grid', type=int, default=20,
                        help="Grid number for P z coordinates analysis. \n"
                             "Default is 20. \n")
    parser.add_argument('-restore', type=bool, default=False,
                        help="Whether restore the protein PBC given specific PBC conditions. \n"
                             "Default is False. ")

    args = parser.parse_args()

    return args

def main():

    args = arguments()

    N = args.num
    layer_step = args.layer

    # extract all files
    files = [ args.pdbf%x for x in range(1, N+1) ]

    tofile = open(args.out, 'w')
    tofile.write("#Thickness Average_APL Upper_APL Lower_APL \n")

    for f in files :
        print("Frame %s " % (f) )

        apl = AreaPerLipid(f)
        lip = LipidThickness(f, args.res, args.head)

        zv = lip.getZvalues()

        # get pbc from files
        if len(args.pbc) :
            pbc = args.pbc

        else :
            pbcpdb = dockml.handlePBC()
            pbc = pbcpdb.getPBCFromPBD(f)
            pbc = list(np.asarray(pbc)[:, 1])

        total_area = apl.totalArea(vectors=pbc)

        thick, middle, up, low = lip.deltaZcoord(zv, args.grid)
        up_num_lip, low_num_lip = lip.lipidsNum(zv, middle)

        up_zrange = [up - layer_step, up + layer_step]
        low_zrange = [low - layer_step, low+ layer_step]

        # calculate protein atom area, restore PBC if required
        up_proarea = apl.proteinArea(apl.selectProteinAtomsCrds(up_zrange), pbc, restorePBC=args.restore)
        low_proarea = apl.proteinArea(apl.selectProteinAtomsCrds(low_zrange), pbc, restorePBC=args.restore)

        up_alp = (total_area - up_proarea) / float(up_num_lip)
        low_alp = (total_area - low_proarea) / float(low_num_lip)

        alp_aver = (total_area * 2 - up_proarea - low_proarea) / float(up_num_lip+low_num_lip)

        print("%5d %5d %8.3f %8.3f "%(up_num_lip, low_num_lip, up_proarea, low_proarea))
        print("%8.3f %8.3f %8.3f %8.3f"%(thick, up_alp, low_alp, alp_aver))

        tofile.write("%6.3f %6.3f %6.3f %6.3f \n"%(thick, alp_aver, up_alp, low_alp))

    tofile.close()

    print("Calculations of Thickness and APL completed!")

