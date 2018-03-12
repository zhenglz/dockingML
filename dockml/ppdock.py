# -*- coding: utf-8 -*-

#!/usr/bin/env python

from dockml import pdbIO
from dockml import algorithms
import sys, os
import argparse
import numpy as np
from argparse import RawTextHelpFormatter

class DockSurfRes :

    def __init__(self):
        pass

    def surfaceRes(self, reslist, recpdb, chain=["A",], atoms=["CA"]):
        """

        :param reslist: list, a list of intergers of residue sequence
               which consist of the pocket
        :param recpdb:str, pdb file name
        :param chain: list of str, chains in consider
        :param atoms: list of str, atom names in consider
        :return: list of list, M * 3, points
        """

        points = []
        with open(recpdb) as lines :
            key_lines = [ x for x in lines
                          if ("ATOM" in x and
                              x.split()[2] in atoms and
                              x[21] in chain and
                              int(x[22:26].strip()) in reslist)
                          ]

            coord = pdbIO.coordinatesPDB()
            points = coord.getAtomCrdFromLines(key_lines)

        return points

    def distanceToPlane(self, recpdb, plane, chain=["A"], atoms=["CA"]):

        pfit = algorithms.PlaneFit()

        with open(recpdb) as lines :
            lines = [ x for x in lines
                      if ("ATOM" in x and
                          x[21] in chain and
                          x.split()[2] in atoms
                          )
                      ]

            coord = pdbIO.coordinatesPDB()
            points = coord.getAtomCrdFromLines(lines)

            distances = []
            for p in points :
                distances.append(pfit.point_distance(params=plane, point=p))

        return list(zip(distances, lines))

def main() :

    d = '''
    Prepare protein-protein docking non-used residues
    '''
    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)

    parser.add_argument('-pdb',type=str, default="receptor.pdb",
                        help="Reference PDB file of the receptor \n")
    parser.add_argument('-surf_res', type=str, nargs="+", default=['1', '2', '3'],
                        help="Surface residue sequence index. \n"
                             "A file or a list of residue index. \n")
    parser.add_argument('-out', type=str, default='non_used_res_list.dat',
                        help="Output not used residue list. \n")
    parser.add_argument('-dist_cutoff', type=float, default=5.0,
                        help="Distance cutoff for non_used residues. \n"
                             "Default is 5.0 Angstrom. \n")
    parser.add_argument('-opdb', type=str, default='non_used_receptor.pdb',
                        help="Not used pdb file content. \n")

    args = parser.parse_args()

    if len(sys.argv) < 2 :
        parser.print_help()

    if os.path.exists(args.surf_res[0]) :
        reslist = []
        with open(args.surf_res[0]) as lines :
            reslist = [ int(x.split()[0]) for x in lines if x[0] != "#" ]
    else :
        reslist = [ int(x) for x in args.surf_res ]

    print(reslist)

    ppdock = DockSurfRes()
    key_points = ppdock.surfaceRes(reslist, args.pdb)
    print(key_points)

    pfit = algorithms.PlaneFit()
    plane = pfit.fitPlane(key_points)
    print(plane)

    dist_lines = ppdock.distanceToPlane(args.pdb, plane)

    opdb = open(args.opdb, 'w')
    nonr = open(args.out,  'w')

    for e in dist_lines :
        d = e[0].tolist()[0][0]
        print(d)

        if (args.dist_cutoff < 0 and d < args.dist_cutoff) or (args.dist_cutoff > 0 and d > args.dist_cutoff):
            opdb.write(e[1])
            nonr.write("%s \n"%(e[1][22:26].strip()))

    opdb.close()
    nonr.close()

    print("Prepare Zdock file completed!")