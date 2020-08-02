#!/usr/bin/env python

import os, sys
from dockml import pdbIO


class GaussianInput :

    def __init__(self, pdbin):
        PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
        self.inpdb = pdbin
        self.sample_in = PROJECT_ROOT + "/../data/sample_gaussian.inp"

    def parsePDB(self):

        pio = pdbIO.parsePDB(inPDB=self.inpdb)
        infor = pio.atomInformation(self.inpdb)

        atomndx = infor.keys()
        elements= [ infor[x][7] for x in atomndx ]

        return (atomndx, elements)

    def getAtomCrds(self):

        atomndx = [ x for x in self.parsePDB()[0] ]

        pio = pdbIO.coordinatesPDB()
        coords = pio.getAtomCrdByNdx(self.inpdb, atomndx)

        return coords

    def pdb2GauInput(self, output, charge, spin):

        with open(self.sample_in) as lines :
            header = lines.readlines()

        with open(output, 'w') as tofile :
            for s in header :
                tofile.write(s)
            tofile.write(" \n")

            tofile.write("%d %d \n"%(charge, spin))

            ndx, element = self.parsePDB()
            coords  = self.getAtomCrds()

            for i in range(len(coords)) :
                tofile.write("%s   %8.4f  %8.4f  %8.4f  \n"%(element[i],
                                                             coords[i][0],
                                                             coords[i][1],
                                                             coords[i][2])
                             )

            tofile.write(" \n")

        return 1

if __name__ == "__main__" :

    if len(sys.argv) < 3 :
        print("usage: gaussian.py input.pdb gau.inp charge spin")
        sys.exit(1)

    inpdb = sys.argv[1]
    output= sys.argv[2]

    if len(sys.argv) < 3 :
        charge = 0
        spin   = 1
    else :
        charge= int(sys.argv[3])
        spin  = int(sys.argv[4])

    gau = GaussianInput(inpdb)
    gau.pdb2GauInput(output, charge, spin)
