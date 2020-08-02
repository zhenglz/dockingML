#!/usr/bin/env python

import sys
import mdtraj as mt
import numpy as np

class RestraintProcess(object):

    def __init__(self, pdbfile, fc=100.0, start_id=1, start_res=1):
        self.pdb = pdbfile

        self.top = None
        self.k = fc
        self.start_atom = start_id
        self.start_res = start_res

    def pdb_parser(self):

        pdb = mt.load_pdb(self.pdb)
        self.top = pdb.topology

    def atomid(self, atomname, resid):

        ids =  self.top.select("name %s and residue %d" % (atomname, int(resid)+self.start_res))

        if len(ids):
            return ids[0] + self.start_atom
        else:
            return None

    def read_dihedral_rst(self, line):

        atomnames = [line.split()[1], line.split()[3], line.split()[5], line.split()[7], ]
        resseq    = [line.split()[2], line.split()[4], line.split()[6], line.split()[8], ]

        atom_ids = []
        for i in range(4):
            id = self.atomid(atomnames[i], resseq[i])
            if id is not None:
                atom_ids.append(id)
            else:
                pass

        angle = float(line.split()[-2])
        # radian to degree
        angle = 180 * angle / 3.1415928

        if len(atom_ids) == 4:
            return "  %4d    %4d    %4d   %4d     1     1  %.2f     0     %.1f    2 \n" \
                   % (atom_ids[0], atom_ids[1], atom_ids[2], atom_ids[3], angle, self.k)
        else:
            return ""

    def read_dist_rst(self, line, fc=100.0):

        atomnames = [line.split()[1], line.split()[3]]
        resseq    = [line.split()[2], line.split()[4]]

        atom_ids = []
        for i in range(2):
            id = self.atomid(atomnames[i], resseq[i])

            if id is not None:
                atom_ids.append(id)
            else:
                pass

        # transform angstrom to nanometer
        distance = float(line.split()[-2]) / 10.0

        if len(atom_ids) == 2:
            info = "%4d  %4d  1  0  1  0.2  0.4  0.6  %.1f  \n" % (atom_ids[0], atom_ids[1], fc)
            return info
        else:
            return ""


if __name__ == "__main__":

    # restraint file name
    fn_rest = sys.argv[1]

    # reference pdb file name
    fn_empdb = sys.argv[2]

    out = sys.argv[3]

    rest = RestraintProcess(fn_empdb, )
    rest.pdb_parser()

    tofile = open(out, 'w')

    # for dihedral
    #tofile.write("[ dihedral_restraints ]  \n; ai   aj    ak    al  type  label  phi  dphi  kfac  power\n")

    # for distance
    tofile.write("[ distance_restraints ]  \n; ai aj type index type’ low up1 up2 fac \n")

    with open(fn_rest) as lines:
        for s in lines:
            #constr = rest.read_dihedral_rst(s)
            constr = rest.read_dist_rst(s, 100)
            if len(constr):
                tofile.write(constr)

    tofile.close()


