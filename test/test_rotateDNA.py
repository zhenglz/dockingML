#!/usr/bin/env python

import argparse
import math
import os
import sys

import numpy as np

import automd.shiftpdb as spdb
import dockml.algorithms as dal
import dockml.pdbIO as pio


#import automd.sumpdb as supdb

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def perpendicular_vector(v):
    if v[1] == 0 and v[2] == 0:
        if v[0] == 0:
            raise ValueError('zero vector')
        else:
            return np.cross(v, [0, 1, 0])
    return np.cross(v, [1, 0, 0])

def replaceXYZ(inpdb, outpdb, newXYZ) :

    with open(inpdb) as lines :
        lines = [ x for x in lines if "ATOM" in x ]

    with open(outpdb, 'w') as tofile :
        for i in range(len(lines)) :
            tofile.write(spdb.shiftPDB(inpdb).xyzReplace(lines[i], newXYZ[i]))

    tofile.close()

    return 1

def main() :
    d = '''
    Usage:
    python rotateDNA.py -inp input.pdb -out output.pdb -dchains C D 
    
    '''
    parser = argparse.ArgumentParser(description=d)

    parser.add_argument("-inp", default="input.pdb", type=str, help="Input PDB file for protein DNA complex")
    parser.add_argument("-out", default="output.pdb", type=str, help="Output PDB file for protein DNA complex")
    parser.add_argument("-dchains", default=['C', 'D'], type=str, nargs="+",
                        help="The chain ids of DNA molecules. ")
    parser.add_argument("-angle", type=int, default=90, help="The rotation angle of the DNA molecule. ")
    args = parser.parse_args()

    if len(sys.argv) < 2 :
        parser.print_help()
        sys.exit(0)

    # dna pdb file name
    fn = args.inp
    output = args.out
    compdb = "com_DNA.pdb"
    angle = np.pi * args.angle / 180.0  # rotation angle, pi/2 = 90 degree

    with open(fn) as lines :
        chains = set([ x[21] for x in lines if "ATOM" in x ])

    DNA_chains = args.dchains

    protein_chains = []
    if len(protein_chains) == 0 :
        protein_chains = list(chains - set(DNA_chains))

    # get DNA coordinates
    with open(fn) as lines :
        lines = [ x for x in lines if ("ATOM" in x and x[21] in DNA_chains) ]
        with open("DNA_tmp.pdb", 'w') as tofile :
            for x in lines :
                tofile.write(x)
        DNA_xyz = pio.coordinatesPDB().getAtomCrdFromLines(lines)

    # get center of mass of a DNA molecule
    com = np.array(DNA_xyz).mean(axis=0)
    trans_xyz = np.array(DNA_xyz) - com

    # create a center of mass shifted DNA and protein molecule
    replaceXYZ("DNA_tmp.pdb", compdb, trans_xyz)

    # create a transformed protein pdb
    with open(fn) as lines :
        lines = [ x for x in lines if ("ATOM" == x.split()[0] and x[21] in protein_chains) ]
        pro_xyz = pio.coordinatesPDB().getAtomCrdFromLines(lines)
        shifted_proxyz = np.array(pro_xyz) - com
        with open("PRO_tmp.pdb", 'w') as tofile :
            for x in lines :
                tofile.write(x)
        replaceXYZ("PRO_tmp.pdb", "transformed_protein.pdb", shifted_proxyz)

    # rotated_xyz, the coordinations of the transformed DNA
    axis = dal.LineFit(trans_xyz).fit_line()
    rotated_xyz  = []
    for p in trans_xyz :
        # p is a point, an atom in a DNA molecule
        rota_p = np.dot(rotation_matrix(axis, angle), p)
        rotated_xyz.append(rota_p)

    replaceXYZ(compdb, "rotated_DNA.pdb", rotated_xyz)

    os.system("cat transformed_protein.pdb rotated_DNA.pdb > {}".format(output))

main()
