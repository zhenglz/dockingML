# -*- coding: utf-8 -*-
#!/usr/bin/env python

import dockml
import numpy as np
import sys
import argparse


class essentialDynamics(object):

    def __init__(self):
        pass

    def transformXYZ(self, xyz, vector, delta):
        """
        transform xyz along a vector
        eg. increase movements of an atom along the PC1 vectors

        Parameters
        ----------
        xyz: list,
            xyz value
        vector: list,
            vectors for xyz
        delta: float,
            stride, unit nano meter

        Returns
        -------
        newxyz: list,
            new list of xyz values

        """

        newxyz = list(map(lambda x, y: x + y*delta, xyz, vector))

        return newxyz

    def pdbIncreaseMotion(self, pdbin, vectors, delta=0.5):
        """
        increase motions of a pdb given its PC component eigenvectors

        Parameters
        ----------
        pdbin
        vectors: ndarray, shape=[ M, 3]
            a M * 3 matrix, M means num of Ca atoms or heavy atoms
        delta: float, default=0.5
            the stride for atom movement

        Returns
        -------
        newxyzs : np.ndarray, shape = [ M, 3]
            the new xyz coordinates, M is number of atoms
        newlines : list, length = M
            the new pdb lines, M is number of atoms

        """

        pdbio = dockml.coordinatesPDB()

        with open(pdbin) as lines:
            lines = [x for x in lines if ("ATOM" in x or "HETATM" in x)]

            coords = pdbio.getAtomCrdFromLines(lines)

            newxyzs =[]
            newlines=[]

            if vectors.shape[0] == len(coords):
                for i in range(len(coords)):
                    newxyz = self.transformXYZ(coords[i], list(vectors[i]), delta)
                    newxyzs.append(newxyz)
                    newlines.append(pdbio.replaceCrdInPdbLine(lines[i], newxyz))

        return newxyzs, newlines

    def genEDAEssemble(self, pdbin, pdbout, vectors, no_files=20, delta=0.5, numres=250):
        """Generate an essemble of pdb files to increase the PC motions

        Parameters
        ----------
        pdbin
        pdbout
        vectors
        no_files
        delta
        numres

        Returns
        -------

        """

        PI = 3.14159

        with open(pdbout, 'w') as tofile :
            for i in range(no_files) :
                length = delta * np.cos(2.0 * PI * (float(i) / float(no_files)) - PI )
                print(length)
                tofile.write("MODEL   %d \n"%i)
                t, nlines = self.pdbIncreaseMotion(pdbin, vector, delta=length)
                for x in nlines :
                    tofile.write(x)
                tofile.write("ENDMDL  \n")

        return self

    def averageVectors(self, vectors, resindex):
        '''
        extract selected (resindex) rows and return their averages
        :param vectors: list of lists,
        :param resindex: list, residue list
        :return:
        '''

        domain_vectors = []
        for res in resindex :
            domain_vectors.append(vectors[res])

        aver_vector = np.sum(np.asarray(domain_vectors), axis=0)
        #print(aver_vector)
        return aver_vector

    def domainWiseEigVec(self, domainf, vectors, scalefactor=1.0, output='aver-vectors.dat'):
        '''
        averaging the vectors on each atoms in a domain
        :param domainf: str, domain information data file
        :param vectors: list of lists, dimension N*3
        :param output: str, output file name
        :return: tuple, (domain_average_vectors, domain names)
        '''

        dom = dockml.parsePDB()
        domains = dom.readDomainRes(domainf)

        # domain_names
        d_name = [ x[0] for x in domains ]

        minResIndx = min(sum([ x[1:] for x in domains ], []))

        aver_vec = []
        for i in range(len(domains)) :
            #print(domains[i])
            # becasue vectors index starting from 0
            resindexlist = []
            for k in range((len(domains[i])-1)/2) :
                #print(domains[i][k*2+1], domains[i][k*2+2])
                resindexlist += list(np.asarray(range(domains[i][k*2+1], domains[i][k*2+2])) - minResIndx)

            v = self.averageVectors(vectors, resindexlist) * scalefactor
            aver_vec.append(v)

        np.savetxt(output, np.asarray(aver_vec),
                   delimiter=' ', fmt='%12.5f',
                   header=" ".join(d_name), comments="#"
                   )

        return (aver_vec, d_name)

    def domainCOM(self, domainf, ref, pdbchain='A', output="com_domains.day", atomNames=["CA"]):
        """
        calculate centor of mass, or geometry center of a domain
        :param domainf:str, domain information data file
        :param ref: str, reference pdb file
        :param pdbchain: chain id
        :param output: str, output file name
        :return: list of list, com of domains, dimension N*3
        """

        dom = dockml.parsePDB()
        domains = dom.readDomainRes(domainf)

        dnames = [x[0] for x in domains ]

        pdbc = dockml.coordinatesPDB()
        ndx = dockml.PdbIndex()

        coms = []

        for i in range(len(domains)) :

            # get atom index of residues in a domain
            atomindex = ndx.res_index(inpdb=ref,
                                      chain=pdbchain,
                                      atomtype="",
                                      atomList=atomNames,
                                      residueNdx=domains[i][1:],
                                      )

            # get crds of a list of atoms
            crds = pdbc.getAtomCrdByNdx(ref, atomindex)

            # calculate geometry center, not centor of mass
            com = np.mean(np.asarray(crds), axis=0)

            coms.append(list(com))

        np.savetxt(output, np.asarray(coms), fmt="%8.3f", delimiter=" ", header=dnames, comments="#")

        return coms


if __name__ == "__main__":

    d = '''
    Generate EDA ensemble of conformations with PCA eigenvectors
    
    Usage:
    python dynamics.py input.pdb output.pdb eigenvectors.dat number_of_files delta_stepsize
    
    Example:
    
    '''

    parser = argparse.ArgumentParser(description=d)
    parser.add_argument("-f", type=str, help="Input. The reference pdb file name.")
    parser.add_argument("-o", type=str, help="Output. The output multiple-frame pdb file name.")
    parser.add_argument("-vect", type=str, help="Input. The eigenvector file name. Assuming that "
                                                "this file is generated and output directly from "
                                                "gmx_pca.py. The shape of the vector matrix is N*M, "
                                                "where M is total number of dimensions and N is number"
                                                "of PC projections.")
    parser.add_argument("-nf", type=int, default=50,
                        help="Output. The number of frames in the output file. ")
    parser.add_argument("-delta", type=float, help="Input. Default is 0.5. \n"
                                                   "The movement stride for the coordinates. ")
    parser.add_argument("-pc", type=int, default=1,
                        help="Input, optional, default is 1. Which PC eigenvectors for essential dynamics move. ")

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()

    vector = np.reshape(np.loadtxt(args.vector)[args.pc -1], (-1, 3))

    dyn = essentialDynamics()
    dyn.genEDAEssemble(args.f, args.o, vector, args.nf, args.delta, vector.shape[0])