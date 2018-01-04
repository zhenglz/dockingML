# -*- coding: utf-8 -*-

import dockml
import numpy as np
import pyemma as pe

class PCA :
    def __init__(self):
        pass

    def pca(self, data, no_ev):
        '''
        perform PCA analysis
        :param data:
        :param no_ev:
        :return:
        '''

        pca_obj = pe.coordinates.pca(data)

        return pca_obj

class essentialDynamics :
    def __init__(self):
        pass

    def transformXYZ(self, xyz, vector, delta):
        '''
        transform xyz along a vector
        eg. increase movements of an atom along the PC1 vectors
        :param xyz: list, xyz value
        :param vector: list, vectors for xyz
        :param delta: stride, unit nano meter
        :return: list, new list of xyz values
        '''

        newxyz = list(map(lambda x, y: x+ y*delta, xyz, vector))

        return newxyz

    def pdbIncreaseMotion(self, pdbin, vectors, delta=0.5):
        '''
        increase motions of a pdb given its PC component eigenvectors
        :param pdbin:
        :param vectors:
        :param delta:
        :return:
        '''

        pdbio = dockml.coordinatesPDB()

        with open(pdbin) as lines :
            lines = [ x for x in lines if ("ATOM" in x or "HETATM" in x)]

            coords = pdbio.getAtomCrdFromLines(lines)

            newxyzs =[]
            newlines=[]

            if len(vectors) == len(coords) :
                for i in range(len(coords)) :
                    newxyz = self.transformXYZ(coords[i], vectors[i], delta)
                    newxyzs.append(newxyz)
                    newlines.append(pdbio.replaceCrdInPdbLine(lines[i], newxyz))

        return newxyzs, newlines

    def genEDAEssemble(self, pdbin, pdbout, vectors, no_files=20, delta=0.5):
        '''
        generate an essemble of pdb files to increase the PC motions
        :param pdbin:
        :param pdbout:
        :param vectors:
        :param no_files:
        :param delta:
        :return:
        '''

        with open(pdbout, 'w') as tofile :
            for i in range(no_files) :
                tofile.write("MODEL   %d \n"%i)
                t, nlines = self.pdbIncreaseMotion(pdbin, vectors, delta=i*delta)
                for x in nlines :
                    tofile.write(x)
                tofile.write("ENDMDL  \n")

        return 1

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
