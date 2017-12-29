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

        newxyz = map(lambda x, y: x+ y*delta, xyz, vector)

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
            if len(vectors) == coords :
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

        with open(pdbout, 'wb') as tofile :
            for i in range(no_files) :
                tofile.write("MODEL   %d \n"%i)
                nlines = self.pdbIncreaseMotion(pdbin, vectors, delta=i*delta)
                for x in nlines :
                    tofile.write(x)

        return 1
