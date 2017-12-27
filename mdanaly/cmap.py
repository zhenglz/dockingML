#!/usr/bin/env python
# -*- coding: utf-8 -*-

#####################################################
# Script for generating contact probability map     #
# Author: ZHENG Liangzhen                           #
# Email: LZHENG002@e.ntu.edu.sg                     #
# Version: V3.1                                     #
# Date: 23 Nov 2017                                 #
#####################################################

import glob, math, sys, os
import numpy as np
import argparse

import dockml

from matplotlib import pyplot as plt
from collections import defaultdict
from argparse import RawTextHelpFormatter
from datetime import datetime
from mpi4py import MPI

import time

class DrawCMap :
    def __init__(self):
        pass

    def drawTimeSeries2D(self, cmapmatrix, refpdb=[],
                         fsize=14,
                         xlabel="", ylabel="",
                         cmaptype="Grey",
                         xlim=[], ylim=[],
                         yticks_loc=[], yticks_labels=[],
                         yticks_showchainid = False,
                         xticks_loc=[], xticks_labels=[],
                         colorbar_label="", colorbar_show=False,
                         savefig="",
                         ):

        """
        plot the ligand protein interactions (time series, x axis)
        :param cmapmatrix: str, the data file, containing time series matrix file
        :param refpdb: list, [ pdbfilename, chain id, residue sequence shift-by ]
        :param fsize:
        :param xlabel:
        :param ylabel:
        :param cmaptype:
        :param xlim:
        :param ylim:
        :param yticks_loc:
        :param yticks_labels:
        :param yticks_showchainid: whether show chainid of protein residues
        :param xticks_loc:
        :param xticks_lables:
        :return:
        """

        # load cmap file
        cmapdata = np.loadtxt(cmapmatrix, delimiter=",")
        cm_sorted = sorted(list(cmapdata), key=lambda x: x[0], reverse=False)

        # get key protein residues involving protein ligand binding
        key_res = []
        true_res = np.sum(np.asarray(cm_sorted)[:, 1:], axis=0) > 0
        for res in range(true_res.shape[0]):
            if true_res[res]:
                key_res.append(res)
        print("KEY RES ", key_res)

        # get full residue name index list
        res_labels = []
        if len(refpdb) == 3 :

            ppdb = dockml.pdbIO.parsePDB("")
            fullreslist = ppdb.getNdxForRes(refpdb[0], [refpdb[1]])

            shortresmap = ppdb.longRes2ShortRes()
            fullreslist = [ x for x in fullreslist if x[2] in refpdb[1] ]

            for resk in key_res :
                if fullreslist[resk][0] in shortresmap.keys() :
                    resseq = str(resk + refpdb[2])
                    resname= shortresmap[fullreslist[resk][0]]
                    chainid= fullreslist[resk][2]

                    if yticks_showchainid :
                        id = resname + resseq + chainid
                    else :
                        id = resname + resseq

                    #print("ID "* 5, id)
                    res_labels.append(id)

        # only keep the important residue cmap
        keyres_cmap = np.asarray(cm_sorted)[:, 1:][:, list(key_res)]

        # get the length of x and y axis
        shapex = len(key_res)
        shapey = cmapdata.shape[0] + 1

        print("Protein residue numbers: %d" % shapex)
        print("Time point numbers: %d"% shapey)

        #x = np.reshape(np.tile(range(shapex), shapey), (shapey, shapex))
        #y = np.asarray(range(shapey))
        #y = np.reshape(np.tile(y, shapex), (shapex, shapey)).T
        z = np.transpose(keyres_cmap[:, 1:]).T

        plt.pcolormesh(z.T, cmap=plt.get_cmap(cmaptype))

        if colorbar_show :
            plt.colorbar(label=colorbar_label)

        plt.xlabel(xlabel, fontsize=fsize)
        plt.ylabel(ylabel, fontsize=fsize)

        if len(yticks_loc) and len(yticks_labels) :
            plt.yticks(yticks_loc, yticks_labels)
        else :
            if len(refpdb) == 3 :
                plt.yticks(np.array(range(shapex))+0.5, res_labels)

        if len(xlim) :
            plt.xlim(xlim)

        if len(ylim) :
            plt.ylim(ylim)
        else :
            plt.ylim([0, shapex+0.5])

        if len(xticks_labels) and len(xticks_loc) :
            plt.xticks(xticks_loc, xticks_labels)

        if len(savefig) :
            plt.savefig(savefig, dpi=2000)

        plt.show()

        return 1


class ContactMap:
    def __init__(self, hugePDBFile) :
        self.mfPDB = hugePDBFile

    def atomDistance(self, atomCrd1, atomCrd2, sqrt=False):
        """
        calculate atom distance
        :param atomCrd1:
        :param atomCrd2:
        :param sqrt:
        :return:
        """

        if sqrt :
            return math.sqrt(sum(map(lambda x, y: (x-y)**2, atomCrd1, atomCrd2)))
        else :
            return sum(map(lambda x, y: (x-y)**2, atomCrd1, atomCrd2))

    def extractReference(self):
        '''
        input: None
        :return: the reference structure of the large pdb file, str
        '''
        reference = self.mfPDB+"_reference.pdb"
        if not os.path.exists(reference) :
            with open(reference, 'wb') as tofile :
                with open(self.mfPDB) as lines :
                    for s in lines :
                        if s[:6].strip() in ["ATOM", "HETATM"] :
                            tofile.write(s)
                        elif s.split()[0] in ["END", "ENDMDL"] :
                            tofile.write(s)
                            break
        return reference

    def splitPdbFile(self, fileHeader="MODEL"):
        '''
        read a large pdb file to obtain the splitted files
        :return: a list of single-frame pdb file
        '''
        """
        This function split the large Multi-frame PDB file
        Only the atoms in the atomNdx will be written to seperated files
        """

        if 'S_1.pdb' in os.listdir('./') :
            filelist = glob.glob('S_*.pdb')

        else :
            filelist = []
            count = 0
            with open(self.mfPDB) as lines :
                for s in lines :
                    if fileHeader == s.split()[0] :
                        count += 1
                        no_model = str(count)
                        tofile = open("S_"+no_model+".pdb",'wb')
                        filelist.append("S_"+no_model+".pdb")
                        tofile.write(s)
                    elif "ATOM" == s.split()[0] :
                        tofile.write(s)

                    elif "ENDMDL" in s :
                        tofile.write(s)

                        tofile.close()
                    else :
                        pass

        print( "PDB File spliting finished! \n")
        return (sorted(filelist))

    def getAtomNdxByAtomName(self, singleFramePDB, atomName):
        '''
        input a pdb file and return related atom index
        :param singleFramePDB: input pdb file, str
        :param resName: residue names, list of str
        :return: atom index, list of str
        '''
        ndx = []
        try :
            with open(singleFramePDB) as lines :
                ndx = [ s.split()[1] for s in lines
                        if (s.split()[0] in ["ATOM", "HETATM"] and s.split()[2] in atomName) ]
        except FileExistsError :
            print("File %s not exists" % (singleFramePDB))

        return ndx

    def getResIndex(self, singleFramePDB, chain, resIdList):
        """
        get the index of the residues from the pdbfile
        indexlist, [ resndx+"_"+chain, ]
        :param singleFramePDB:
        :param chain:
        :param resIdList:
        :return: list,
        """
        indexlist = []
        with open(singleFramePDB) as lines :
            indexlist = [ s[22:28].strip() + '_' + chain
                          for s in lines
                          if "ATOM" in s and s[21] == chain and
                          int(s[22:26].strip()) in resIdList and
                          s[22:28].strip() + '_' + chain not in indexlist
                          ]

            '''for s in lines:
                if "ATOM" in s and s[21] == chain and int(s[22:26].strip()) in resIdList:
                    resIndex = s[22:28].strip() + '_' + chain
                    if resIndex not in indexlist:
                        indexlist.append(resIndex)
            '''
        return indexlist

    def findAtomTypePerEle(self, elements, singleFramePDB):
        '''
        input a list of elements and a reference pdb file
        return the related atom names
        :param elements:
        :param singleFramePDB:
        :return:
        '''
        ndx = dockml.index.PdbIndex()
        atominfor = ndx.atomInformation(singleFramePDB)
        atomType = []
        keys = atominfor.keys()
        for key in keys :
            if atominfor[key][-1] in elements and atominfor[key][0] not in atomType :
                atomType.append(atominfor[key][0])
        return atomType

    def findAtomType(self, information, singleFramePDB):
        '''
        find atomtype of the atom
        :param information: a list of subgroup names for protein reisudes
        :param singleFramePDB: a reference pdb file
        :return:
        '''
        atomType = []

        if information in ['Alpha-Carbon', 'CA', 'alpha', 'Ca', 'Alpha']:
            atomType = ['CA']
        elif information in ['main', 'mainchain', 'Mainchain', 'MainChain']:
            atomType = ['CA', 'N', 'C', 'O']
        elif information in ['back', 'backbone', 'Backbone', 'BackBone']:
            atomType = ['CA', 'N']

        elif information in ['noH', 'non-H', 'non-hydrogen', 'Non-H', 'no-h',
                             'heavy', 'lig-all_heavy', 'lig-all_noH']:
            with open(singleFramePDB) as lines:
                for s in lines:
                    if "ATOM" in s and s.split()[2] not in atomType and s[13] != "H" and s[-1] != "H":
                        atomType.append(s.split()[2])

        elif information in ['sidechain', 'side'] :
            with open(singleFramePDB) as lines:
                for s in lines:
                    if "ATOM" in s and s.split()[2] not in atomType and s[13] != "H" :
                        if s.split()[2] not in ['CA', 'N', 'C', 'O'] :
                            atomType.append(s.split()[2])

        elif information in ['all', 'all-atom', 'All-Atom', 'ALL', 'lig-all', ]:
            with open(singleFramePDB) as lines:
                for s in lines:
                    if "ATOM" in s and s.split()[2] not in atomType:
                        atomType.append(s.split()[2])
                        
        elif len(information):
            atomType = [information]
        else:
            print( "Error! AtomType not correct. Exit Now! \n")
            sys.exit(1)

        return atomType

    def findAtomNdx(self, pdbfile, resList, chain, atomType, verbose=False):
        '''
        give a pdb file, return the atom ndx needed
        :param pdbfile:
        :param resList: a list of index of residues
        :param chain: chains
        :param atomType: atom names
        :param verbose:
        :return: a list of atom ndx, string
        '''
        if verbose :
            print( pdbfile, resList, chain, atomType)

        atomndx = []
        for key in resList.keys() :
            with open(pdbfile) as lines :
                atomndx += [ s.split()[1]
                            for s in lines
                            if len(s.split()) and
                            s[:6].strip() in ["ATOM", "HETATM"] and
                            s.split()[2] in atomType and
                            s[21] in list(chain) and
                            int(s[22:26].strip()) in resList[key] and
                            s.split()[1] not in atomndx
                            ]
                '''
                for s in lines :
                    if "ATOM" in s or "HETATM" in s :
                        if s.split()[2] in atomType and s[21] in list(chain) and
                            int(s[22:26].strip()) in resList[key] :
                            if s.split()[1] not in atomndx :
                                append(s.split()[1])
                            else :
                                pass '''
        return atomndx

    def getPdbCrd(self, singleFramePDB, atomList, perAtom=False):
        '''
        input a pdb file return the atom crds in each residue
        in a dictionary format
        :param singleFramePDB:
        :param atomList:
        :return:
        '''
        resCrds = []
        resRecord = []

        with open(singleFramePDB) as lines:
            crdPerRes = []
            for s in lines:
                if s[:5].strip() in ['ATOM', 'HETATM'] and "TER" not in s:

                    if s.split()[1] in atomList :
                        if perAtom :
                            if s[21] + s.split()[1] not in resRecord :
                                resRecord.append(s[21] + s.split()[1])
                                if len(crdPerRes):
                                    resCrds.append(crdPerRes)
                                crdPerRes = []
                        else :
                            if (s[21] + s[22:26].strip()) not in resRecord :
                                resRecord.append(s[21] + s[22:26].strip())
                                if len(crdPerRes) :
                                    resCrds.append(crdPerRes)
                                crdPerRes = []
                        # resIndex = s.split()[4 + len(s[21].strip())] + '_' + s[21]
                        crd_list = []
                        crd_list.append(float(s[30:38].strip()))
                        crd_list.append(float(s[38:46].strip()))
                        crd_list.append(float(s[46:54].strip()))

                        crdPerRes.append(crd_list)

            resCrds.append(crdPerRes)

        ## resCrdDist format : {'123':[[0.0,0.0,0.0],[]]}
        return resCrds

    def getResidueName(self, singleFramePDB, resList, chains, perAtom=False, verbose=False):
        '''
        input a pdb file, return the residue names, NDX_Chain
        in a list format
        :param singleFramePDB: str, reference pdb file name
        :param resList: list, residue sequence list
        :param chains: list, a list of chains
        :return:
        '''
        used = []
        for key in resList.keys():
            with open(singleFramePDB) as lines:
                lines = list(filter(lambda x: ("ATOM" in x or "HETATM" in x) and x[21] in chains, lines))
                if perAtom:
                    resNames = [x.split()[1] + '_' + x[21]
                                for x in lines
                                if int(x[22:26].strip()) in resList[key]]
                else:
                    resNames = [x[22:26].strip() + '_' + x[21]
                                for x in lines
                                if int(x[22:26].strip()) in resList[key]]

                for item in resNames:
                    if item not in used:
                        used.append(item)
        if verbose:
            print(used, len(used))
        return used

    def residueContacts(self, resCrd1, resCrd2,
                        distcutoff, countcutoff=1.0,
                        switch=False, verbose=False,
                        rank=0, NbyN=False,
                        ):
        '''
        read the coordinates of two residue,
        return whether contacts
        :param resCrd1: 2d list
        :param resCrd2: 2d list
        :param distcutoff: the squared the distance cutoff
        :param countcutoff: number of contacts within a residue
        :param switch: apply a switch function
        :param verbose: verbose
        :return: contact if 1, no contact if zero
        '''
        newlist = []

        for item1 in resCrd1:
            newlist += [[item1, item2] for item2 in resCrd2]

        distances = [sum(map(lambda x, y: (x - y) ** 2, pair[0], pair[1])) for pair in newlist]
        if verbose :
            print( rank, " DISTANCES ", distances)

        if NbyN :
            '''
            for community analysis
            fij = N / sqrt(N_i * N_j)
            '''
            counts = len(filter(lambda x: x <= distcutoff, distances))
            return float(counts) / np.sqrt(float(len(distances)))
        else :
            if switch:
                dc = 2 * math.sqrt(distcutoff)
                swf= dockml.algorithms.BasicAlgorithm()
                counts = swf.switchFuction(math.sqrt(distances[0]), dc)
                if verbose :
                    print( rank, " COUNTS: ", counts)
                return counts
            else :
                counts = len(filter(lambda x: x <= distcutoff, distances))

                if verbose :
                    print( "Counts here ", counts)
                if counts >= countcutoff :
                    return 1.0
                else:
                    return 0.0

    def subgroupCmap(self, pdbin, cutoff, atomNdx=[], verbose=False, logifle='log.log'):
        """
        calculate sub-group based contactmap
        :param pdbin:
        :param cutoff:
        :param atomNdx:
        :param verbose:
        :param logifle:
        :return:
        """
        if not os.path.exists(pdbin):
            # raise Exception("Boo! \nNot find PDB files for calculation! ")
            print( 'Exit Now!')
            sys.exit(1)

        tofile = open(logifle, 'w')

        # elements in atom infor
        # key: atom index
        # value: [atomname, molecule type, is_hydrogen,  resname, resndx, chainid,(mainchian, sidechain,
        #           sugar ring, phosphate group, base ring)]
        # atominfor[atomndx] = [atomname, moltype, is_hydrogen, resname, resndx, chain, subgroup]
        ndx = dockml.index.PdbIndex()
        detailatomInfor = ndx.atomInformation(pdbin)

        atomndx_1 = atomNdx[0]
        atomndx_2 = atomNdx[-1]

        pdbio = dockml.pdbIO.coordinatesPDB()
        crds1 = pdbio.getAtomCrdByNdx(pdbin, atomndx_1)
        crds2 = pdbio.getAtomCrdByNdx(pdbin, atomndx_2)

        distance_cutoff = cutoff ** 2

        for y in range(len(atomndx_2)) :
            atom2 = atomndx_2[y]
            infor2 = detailatomInfor[atom2]
            if infor2[-1] in ['P','S', 'N','O'] :
                for x in range(len(atomndx_1)) :
                    atom1 = atomndx_1[x]
                    infor1 = detailatomInfor[atom1]
                    if infor1[-1] in ['P','S', 'N','O'] :
                        crd1 = crds1[x]
                        crd2 = crds2[y]
                        if verbose :
                            print(infor1, infor2)
                            print(crd1, crd2)

                        distance =  sum(map(lambda x, y: (x - y) ** 2, crd1, crd2))
                        if verbose :
                            print( distance)

                        if distance <= distance_cutoff :

                            resid1 = infor1[3] + "_" + str(infor1[4]) + "_" + infor1[6]
                            resid2 = infor2[3] + "_" + str(infor2[4]) + "_" + infor2[6]

                            print( resid1, resid2)

                            tofile.write("%16s   %-20s   %-4s  %-4s  %8.3f\n"%
                                         (resid2, resid1, infor2[0], infor1[0], math.sqrt(distance)))

            if verbose :
                print("complete atom %s "% atom2)

        tofile.close()

        return 1

    def cmap_timeseries(self, pdbFileList, cutoff=3.5, deltaT=0.001,
                        switch=False, atomNdx=[],
                        rank=0, perAtom=[False, False],
                        ccutoff = 0.0, NbyN=False,
                        nresidues=0, verbose=False,
                        ):

        """
        calculate time series based cmap
        :param pdbFileList: list, a list of pdb files
        :param cutoff: float, distance cutoff
        :param deltaT: float, time delta
        :param switch: bool, use switch function
        :param atomNdx: list, [ receptor atom list, ligand atom list]
        :param rank: int, the process cpu thread id
        :param perAtom: bool, whether calculate atom-atom matrix
        :param ccutoff: float, a cutoff for contacts between residues
        :param NbyN: bool, whether use N by N to normalize contact
        :param nresidues: number of total pairs
        :param verbose: bool, detail output
        :return: list, list of lists
        """

        if len(pdbFileList) == 0:
            # raise Exception("Boo! \nNot find PDB files for calculation! ")
            print('Exit Now!')
            sys.exit(1)

        atomndx_1 = atomNdx[0]
        atomndx_2 = atomNdx[-1]

        distance_cutoff = cutoff ** 2

        if ccutoff > 0:
            countCutoff = ccutoff
        else:
            if len(atomNdx[0]) * len(atomNdx[-1]) > nresidues :
                countCutoff = 2.0
            else:
                countCutoff = 1.0

        progress = 0
        tscmap = []

        for pdbfile in pdbFileList:

            progress += 1

            # time stamp of the pdb file
            time = int(pdbfile.split("_")[1].split(".")[0]) * deltaT
            print("Rank %d Progress: The %dth File %s out of total %d files" % \
                  (rank, progress, pdbfile, len(pdbFileList)))
            if verbose:
                print(rank, atomndx_1, atomndx_2)

            # get residue based list of xyz coordinates
            resCrdDict1 = self.getPdbCrd(pdbfile, atomndx_1, perAtom=perAtom[0])
            resCrdDict2 = self.getPdbCrd(pdbfile, atomndx_2, perAtom=perAtom[1])

            nmax = len(resCrdDict2)

            contCountMap = []

            # looping over all residues
            for m in range(len(resCrdDict1)):
                for n in range(len(resCrdDict2)):

                    if verbose:
                        print(rank, " RESIDUES ", m, n, resCrdDict1[m], resCrdDict2[n])

                    # start calculate the contact map for each frame
                    contacter = self.residueContacts(resCrd1=resCrdDict1[m],resCrd2=resCrdDict2[n],
                                                     distcutoff=distance_cutoff,countcutoff=countCutoff,
                                                     switch=switch,verbose=verbose,
                                                     rank=rank, NbyN=False,)
                    # append the cmap vector to the dataset
                    contCountMap.append(contacter)

            tscmap.append([time] + contCountMap)

        # sort the 2d list, time stamp is re-ordered in accending order
        #tscmap = sorted(tscmap, key=lambda x: x[0], reverse=False)

        return tscmap

    def cmap_ca(self, pdbFileList, cutoff, switch=False, atomNdx=[],
                rank=0, verbose=False, nresidues=0,
                perAtom=[False, False], ccutoff = 0.0, NbyN=False,
                ):
        '''
        construct the contact map
        :param pdbFileList:
        :param cutoff: distance cutoff in angstrom
        :param switch: appy a switch function
        :param atomNdx: atom index
        :param rank: number of cpu threads
        :param verbose: bool, detail output
        :param nresidues: int, number of residues
        :param perAtom: bool, caculate per atom cmap
        :param ccutoff: int, count cutoff
        :param NbyN: boo, an NbyN contact scheme
        :return: list, a list of contacts,
        '''

        if len(pdbFileList) == 0:
            # raise Exception("Boo! \nNot find PDB files for calculation! ")
            print( 'Exit Now!')
            sys.exit(1)

        atomndx_1 = atomNdx[0]
        atomndx_2 = atomNdx[-1]

        distance_cutoff = cutoff ** 2

        if ccutoff > 0 :
            countCutoff = ccutoff
        else :
            if len(atomNdx[0]) * len(atomNdx[-1]) > nresidues :
                countCutoff = 2.0
            else:
                countCutoff = 1.0

        ### start calculate all the pdbfile residue c alpha contact
        #contCountMap = [0] * nresidues
        progress = 0
        contCountMap = [0.0] * nresidues
        oldtime = datetime.now()

        for pdbfile in pdbFileList:
            progress += 1
            print( "Rank %d Progress: The %dth File %s out of total %d files" % \
                  (rank, progress, pdbfile, len(pdbFileList)))

            delta = datetime.now() - oldtime
            print("Rank %d Progress: Time Usage for 1 frame %d seconds" % (rank, delta.total_seconds))
            oldtime = datetime.now()

            if verbose :
                print( rank, atomndx_1, atomndx_2)

            resCrdDict1 = self.getPdbCrd(pdbfile, atomndx_1, perAtom=perAtom[0])
            resCrdDict2 = self.getPdbCrd(pdbfile, atomndx_2, perAtom=perAtom[1])

            nmax = len(resCrdDict2)

            for m in range(len(resCrdDict1)):
                for n in range(len(resCrdDict2)):
                    if verbose :
                        print( rank, " RESIDUES ", m, n, resCrdDict1[m], resCrdDict2[n])
                    contCountMap[n + m * nmax] += self.residueContacts(resCrdDict1[m],
                                                                       resCrdDict2[n],
                                                                       distance_cutoff,
                                                                       countCutoff,
                                                                       switch,
                                                                       verbose,
                                                                       rank,
                                                                       NbyN=NbyN
                                                                       )

            print( "PDB file " + str(pdbfile) + " Finished!")
            del resCrdDict1, resCrdDict2
        if verbose :
            print(contCountMap)

        return contCountMap

    def writeFiles(self, cmap, nFiles, cmapName, resNames1, resNames2, rank):
        """
        output the contact map to a file
        :param cmap: list, vector of contact map data for all pairs
        :param nFiles: int, dividing factor to averaging the cmap
        :param cmapName: str, output file name
        :param resNames1: list, list of residue names
        :param resNames2: list, list of residue names
        :param rank: int, process rank
        :return:
        """
        if len(cmap) != len(resNames1) * len(resNames2) :
            print("Fatal Error! Write data to file failed!")
            sys.exit(0)
        else :
            nmax = len(resNames2)
            # generate contact matrix file
            result = open(str(rank) + "_" + cmapName, 'wb')
            result.write("NDX  ")
            for i in resNames2 :
                result.write('%5s ' % str(i))

            # print contCountMap.keys()
            for m in range(len(resNames1)):
                result.write('\n%5s ' % resNames1[m])
                for n in range(len(resNames2)):
                    result.write('%8.1f ' % (cmap[n + m * nmax ]))
            result.close()

            result = open(str(rank) + "_" + cmapName + '.xyz', 'wb')
            result.write("# Receptor Ligand Contact_probability \n")
            for i in range(len(resNames1)):
                for j in range(len(resNames2)) :
                    result.write('%5s  %5s  %8.4f \n' % (resNames1[i], resNames2[j], cmap[j+i*nmax]/nFiles ))
            result.close()

        return 1

    def arguements(self):
        d = '''
        ########################################################################
        #  Script for generating contact probability map                       #
        #  Author:  ZHENG Liangzhen                                            #
        #  Email:   LZHENG002@e.ntu.edu.sg                                     #
        #  Version: V2.2                                                       #
        #  Date:    27 Dec 2017                                               #
        ########################################################################

        Generating contact probability Map (Cmap)

        Input a multi-frame pdb (MFPDB file) file to construct a contact probability map.
        This MFPDB have multiple models in a single file, and all the structures should
        stay whole, broken structures will cause inaccurate results.
        All the frames in MFPDB do not consider PBC conditions, you should keep structures
        as a whole.

        If some arguements not given, default values would be used.

        Usage:
        1. Show help information
        python ContactMap.py -h

        2. Construct a Ca-Ca Cmap for a protein chain
        python ContactMap.py -inp MF_pro.pdb -out Cmap.dat -rc A 1 250
        -lc A 1 250 -cutoff 3.5 -switch T -atomtype CA

        3. Generate a full-atom Cmap for a poly-peptide chain
        python ContactMap.py -inp MF_pro.pdb -out Cmap.dat -rc A 1 250
        -lc A 1 250 -cutoff 3.5 -atomtype all all

        4. Construct a Cmap between a small ligand and a protein
        python ContactMap.py -inp MF_pro.pdb -out Cmap.dat -rc A 1 250
        -lc A 251 251 -cutoff 3.5 -atomtype all all

        5. Construct a Cmap between a small ligand and a protein, Ca-allatom
        python ContactMap.py -inp MF_pro.pdb -out Cmap.dat -rc A 1 250
        -lc A 251 251 -cutoff 3.5 -atomtype CA all

        6. Construct a cmap between a protein chain with MPI
        mpirun -np 4 python ContactMap.py -inp MF_pro.pdb -out Cmap.dat -rc A 1 250
        -lc A 251 251 -cutoff 3.5 -atomtype CA all -np 4

        '''
        # parse arguments
        parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
        parser.add_argument('-inp', type=str, help="The input huge PDB file with multiple frames")
        parser.add_argument('-out', type=str, default='ContactMap.dat',
                            help="The output file name. Default name is ContactMap.dat \n")
        parser.add_argument('-rc', type=str, nargs='+', default=['A', '1', '250'],
                            help="The receptor chains and residue index for Cmap construction.\n"
                                 "You must enter a chain name, start residue index, and end chain index.\n"
                                 "Default is: A 1 250 \n")
        parser.add_argument('-lc', type=str, nargs='+', default=['A', '1', '250'],
                            help="The ligand chains and residue index for Cmap construction.\n"
                                 "You must enter a chain name, start residue index, and end chain index.\n"
                                 "Default is: B 1 250 \n")
        parser.add_argument('-cutoff', type=float, default=0.35,
                            help="Distance Cutoff for determining contacts. \n"
                                 "Default is 3.5 (angstrom). \n")
        parser.add_argument('-atomtype', type=str, nargs='+', default=[],
                            help="Atom types for Receptor and Ligand in Contact Map Calculation. \n"
                                 "Only selected atoms will be considered.\n"
                                 "Options: CA, Backbone, MainChain, All, non-H(All-H), lig-all. \n"
                                 "CA, alpha-carbon atoms. Backbone, backbone atoms in peptides. \n"
                                 "MainChain, including CA and N atoms. All, means all atoms.\n"
                                 "non-H, non-hydrogen atoms, all the heavy atoms. \n"
                                 "lig-all trys to consider all the atoms of a ligand (H atoms not considered). \n"
                                 "Two choices should be provided for receptor and ligand respectively. \n"
                                 "If only one atomtype given, the 2nd will be the same as 1st.\n"
                                 "Default is: [] \n")
        parser.add_argument('-atomname1', type=str, nargs='+', default=[],
                            help="Atom names for Recetpor in Contact Map. \n"
                                 "Default is []. ")
        parser.add_argument('-atomname2', type=str, nargs='+', default=[],
                            help="Atom names for Ligand in Contact Map. \n"
                                 "Default is []. ")
        parser.add_argument('-eletype', type=str, nargs="+", default=[],
                            help="Choose the specific elements for atom indexing to construct the cmap."
                                 "Default is empty.")
        parser.add_argument('-switch', type=str, default='True',
                            help="Apply a switch function for determing Ca-Ca contacts for a smooth transition. \n"
                                 "Only work with atomtype as CA. Options: True(T, t. TRUE), False(F, f, FALSE) \n"
                                 "Default is False. \n")
        parser.add_argument('-np', default=0, type=int,
                            help='Number of Processers for MPI. Interger value expected. \n'
                                 'If 4 is given, means using 4 cores or processers.\n'
                                 'If 1 is given, means not using MPI, using only 1 Core.\n'
                                 'Default is 1. ')
        parser.add_argument('-test', default=0, type=int,
                            help="Do a test with only a number of frames. For example, 4 frames. \n"
                                 "Default value is 0. ")
        parser.add_argument('-NbyN', type=bool, default=False,
                            help="For community analysis, calculate atom contact number, normalized. \n"
                                 "Default is False.")
        parser.add_argument('-verbose', default=False, type=bool,
                            help="Verbose. Default is False.")
        parser.add_argument('-details', default=None, type=str,
                            help="Provide detail contact information and write out to a file. \n"
                                 "Default is None."
                            )

        parser.add_argument('-opt', default="Separated", type=str,
                            help="Optional setting controls. \n"
                                 "Separated, using separated files instead of splitting files.\n")

        args = parser.parse_args()

        # decide to print help message
        if len(sys.argv) < 2:
            # no enough arguements, exit now
            parser.print_help()
            print("You chose non of the arguement!\nDo nothing and exit now!\n")
            sys.exit(1)

        return args

def main() :
    ## change to working directory
    pwd = os.getcwd()
    os.chdir(pwd)

    startTime = datetime.now()

    d = '''
    ########################################################################
    #  Generating contact probability map                                  #
    #  Author:  ZHENG Liangzhen & Mu Yuguang                               #
    #  Email:   LZHENG002@e.ntu.edu.sg                                     #
    #  Version: V2.2                                                       #
    #  Date:    27 Dec 2017                                                #
    ########################################################################

    Generating contact probability Map (Cmap)

    Input a multi-frame pdb (MFPDB file) file to construct a contact probability map.
    This MFPDB have multiple models in a single file, and all the structures should
    stay whole, broken structures will cause inaccurate results.
    All the frames in MFPDB do not consider PBC conditions, you should keep structures
    as a whole.

    If some arguements not given, default values would be used.

    Usage:
    1. Show help information
    python ContactMap.py -h

    2. Construct a Ca-Ca Cmap for a protein chain
    python ContactMap.py -inp MF_pro.pdb -out Cmap.dat -rc A 1 250
    -lc A 1 250 -cutoff 3.5 -switch T -atomtype CA

    3. Generate a full-atom Cmap for a poly-peptide chain
    python ContactMap.py -inp MF_pro.pdb -out Cmap.dat -rc A 1 250
    -lc A 1 250 -cutoff 3.5 -atomtype all all

    4. Construct a Cmap between a small ligand and a protein
    python ContactMap.py -inp MF_pro.pdb -out Cmap.dat -rc A 1 250
    -lc A 251 251 -cutoff 3.5 -atomtype all all

    5. Construct a Cmap between a small ligand and a protein, Ca-allatom
    python ContactMap.py -inp MF_pro.pdb -out Cmap.dat -rc A 1 250
    -lc A 251 251 -cutoff 3.5 -atomtype CA all

    6. Construct a cmap between a protein chain with MPI
    mpirun -np 4 python ContactMap.py -inp MF_pro.pdb -out Cmap.dat -rc A 1 250
    -lc A 251 251 -cutoff 3.5 -atomtype CA all -np 4

    '''
    # parse arguments
    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-inp', type=str, help="The input huge PDB file with multiple frames")
    parser.add_argument('-out',type=str, default='ContactMap.dat',
                        help="The output file name. Default name is ContactMap.dat \n")
    parser.add_argument('-rc',type=str,nargs='+', default=['A','1','250'],
                        help="The receptor chains and residue index for Cmap construction.\n"
                             "You must enter a chain name, start residue index, and end chain index.\n"
                             "Default is: A 1 250 \n")
    parser.add_argument('-lc', type=str, nargs='+', default=['A','1','250'],
                        help="The ligand chains and residue index for Cmap construction.\n"
                             "You must enter a chain name, start residue index, and end chain index.\n"
                             "Default is: B 1 250 \n" )
    parser.add_argument('-cutoff',type=float,default=0.35,
                        help="Distance Cutoff for determining contacts. \n"
                             "Default is 3.5 (angstrom). \n")
    parser.add_argument('-atomtype',type=str,nargs='+', default=[],
                        help="Atom types for Receptor and Ligand in Contact Map Calculation. \n"
                             "Only selected atoms will be considered.\n"
                             "Options: CA, Backbone, MainChain, All, non-H(All-H), lig-all. \n"
                             "CA, alpha-carbon atoms. Backbone, backbone atoms in peptides. \n"
                             "MainChain, including CA and N atoms. All, means all atoms.\n"
                             "non-H, non-hydrogen atoms, all the heavy atoms. \n"
                             "lig-all trys to consider all the atoms of a ligand (H atoms not considered). \n"
                             "Two choices should be provided for receptor and ligand respectively. \n"
                             "If only one atomtype given, the 2nd will be the same as 1st.\n"
                             "Default is: [] \n")
    parser.add_argument('-atomname1', type=str, nargs='+', default=[],
                        help="Atom names for Recetpor in Contact Map. \n"
                             "Default is []. ")
    parser.add_argument('-atomname2', type=str, nargs='+', default=[],
                        help="Atom names for Ligand in Contact Map. \n"
                             "Default is []. ")
    parser.add_argument('-eletype', type=str, nargs="+", default=[],
                        help="Choose the specific elements for atom indexing to construct the cmap."
                             "Default is empty.")
    parser.add_argument('-switch', type=str, default='True',
                        help="Apply a switch function for determing Ca-Ca contacts for a smooth transition. \n"
                             "Only work with atomtype as CA. Options: True(T, t. TRUE), False(F, f, FALSE) \n"
                             "Default is False. \n")
    parser.add_argument('-np', default=0, type=int,
                        help='Number of Processers for MPI. Interger value expected. \n'
                             'If 4 is given, means using 4 cores or processers.\n'
                             'If 1 is given, means not using MPI, using only 1 Core.\n'
                             'Default is 1. ')
    parser.add_argument('-test', default=0, type=int,
                        help="Do a test with only a number of frames. For example, 4 frames. \n"
                             "Default value is 0. ")
    parser.add_argument('-NbyN', type=bool, default=False,
                        help="For community analysis, calculate atom contact number, normalized. \n"
                             "Default is False.")
    parser.add_argument('-verbose', default=False , type=bool,
                        help="Verbose. Default is False.")
    parser.add_argument('-details', default=None, type=str,
                        help="Provide detail contact information and write out to a file. \n"
                             "Default is None."
                        )
    parser.add_argument('-opt', default="TimeSeries", type=str,
                        help="Optional setting controls. \n"
                             "Separated, using separated files instead of splitting files."
                             "Options: S(Separated), TS (TimeSeries).\n"
                             "Default is TimeSeries. \n")

    args = parser.parse_args()

    if args.np > 0 :

        # load the large mutiple frame pdb file
        cmap = ContactMap(args.inp)

        reference = cmap.extractReference()

        # decide to print help message
        if len(sys.argv) < 2:
            # no enough arguements, exit now
            parser.print_help()
            print( "You chose non of the arguement!\nDo nothing and exit now!\n")
            sys.exit(1)

        # setup MPI environment
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.rank

        # define the atoms for contact map constructions
        atomType = []

        atomName1 = args.atomname1
        atomName2 = args.atomname2

        if len(args.atomtype) not in [1, 2]:
            if len(args.eletype) != 0 :
                atomType = [ cmap.findAtomTypePerEle(args.eletype, reference),
                             cmap.findAtomTypePerEle(args.eletype, reference) ]
            elif len(atomName1) and len(atomName2) :
                atomType = [ args.atomname1, args.atomname2 ]
            else :
                print( "Define atom indexing failed")
                sys.exit(0)
        else :
            atomType.append(cmap.findAtomType(args.atomtype[0], reference))
            atomType.append(cmap.findAtomType(args.atomtype[0 + len(args.atomtype) - 1], reference))

        if atomType == [['CA'],['CA']] :
            switch = args.switch
        else :
            switch = False

        # receptor information about residues
        rcResNdx = defaultdict(list)
        round = len(args.rc) / 3
        rcChains = []
        for i in range(round):
            rcResNdx[args.rc[i * 3]] = range(int(args.rc[i * 3 +1]),int(args.rc[i *3+2]) + 1)
            rcChains.append(args.rc[(i + 1) * 3 - 3])

        # ligand information about residues
        lcResNdx = defaultdict(list)
        round = len(args.lc) / 3
        lcChains = []
        for i in range(round):
            lcChains.append(args.lc[i*3])
            lcResNdx[args.lc[i * 3]] = range(int(args.lc[i * 3 + 1]), int(args.lc[i * 3 + 2]) + 1)

        ## start to construct map
        receptorAtomNdx = cmap.findAtomNdx(reference, rcResNdx, rcChains, atomType[0], args.verbose)
        ligandAtomNdx   = cmap.findAtomNdx(reference, lcResNdx, lcChains, atomType[-1], args.verbose)

        if args.verbose :
            print( "ATOM RES AND CHAIN " * 5)
            print( lcResNdx, rcResNdx, lcChains, rcChains)
            print( "ATOM NDX")
            print( receptorAtomNdx, ligandAtomNdx)

        # report detail interactions, verbose reports
        if args.details :
            cmap.subgroupCmap(reference, args.cutoff,
                              [receptorAtomNdx, ligandAtomNdx],
                              args.verbose, args.details
                              )
        else :
            if rank == 0:
                if args.test:
                    pdbFileList = sorted(cmap.splitPdbFile())[: args.test]
                else:
                    pdbFileList = sorted(cmap.splitPdbFile())
            else:
                pdbFileList = None

            # board cast the list of files into different threads
            pdbFileList = comm.bcast(pdbFileList, root=0)
            totalNumOfFiles = len(pdbFileList)

            if rank == 0 :
                load4each = int(math.ceil(float(totalNumOfFiles) / float(args.np)))
                filesList = []

                for i in range(args.np - 1) :
                    filesList.append(pdbFileList[i * load4each : load4each * (i+1)])
                filesList.append(pdbFileList[(args.np-1)*load4each:])

                if args.verbose:
                    print( "Full File List " * 10, pdbFileList, filesList)

                if args.opt in ["S", "Separated", "separated"] :
                    filesList = [[args.inp,]]
            else :
                filesList = None

            ## scatter data to sub-processers/threads
            filesList = comm.scatter(filesList, root=0)

            if "lig-all" in args.lc[0] :
                allatoms = [False, True]
            else :
                allatoms = [False, False]

            recNames = cmap.getResidueName(filesList[0], rcResNdx, rcChains, perAtom=allatoms[0])
            ligNames = cmap.getResidueName(filesList[0], lcResNdx, lcChains, perAtom=allatoms[1])

            Cmap = cmap.cmap_ca(filesList, args.cutoff,switch,
                                [receptorAtomNdx, ligandAtomNdx],
                                rank, args.verbose, len(recNames)*len(ligNames),
                                perAtom=allatoms, ccutoff=1.0, NbyN=args.NbyN
                                )

            cmap.writeFiles(Cmap, len(filesList), args.out, recNames, ligNames, rank)

            overallValuesList = comm.gather(Cmap, root=0)

            ## once calculation done, wrap up and write data out
            if rank == 0 :
                # "Wrap Up and write data to files "
                final = np.zeros(len(Cmap))
                for rank_values in overallValuesList :
                    final += np.asarray(rank_values)

                cmap.writeFiles(final, totalNumOfFiles, args.out, recNames, ligNames, rank='all')

        if rank == 0:
            print( "Total Time Usage: ")
            print( datetime.now() - startTime)

        MPI.Finalize()
        sys.exit(1)
