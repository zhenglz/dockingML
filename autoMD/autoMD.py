#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import urllib
import time
import argparse
import subprocess as sp
from argparse import RawTextHelpFormatter
import math
from glob import glob
from collections import *
import linecache
import numpy as np
import glob
import time

from dockingML import extract

# import modeller for loop refinement
try:
    from modeller import *
    from modeller.automodel import *
    MODELLER_EXIST = True
except ImportError :
    print("Warning: Modeller is not well imported, some features may not work. ")
    MODELLER_EXIST = False

class RewritePDB :
    def __init__(self, filename):
        self.pdb = filename

    def pdbRewrite(self, output, chain, atomStartNdx, resStartNdx):
        # change atom id, residue id and chain id
        resseq = resStartNdx
        atomseq = int(atomStartNdx)
        chainname = chain
        lines = open(self.pdb)
        newfile = open(output,'w')
        resseq_list = []

        for s in lines :
            if "ATOM" in s and len(s.split()) > 6 :
                atomseq += 1
                newline = s
                newline = self.atomSeqChanger(newline, atomseq)
                newline = self.chainIDChanger(newline, chainname)
                if len(resseq_list) == 0 :
                    newline = self.resSeqChanger(newline, resseq)
                    resseq_list.append(int(s[22:26].strip()))
                else :
                    if resseq_list[-1] == int(s[22:26].strip()) :
                        newline = self.resSeqChanger(newline, resseq)
                    else :
                        resseq += 1
                        newline = self.resSeqChanger(newline, resseq)
                    resseq_list.append(int(s[22:26].strip()))
                newfile.write(newline)
            else :
                newfile.write(s)
        print "Completed!"

    # Change the ID of residue, or index of residue
    def resSeqChanger(self, inline, resseq):
        resseqstring = " "*(4- len(str(resseq)))+str(resseq)
        newline = inline[:22] + resseqstring + inline[26:]
        return newline

    def atomSeqChanger(self, inline, atomseq) :
        atomseqstring = " " * (5 - len(str(atomseq))) + str(atomseq)
        newline = inline[:6] + atomseqstring + inline[11:]
        return newline

    def resNameChanger(self, inline, resname) :
        resnamestr = " " * ( 4 - len(str(resname) ) ) + str(resname)
        newline = inline[:17] + resnamestr + inline[20:]
        return newline

    def chainIDChanger(self, inline, chainid) :
        newline = inline[:21] + str(chainid) + inline[22:]
        return newline




class PrepScripts :
    def __init__(self, qsubScriptSample):
        self.qsubMD = qsubScriptSample

    def prepareQsub(self,  outSH, cmds):
        tofile = open(outSH, 'wb')

        with open(self.qsubMD) as lines:
            for s in lines:
                if "#define_commends" in s:
                    for cmd in cmds :
                        tofile.write(cmd+" \n")
                else:
                    tofile.write(s)
        tofile.close()

    def RewriteAntechamber(self, nc, antechamber="sample_antechamber.sh"):
        tofile = open("antechamber.sh", 'wb')
        with open(antechamber) as lines:
            for s in lines:
                if len(s.split()) > 0 and s.split()[0] == "antechamber":
                    tofile.write(
                        "antechamber -i $1 -fi mol2 -o amberff.prep.$2 -fo prepi -at $2 -pf y -s 2 -c bcc -nc %d \n" % nc)
                else:
                    tofile.write(s)
        return 1

class AutoRunMD :
    def __init__(self, topFile, taskName, grompp, mdrun, verbose=True, qsub=False):
        self.top = topFile
        self.task= taskName
        self.verbose = verbose
        self.grompp = grompp
        self.mdrun = mdrun

    def modifyMDP(self, inputMDPFile, outputMDPFile, parameters):
        tofile = open(outputMDPFile, 'wb')
        with open(inputMDPFile) as lines :
            for s in lines :
                if len(s.split()) > 0 and s[0] != ";" :
                    if s.split()[0] in parameters.keys() \
                            and len(parameters[s.split()[0]]) > 0 :
                        tofile.write("%s    = ")
                        for item in parameters[s.split()[0]] :
                            tofile.write("%s "%str(item))
                        tofile.write(" \n")
                    else:
                        tofile.write(s)
                else :
                    tofile.write(s)
        return(outputMDPFile)

    def energyMinimization(self, emMDPFile, groFile, outTprFile, np=4, qsub=False):
        if self.verbose :
            print("Start energy minimization for task %s."%self.task)

        # generate GMX tpr file
        job = sp.check_output("%s -f %s -c %s -o %s -p %s "%
                              (self.grompp, emMDPFile, groFile, outTprFile, self.top),
                              shell=True )
        if self.verbose :
            print(job)
            print("Generating TPR file for EM completed.")

        # run MD with qsub or on the terminal
        cmd = "mpirun -np %d %s -deffnm %s " % (np, self.mdrun, outTprFile)
        if qsub :
            pscript = PrepScripts('sample_qsub.sh')
            pscript.prepareQsub('em_qsub.sh', cmd)
            job = sp.Popen("qsub em_qsub.sh", shell=True)
            job.communicate()

        else :
            job = sp.Popen(cmd, shell=True)
            job.communicate()

        return(outTprFile+".gro")

    def runMD(self, MDPFile, groFile, outTprFile, mdpOptions,
              sampleScript='sample_qsub.sh',
              qsub=False, np=4,
              index="index.ndx",
              ):

        # prepare mdp file for MD
        current_mdp = MDPFile
        if len(mdpOptions.keys()) :
            self.modifyMDP(inputMDPFile=MDPFile,
                           outputMDPFile="modified_"+MDPFile,
                           parameters=mdpOptions)
            current_mdp = "modified_"+MDPFile

        # prepare tpr file for md
        job = sp.check_output("%s -f %s -c %s -o %s -p %s -n %s"
                              %(self.grompp, current_mdp, groFile, outTprFile, self.top, index),
                              shell=True)
        if self.verbose :
            print(job)

        # run md
        cmd = "mpirun -np %d %s -deffnm %s " % (np, self.mdrun, outTprFile[-4:])
        if qsub :
            pscript = PrepScripts(sampleScript)
            oscript = str(len(glob("./*_"+sampleScript)))+"_"+sampleScript
            pscript.prepareQsub(oscript, cmd)
            job = sp.Popen("qsub %s"%oscript, shell=True)
            job.communicate()
        else :
            job = sp.Popen(cmd+" &", shell=True)
            job.communicate()

        return(outTprFile[-4:]+".gro")

class ContactMap:
    def __init__(self, hugePDBFile) :
        self.mfPDB = hugePDBFile

    def splitPdbFile(self):

        if 'S_1.pdb' in os.listdir('./') :
            filelist = glob.glob('S_*.pdb')

        else :
            filelist = []
            with open(self.mfPDB) as lines :
                for s in lines :
                    if "MODEL" == s.split()[0] :
                        no_model = s.split()[1]
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

        print("PDB File spliting finished! \n")
        return filelist

    def switchFuction(self, x, d0, m=12, n=6):
        # for countting, implement a rational switch function to enable a smooth transition
        # the function is lik  s= [1 - (x/d0)^6] / [1 - (x/d0)^12]
        # d0 is a cutoff
        y = (1.0 - math.pow((x / d0), n)) / (1.0 - math.pow((x / d0), m))
        return y

    def getResIndex(self, singleFramePDB, chain, resIdList) :
        '''
        read in a pdb file, return the residueNdx_chain information
        :param singleFramePDB:
        :param chain:
        :param resIdList:
        :return:
        '''
        lines = open(singleFramePDB)
        indexlist = []
        for s in lines :
            if "ATOM" in s and s[21] == chain and int(s[22:26].strip()) in resIdList :
                resIndex = s[22:28].strip() + '_'+chain
                if resIndex not in indexlist :
                    indexlist.append(resIndex)
        lines.close()
        return indexlist

    def findAtomType(self, information, singleFramePDB):
        atomType = []

        if information in ['Alpha-Carbon','CA','alpha','Ca','Alpha'] :
            atomType = ['CA']
        elif information in ['main','mainchain', 'Mainchain','MainChain'] :
            atomType = ['CA', 'N', 'C', 'O']
        elif information in ['back','backbone', 'Backbone','BackBone'] :
            atomType = ['CA', 'N']

        elif information in ['noH','non-H', 'non-hydrogen', 'Non-H', 'no-h', 'heavy', 'heavy-atom'] :
            with open(singleFramePDB) as lines :
                for s in lines :
                    if len(s.split()) and s.split()[0] in ['ATOM','HETATM'] and s.split()[2] not in atomType and s[13] != "H" and s[-1] != "H" :
                        atomType.append(s[12:16].strip())

        elif information in ['all','all-atom','All-Atom','ALL'] :
            with open(singleFramePDB) as lines :
                for s in lines :
                    if len(s.split()) and s.split()[0] in ['ATOM','HETATM'] and s.split()[2] not in atomType :
                        atomType.append(s[12:16].strip())
        else :
            print("Error! AtomType not correct. Exit Now! \n")
            sys.exit(1)

        return atomType

    def getPdbCrd(self, singleFramePDB, atomType, resList) :
        # input a pdb file return the atom crds in each residue
        # in a dictionary format
        resCrdDict = defaultdict(list)

        with open(singleFramePDB) as lines :
            for s in lines:
                if "ATOM" in s and "TER" not in s :
                    resIndex = s[22:28].strip() + '_' + s[21]
    #                #if s.split()[2] in atomType and resIndex in resList :
                    # non-hydrogen bonds
                    if s.split()[2] in atomType and resIndex in resList :
                        #resIndex = s.split()[4 + len(s[21].strip())] + '_' + s[21]
                        if resIndex not in resCrdDict.keys() :
                            resCrdDict[resIndex] = []
                        crd_list = []
                        crd_list.append(float(s[30:38].strip()))
                        crd_list.append(float(s[38:46].strip()))
                        crd_list.append(float(s[46:54].strip()))

                        resCrdDict[resIndex].append(crd_list)

        lines.close()
        ## resCrdDist format : {'123':[[0.0,0.0,0.0],[]]}
        return resCrdDict

    def coordinationNumber(self, pdbfiles, cutoff,
                           rResNdxs, lResNdxs,
                           chains,
                           atomTypes=['heavy', 'all'],
                           switch=False,
                           maxCoordinate=100,
                           ):
        '''
        check coordinate number (contact number) between receptor and ligand
        :param pdbfiles:
        :param cutoff:
        :param rResNdxs:
        :param lResNdxs:
        :param chains:
        :param atomTypes:
        :param switch:
        :param maxCoordinate:
        :return:
        '''
        coordinate = 0.0

        rec_atomtype = self.findAtomType( atomTypes[0],pdbfiles[0])
        lig_atomtype = self.findAtomType( atomTypes[-1],pdbfiles[-1])

        rec_Crds = self.getPdbCrd(pdbfiles[0], rec_atomtype, rResNdxs)
        lig_Crds = self.getPdbCrd(pdbfiles[-1], lig_atomtype, lResNdxs)

        for rec_Res in rec_Crds.keys() :
            for lig_Res in lig_Crds.keys() :
                for crd1 in rec_Crds[rec_Res] :
                    for crd2 in lig_Crds[lig_Res] :
                        distance = 0.0
                        crd1 = np.asarray(crd1)
                        crd2 = np.asarray(crd2)
                        distance = np.sum((crd1 - crd2) ** 2)
                        count = 1.0

                        if switch :
                            count = self.switchFuction(distance, cutoff*2)

                        if distance <= cutoff :
                            coordinate += count

                        if coordinate >= maxCoordinate :
                            return(coordinate)

        return(coordinate)

    def contactMap(self, pdbFileList, cutoff, cmapName,
                   rchainRes, lchainRes, switch,
                   atomType = [["CA",]],
                   rank=0
                   ):
        distCutOff = cutoff * cutoff   ###  square of cutoff
        MaxDistCutoff = distCutOff * 16

        if len(atomType[0]) > 1 or len(atomType[1]) > 1 :
            condition = True
        else :
            condition = False

        #print "Condition ", condition

        if len(pdbFileList) == 0 :
            #raise Exception("Boo! \nNot find PDB files for calculation! ")
            print('Exit Now!')
            sys.exit(1)

        ### identify all the residue index
        resIndexList1 = []
        for chain in rchainRes.keys() :
            resIndexList1 += self.getResIndex(pdbFileList[0], chain, rchainRes[chain])
        #print resIndexList1
        resIndexList2 = []
        for chain in lchainRes.keys() :
            resIndexList2 += self.getResIndex(pdbFileList[0], chain, lchainRes[chain])
        #print resIndexList2

        ### start calculate all the pdbfile residue c alpha contact
        contCountMap = OrderedDict()
        progress = 0
        for pdbfile in pdbFileList :
            progress += 1
            print("Rank %d Progress: The %dth File %s out of total %d files" %
                  (rank, progress, pdbfile, len(pdbFileList)))
            resCrdDict1 = self.getPdbCrd(pdbfile, atomType[0], resIndexList1)
            #print resCrdDict1

            resCrdDict2 = self.getPdbCrd(pdbfile, atomType[1], resIndexList2)

            #print resCrdDict1

            for m in resIndexList1 :
                for n in resIndexList2 :
                    id = m +"+" + n
                    if id not in contCountMap.keys() :
                        contCountMap[id] = 0.0
                    ### compute the distance between any combination of
                    ### of residue indexes and compare with cutoff square
                    count = 0.0
                    for crd1 in resCrdDict1[m] :
#                        distance = 0.0
                        for crd2 in resCrdDict2[n] :
                            distance = 0.0
                            crd1 = np.asarray(crd1)
                            crd2 = np.asarray(crd2)
                            distance = np.sum((crd1 - crd2)**2)
                            #print id, "distance ", distance
                            #print "distance distance ", distance

                            if distance > MaxDistCutoff * 25 :
                                #print "Break ", n, m
                                break

                            else :
                                #print "ATOMTYPE", atomType[0], atomType
                                #if  1 : #atomType[0] == "CA" : #and atomType[1] in ['CA', 'Ca', 'ca', 'alpha-C', 'Calpha'] :
                                #print "CA " * 10

                                if not condition :
                                    if switch in ['T', 'True', 't', 'TRUE', 'true'] and distance <= MaxDistCutoff :
                                        contCountMap[id] += \
                                            self.switchFuction(math.sqrt(distance), cutoff*2.0)
                                        #print "Counting", contCountMap[id],  distance

                                    elif distance <= distCutOff :
                                        contCountMap[id] += 1.0

                                else :
                                    if distance <= distCutOff:
                                        count += 1.0

                                if condition and count >= 2.0 :
                                    contCountMap[id] += 1.0
                                    #print "id ", id, " counts 1"
                                    break

                        if distance > MaxDistCutoff * 25 :
                            #print "MAX Distance", distance
                            break

                        if count >= 2.0 :
                            break

            print "PDB file "+str(pdbfile)+" Finished!"
            #print contCountMap[id]

        # generate contact matrix file
        result = open(str(rank)+"_"+cmapName,'wb')
        result.write("NDX  ")
        for i in resIndexList1 :
            result.write('%5s '%str(i))

        #print contCountMap.keys()
        for j in resIndexList2 :
            result.write('\n%5s '%str(j))
            for k in resIndexList1 :
                result.write('%8.1f '% (contCountMap[k+'+'+j] / len(pdbFileList)))

        result = open(str(rank)+"_"+cmapName+'.xyz', 'wb')
        result.write("# Receptor Ligand Contact_probability \n")
        for i in resIndexList1:
            for j in resIndexList2:
                result.write('%5s  %5s  %8.1f \n' % (i, j, contCountMap[i+'+'+j]/ len(pdbFileList)))

        return contCountMap.keys(), contCountMap.values(), resIndexList1, resIndexList2



def runGenTop() :
    '''
    instant a gmxTopBuilder class, and input some parameters,
        to generate amber and gromacs topology and coordinates

    if a small ligand provide, if not prep and frcmod exist,
        then am1/bcc charge based on amber antechamber sqm method
        will be calculated. in this process, tleap will be used to
        get bonded and non-bonded parameters
    :parameter: None
    :return: None
    '''
    top = GenerateTop()
    args = top.arguments()

    try :
        os.environ["AMBERHOME"] = "/home/liangzhen/software/amber14/"
    except :
        os.system("export AMBERHOME=/home/liangzhen/software/amber14/")

    # define the input coordinates information, a pdb, or mol2, or a sequence
    # of amino acids
    if args.aaseq:
        structure = args.aaseq
    elif args.inp:
        structure = args.inp
    else:
        structure = ''
        print('Error: No input structure defined.\n Exit Now!')
        sys.exit(0)

    if not args.resp:
        # not calculate atomic charges
        top.gmxTopBuilder(args.frcmod, args.prep,
                          structure, args.out,
                          args.ion, args.nion,
                          args.bt, args.size,
                          args.ff)
    else:
        # prepare resp charges for small ligand with antechamber
        # give a netcharge: args.nc
        sh = top.runAntechamber(args.nc)

        try:
            # set env for AMBERHOME
            os.environ["AMBERHOME"] = "/home/liangzhen/software/amber14/"
            job = sp.Popen("sh %s %s %s" % (sh, args.inp, args.out), shell=True)
            job.communicate()
        except :
            print("Antechamber generateing RESP and atomic parameters error! \n "
                  "Double check your input structure.")
            sys.exit(1)

        top.gmxTopBuilder("frcmod." + args.out,
                          "prep." + args.out,
                          structure, args.out,
                          args.ion, args.nion,
                          args.bt, args.size,
                          args.ff)

    if args.ion[0] == "X+" or args.bt == "NA":
        # parse a gromacs .top to a .itp file
        # print "ITP file created!"
        top.removeMolInfor(args.out)

if __name__ == "__main__" :
    # cd to PWD
    os.chdir(os.getcwd())

    if len(sys.argv) <= 2 :
        help_information = '''
        The tool is used for automating Gromacs based MD simulation.
        Several independent tools are provided.

        Separated tools including index, gentop, autoMD etc.

        Usage:

        python autoMD.py -h
        python autoMD.py gentop -h
        python autoMD.py index -h
        python autoMD.py extract -h
        '''
        print(help_information)
        sys.exit(1)
    else :
        if str(sys.argv[1]) == "index" :
            #print("Show help message: ")
            index = PdbIndex()
            index.genGMXIndex()
        #if 1 :
        elif str(sys.argv[1]) == "gentop" :
            runGenTop()

        elif str(sys.argv[1]) == "extract" :
            runExtract()
