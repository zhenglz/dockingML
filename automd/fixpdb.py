#!/usr/bin/env python
# -*- coding: utf-8 -*-

import mdanaly
import math
import os, sys
import urllib
from collections import OrderedDict
from os import environ
import subprocess as sp
from collections import defaultdict

from dockml import convert
#from rdkit import Chem

# import modeller for loop refinement
try:
    from modeller import *
    from modeller.automodel import *
    MODELLER_EXIST = True
except ImportError:
    print("Warning: Modeller is not well imported, some features may not work. ")
    MODELLER_EXIST = False


class SummaryPDB(object):

    def __init__(self, pdbfile, aminoLibFile="amino-acid.lib"):
        self.pdbfile = pdbfile

        if not os.path.exists(aminoLibFile) :
            PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
            DEFINITIONS_ROOT = os.path.join(PROJECT_ROOT,
                                            '../data/amino-acid.lib')
            aminoLibFile = DEFINITIONS_ROOT

        resShortName = {}
        with open(aminoLibFile) as lines:
            for s in [x for x in lines if "#" not in x ]:
                resShortName[s.split()[2]] = s.split()[3]
        self.resShortName = resShortName

    def centerOfMass(self, inputMol, atomNdx,
                     obabelexe='obabel', molBox=False):
        """Given a file (preferable PDB file format), if not,
        the file will be convert into a pdb file,
        and selected atom sequence number, determine the
        Center Of Mass of the coordinates

        Parameters
        ----------
        inputMol: str,
            input pdb file, or mol2 file, or pdbqt file
        atomNdx: list,
            the list of atom identifier
        obabelexe: str,
            the obabel executable command
        molBox: bool, default is False
            the pbc information of a box

        Returns
        -------
        com: list,
            the x, y, z coordinate of the com
        [xsize, ysize, zsize]
        """

        from dockml import pdbIO
        import numpy as np

        pdb = inputMol
        # convert the molecule file to pdb where necessary
        if inputMol.split(".")[-1] not in ['pdb', 'pdbqt']:
            converter = convert.Convert(obabel=obabelexe)
            converter.convert(inputMol, inputMol+".pdb")
            pdb = inputMol + ".pdb"

        coordinates = []

        with open(pdb) as lines:

            lines = [ x for x in lines if x.split()[0] in ['ATOM', 'HETATM']]
            coordinates = pdbIO.coordinatesPDB().getAtomCrdFromLines(lines)

        coordinates = np.asarray(coordinates)
        # get the mean values of x, y and z coordinates
        com = np.array(coordinates).mean(axis=0)
        xcenter, ycenter, zcenter = com[0], com[1], com[2]

        # when PBC box exists, add the PBC information
        if molBox:
            xsize = 2 * max(np.max(coordinates[:, 0]) - xcenter,
                        np.abs(np.min(coordinates[:,0])-xcenter))
            ysize = 2 * max(np.max(coordinates[:, 1]) - ycenter,
                        np.abs(np.min(coordinates[:, 1]) - ycenter))
            zsize = 2 * max(np.max(coordinates[:, 2]) - zcenter,
                        np.abs(np.min(coordinates[:, 2]) - zcenter))
        else:
            xsize, ysize, zsize = 100, 100, 100

        return com, [xsize, ysize, zsize]

    def netCharges(self, inputMol, ligName=None):
        """Deduce the total net charge of a molecule (mostly a
        small ligand).
        netCharges determine the net charge of a molecule,
        given a pdb file.
        if the input molecule is a pdbqt (Vina input), or pqr
        (APBS input type),
           a molecule name (residue code) is need.
        if a mol2 file provide, only the @<TRIPOS>ATOM field
        will be used, no ligand name required.

        other formats are not supported.

        last column of a pdbqt, pqr and mol2 file generally
        will be the atomic charge field,
          otherwise, a ValueError exception will be rasied.

        Parameters
        ----------
        inputMol: str,
            input file with atomic charges in the last column
        ligName: str,
            ligand code

        Returns
        -------
        netCharge: int,
            the net charge of a compound
        """

        extension = inputMol.split(".")[-1]
        netCharge = 0.0

        try:
            from rdkit import Chem

            m = Chem.MolFromPDBFile(inputMol)
            netCharge = Chem.GetFormalCharge(m)
            #netCharge = Chem.rdmolops.GetFormalCharge(m)

        except ModuleNotFoundError:
            if extension in ['pdbqt', 'pqr']:
                with open(inputMol) as lines:
                    for s in lines:

                        if s.split()[0] in ["ATOM", "HETATM"] and len(s.split()) > 5:
                            try:
                                if ligName and ligName in s:
                                    netCharge += float(s.split()[-1])
                                elif not ligName:
                                    netCharge += float(s.split()[-1])
                            except ValueError:
                                netCharge += 0.0
                                print("Last column in %s is not a float point charge value."%inputMol)
            elif extension in ['mol2']:

                with open(inputMol) as lines:
                    condition = 0
                    for s in lines:
                        if len(s.split()) and "@" in s:
                            if "<TRIPOS>ATOM" in s:
                                condition += 1
                            elif "@<TRIPOS>BOND" in s:
                                condition = 0
                            else:
                                pass
                        elif condition and len(s.split()):
                            try:
                                netCharge += float(s.split()[-1])
                            except ValueError:
                                netCharge += 0.0
                                print("Last column in %s is not a float point charge value." % inputMol)
                    # print("NET CHARGE " * 10)

        return int(netCharge)

    def details(self, verbose=False):
        '''
        obtain details of the pdbfile
        :param verbose:
        :return:
        '''
        chains = []
        resNdx = defaultdict(list)
        resName= defaultdict(list)
        resAtom= defaultdict(list)
        resNameNdx = {}
        with open(self.pdbfile) as lines :
            for s in [ x for x in lines if
                       ( len(x.split()) and  x.split()[0] in ["ATOM","HETATM"] and "TER" not in x)
                       ]  :
                # get chains
                if s[21] not in chains:
                    chains.append(s[21])

                ## residue index information
                if s[21] not in resNdx.keys() :
                    resNdx[s[21]] = []
                if s[22:26].strip() not in resNdx[s[21]] :
                    resNdx[s[21]].append((s[22:26].strip()))

                ## residue name information
                if s[21] not in resName.keys():
                   resName[s[21]] = []
                if s[17:20].strip() + "_" + s[22:26].strip() not in resName[s[21]] :
                   resName[s[21]].append(s[17:20].strip() + "_" + s[22:26].strip())

                # pdb file atoms name in each residue
                resId = (s[22:26].strip() + '_' + s[21]).strip()
                if resId not in resAtom.keys() :
                    resAtom[resId] = []
                resAtom[resId].append(s[12:16].strip())

                # residue index and name hash map
                resNameNdx[s[22:26].strip()+'_' + s[21].strip()] = s[17:20].strip()

        ## print some information
        if verbose :
            print("\nNumber of chains in this pdb file : ", len(chains), ". They are ", chains)

            print("\nNumber of residues in each chain: ")
            for chain in chains :
                print("Chain ", chain, " ",len(resNdx[chain]))

            print("\nResidue Names for each chain are : ")
            for chain in chains :
                print("For chain ", chain)
                if len(resNdx[chain]) > 10 :
                    for i in range(10) :
                        print(resNdx[chain][i], "  ", resName[chain][i])
                    print("......")
                    for j in range(10) :
                        print(resNdx[chain][-10+j], "  ", resName[chain][-10+j])
                else :
                    for i in range(len(resNdx[chain])) :
                        print(resNdx[chain][i], "  ", resName[chain][i])

        return chains, resNdx, resName, resAtom, resNameNdx

    def getFastaSeq(self, fastaFile):
        """obtain fasta sequence from the input file

        Parameters
        ----------
        fastaFile: str,
            the input fasta file

        Returns
        -------
        fastaseq: str,
            the result fasta sequence, or the plain text sequence

        """

        fastaseq = ''

        if os.path.exists(fastaFile):
            with open(fastaFile) as lines:
                for s in lines:
                    if '>' not in s:
                        #print s #strip()
                        fastaseq += s.strip("\n")

        return fastaseq

    def summary(self, chain, verbose=False):

        chains, resNdx, resName, resAtom, resNameNdx = self.details()

        proteinSequence = {}
        noProteinResNdx = defaultdict(list)
        noProteinResName= defaultdict(list)
        missingResNdx = defaultdict(list)
        fullResNdx = defaultdict(list)
        #for chain in chains :
        #resNdxNameDict = {}
        if len(resNdx[chain]) != len(resName[chain]) :
            print("Error in function SummaryPDB.details. \nNumber of index is different with number of residues.")
            print("Exit Now!")
            sys.exit(1)

        proResNdx = []
        for i in range(len(resNdx[chain])) :
            ## get index for all the protein residues in a specific chain
            if resName[chain][i].split("_")[0] in self.resShortName.keys() :
                proResNdx.append(int(resNdx[chain][i]))

                #resNdxNameDict[int(resNdx[chain][i])] = resName[chain][i].split("_")[0]
            else :
                ## non-protein residues information
                if chain not in noProteinResName.keys() :
                    noProteinResName[chain] = []
                if chain not in noProteinResNdx.keys() :
                    noProteinResNdx[chain]  = []

                noProteinResName[chain].append(resName[chain][i])
                noProteinResNdx[chain].append(resNdx[chain][i])

        # get protein chain sequence
        startNdx = proResNdx[0]
        finalNdx = proResNdx[-1]

        fullNdx = range(startNdx, finalNdx+1)
        fullResNdx[chain] = fullNdx

        #print fullNdx
        for i in range(len(fullNdx) ):
            residueid = str(fullNdx[i]) + "_" + chain

            if chain not in proteinSequence.keys() :
                proteinSequence[chain] = ''

            if residueid in resNameNdx.keys() :
                if resNameNdx[residueid] not in self.resShortName.keys() :
                    proteinSequence[chain] += "_" + resNameNdx[residueid]
                else :
                    proteinSequence[chain] += self.resShortName[resNameNdx[residueid]]
            else :
                proteinSequence[chain] += "-"

                ## find missing residues' index
                if chain not in missingResNdx.keys() :
                    missingResNdx[chain] = []
                missingResNdx[chain].append(str(fullNdx[i]))

        # print some information
        if verbose :
            print("\nFull sequence of protein chains are: ")
            #for chain in chains :
            print(chain, "  ", proteinSequence[chain])

            print("\nSome missing protein residues in each chain ")
            #for chain in chains :
            print(chain, "  ", missingResNdx[chain])

            print("\nThe information of the non-protein residues here ")
            #for chain in chains :
            print(chain, "  ", noProteinResName[chain])

            print("\nThe sequence of the full protein chain is: ")
            print("Chain  ", chain)
            sections = math.ceil(float(len(proteinSequence[chain])) / 20.0)

            for i in range(int(sections - 1)):
                print(fullResNdx[chain][i * 20], " " * (-len(str(fullResNdx[chain][i * 20])) + 11), fullResNdx[chain][i * 20 + 10])
                print(proteinSequence[chain][i * 20: i * 20 + 10], " ", proteinSequence[chain][i * 20 + 10: i * 20 + 20])
            print(fullResNdx[chain][i * 20 + 20])
            print(proteinSequence[chain][i * 20 + 20:])

        return proteinSequence, missingResNdx, noProteinResNdx, noProteinResName, fullResNdx

    def missingRes(self, chain, fastaSeq, verbose=False):
        ## find the missing residues sequences in the pdb file
        #chains, resNdx, resName, resAtom, resNameNdx = self.details()
        proteinSequence, missingResNdx, noProteinResNdx, noProteinResName, fullResNdx = self.summary(chain)
        #print fullResNdx[chain]
        trueResName = defaultdict(list)
        ## full fasta sequences here :
        fastaseq = ''
        if os.path.isfile(fastaSeq):
            fastaseq = self.getFastaSeq(fastaFile=fastaSeq)

        else :
            fastaseq = fastaSeq

        # print protein missing residue information
        startndx = fullResNdx[chain][0]
        for index in missingResNdx[chain] :
            if chain not in trueResName.keys() :
                trueResName[chain] = []
            trueResName[chain].append(fastaseq[int(index)-startndx]+'_'+index)

        ranges = []
        siteRange = []
        missedSeqsList = defaultdict(list)
        missedRSeq = ''
        brokenSeq = proteinSequence[chain]
        for i in range(len(brokenSeq)) :
            currentResNdx = i + int(startndx)
            if brokenSeq[i] != '-' and len(siteRange) == 0 :
                pass
            elif brokenSeq[i] == '-' :
                siteRange.append(currentResNdx)
                missedRSeq += fastaseq[i]
                #print missedRSeq
            elif brokenSeq[i] != '-' and len(siteRange) > 0 :
                ranges.append([siteRange[0],siteRange[-1]])
                #missedRSeq.append(seq[i])

                missedSeqsList[str(siteRange[0])+"_"+str(siteRange[-1])] = missedRSeq
                missedRSeq = ''
                siteRange = []
            else :
                pass

        if verbose :
            ## print full sequence information
            print("The fasta sequence of the full protein chain is: ")
            print("Chain  ", chain)
            sections = math.ceil(float(len(fastaseq)) / 20.0)

            for i in range(int(sections-1)) :
                print(fullResNdx[chain][i*20], " "*(-len(str(fullResNdx[chain][i*20]))+11), fullResNdx[chain][i*20+10])
                print(fastaseq[i*20 : i*20+10], " ",fastaseq[i*20+10 : i*20+20])
            print(fullResNdx[chain][i*20+20])
            print(fastaseq[i*20+20:])

            ## print the true residue name
            print("\nThe missing residues are :")
            print(chain,"  ", trueResName[chain])

            ## information about the missing sequences
            for ndxrange in missedSeqsList.keys() :
                print (("%12s   %-s ")%(ndxrange, missedSeqsList[ndxrange]))

        return fastaseq, missedSeqsList, fullResNdx

    def extendMissRange(self, chain, fastaSeq, extend=10, verbose=False):
        extMissSeq = {}
        fastaseq, missedSeqs, fullResNdx = self.missingRes(chain,fastaSeq)
        startndx = int(fullResNdx[chain][0])
        for ndxrange in missedSeqs.keys() :
            fixSeq = fastaseq[(int(ndxrange.split("_")[0])-extend-startndx):(int(ndxrange.split("_")[1]) + extend-startndx+1)]
            newrange = str((int(ndxrange.split("_")[0])-extend)) + "_" +\
                str(int(ndxrange.split("_")[1]) + extend)
            extMissSeq[newrange] = fixSeq

        if verbose :
            print( "\nExtend %d residues at broken part for chain %s : "%(extend,chain))
            for newrange in extMissSeq.keys() :
                print(("%12s   %-s")%(newrange, extMissSeq[newrange]))

        return extMissSeq


class FixPDB(object):
    def __init__(self):
        pass

    def pdbDownloader(self, pdbCode):
        if not os.path.exists(pdbCode + '.pdb'):
            try :
                source = urllib.urlopen("http://www.rcsb.org/pdb/files/" + pdbCode + ".pdb")
                with open(pdbCode + '.pdb', 'wb') as target:
                    target.write(source.read())
            except urllib.URLError :
                print( "URL error")

        return 1

    def saveHETATM(self, pdbin, chain=['A'], waterOnly=False):
        '''

        :param pdbin:
        :param chain:
        :param waterOnly: only save water molecule
        :return:
        '''
        with open(pdbin+"_HETATM", 'wb') as tofile :
            try :
                with open(pdbin) as lines :
                    for s in lines :
                        if len(s.split()) \
                                and s.split()[0] == "HETATM" \
                                and s[21] in chain :
                            if waterOnly and s[17:20].strip() in ["WAT","HOH"]:
                                tofile.write(s)
                            else :
                                tofile.write(s)
            except IOError :
                print("File %s not exist \n"% pdbin)

        return (pdbin + "_HETATM")

    def saveLigand(self, pdbin, chain, ligCode) :
        with open(pdbin + "_Ligand", 'wb') as tofile :
            try:
                with open(pdbin) as lines:
                    for s in lines:
                        if len(s.split()) \
                                and s.split()[0] in ["HETATM","ATOM"] \
                                and s[21] in chain\
                                and s[17:20].strip() == ligCode :
                            tofile.write(s)

            except IOError:
                print("File %s not exist \n" % pdbin)

        return (pdbin + "_Ligand")

    def checkCrash(self, pdbfiles, cutoff,
                   loopRange, ligandRange,
                   chainLoop, chainLigand,
                   maxCoordinate=10
                   ):
        '''
        check number of contacts between the newly modeled loop with the ligand
        if too many contacts, crashes are existed

        :param pdbfiles:
        :param cutoff:
        :param loopRange:
        :param ligandRange:
        :param chainLoop:
        :param chainLigand:
        :param maxCoordinate:
        :return: coordination number ( number of contacts)
        '''
        cmap = mdanaly.ContactMap(pdbfiles[0])
        coordinate = 0

        pdblist = []
        for item in pdbfiles :
            if '.pdb' not in item :
                if os.path.exists(item+'.pdb') :
                    pdblist.append(item+'.pdb')
                elif os.path.exists(item + '.pdbqt'):
                    pdblist.append(item + '.pdbqt')
                else :
                    print("File %s not found! "%item)
            else :
                pdblist.append(item)

        try :
            coordinate = cmap.coordinationNumber(pdblist, cutoff,
                                                range(loopRange[0], loopRange[-1]+1),
                                                range(ligandRange[0], ligandRange[-1]+1),
                                                [chainLoop, chainLigand],
                                                ['heavy','all'],
                                                False, maxCoordinate
                                                )
        except IOError:
            print("Calculate coordination number between receptor and ligand failed!")

        return (coordinate)

    def addLoopLOOPY(self, pdbin,
                          loopRange,
                          loopSeq,
                          loopchain,
                          loopyexe='loopy',
                          num_mode = 10,
                          iteration=200,
                          ligandName='',
                          ligandNdx=1,
                          ligandchain='',
                          verbose=False
                          ):
        '''
        perform loop building with LOOPY
        LOOPY: https://honiglab.c2b2.columbia.edu/software/Jackal/Jackalmanual.htm#loopy
        Usage: ./loopy -o=900-908 -r=ALSVEFPEM -n=10 -i=200 2ovh.pdb -c=A > 2ovh_899_909

        Notes: in the process of loop modeling, the Ligand and the HETATM atom records are
        all lost. thus if you intend to keep them, you could save them ahead using the
        self.saveLigand and self.saveHETATM methods

        :param pdbin:
        :param loopRange:
        :param loopSeq:
        :param loopyexe:
        :param verbose:
        :return:
        '''

        cmd = "%s -o=%d-%d -r=%s -n=%d -i=%d -c=%s %s > %s_adding_loop_%d_%d.log" % \
              (
                  loopyexe, loopRange[0], loopRange[1], loopSeq,
                  num_mode, iteration,
                  loopchain, pdbin, pdbin,
                  loopRange[0], loopRange[1]
              )
        if not verbose :
            cmd += " > /dev/null 2>&1 "

        output = sp.check_output(cmd, shell=True)
        finalModel = pdbin + "_looper_0.pdb"
        if verbose :
            print(output)

        model_crash = OrderedDict()
        for i in range(num_mode) :
            model =  pdbin + "_looper_" + str(i) + ".pdb"
            # check crash. Search for any crash for the loop and the ligand.
            crash = self.checkCrash([model, pdbin], 1.5, loopRange, [ligandNdx,], loopchain, ligandchain, 10)
            model_crash[model] = crash

        #model_crash = OrderedDict(sorted(model_crash.items(), key=lambda x: x[1], reverse=False))
        #crashIndx = 0
        key = pdbin + "_looper_" + str(i) + ".pdb"
        for key in model_crash.keys():
            if model_crash[key] <= 1.0 :
                finalModel = key
                break
            else :
                pass
        if key == len(model_crash.keys()[-1]) :
            model_crash = OrderedDict(sorted(model_crash.items(), key=lambda x: x[1], reverse=False))
            finalModel = model_crash.keys()[0]

        if verbose :
            print("Checking crash in the LOOPY models ")
            print(model_crash)

        if model_crash.values()[0] > 10.0 :
            print("Warning: the LOOPY model %s is not perfect, crashes occur! ")

        #finalModel = model_crash.keys()[0]

        return (finalModel)

    def removeRegions(self, filename, residuesNdx, chain, pdbout="temp_1.pdb", verbose=False):
        '''
        input a pdbfile, remove the selected residues

        :param filename: input pdbfile
        :param residuesNdx: the range of residue index for deleting
        :param chain: which chain to perform delete
        :param pdbout: the output pdb file name
        :return: 1
        '''
        tofile = open(pdbout,'w')
        with open(filename) as lines :
            for s in lines :
                if s.split()[0] in ["ATOM", "HETATM" ] and s[21] == chain and int(s[22:26].strip()) in residuesNdx :
                    pass
                else :
                    tofile.write(s)
        tofile.close()

        return(1)

    def addModeledRegions(self, basepbd, modeledpdb,
                          modelNdx, removeNdx, chain,
                          pdbout="temp.pdb",
                          verbose=False
                          ):
        tofile = open(pdbout, 'w')
        with open(basepbd) as lines:
            for s in lines :
                if s.split()[0] in ["ATOM", "HETATM" ] and s[21] == chain and int(s[22:26].strip()) < removeNdx[0] :
                    tofile.write(s)
                else :
                    pass

        ## addd modeled pdb here
        with open(modeledpdb) as lines :
            for s in lines :
                if s.split()[0] in ["ATOM", "HETATM"] and int(s[22:26].strip()) in modelNdx :
                    tofile.write(s)
                else :
                    pass

        ## add the following part of the original PDB
        with open(basepbd) as lines:
            for s in lines :
                if s.split()[0] in ["ATOM", "HETATM"] and s[21] == chain and int(s[22:26].strip()) > removeNdx[-1]:
                    tofile.write(s)
                else :
                    pass

        tofile.close()
        return(1)

    def addhydrogenReduce(self, pdbin, pdbout='outputH.pdb',reduce='reduce', flipH=True,verbose=True):
        if flipH:
            cmd = "%s -FLIP %s > %s " % (reduce, pdbin, pdbout)
        else:
            cmd = "%s -NOFLIP %s > %s " % (reduce, pdbin, pdbout)
        job = sp.check_output(cmd, shell=True)

        if verbose:
            print(job)

        return 1

    def addMissingAtoms(self, pdbin, pdb2pqrPy,
                        ff='amber', pqrout='output.pqr',
                        addhydrogen=True,
                        ph=7.0, verbose=False ):
        '''
        add missing heavy atoms, and assign hydrogens according to the ph

        :param pdbin:
        :param pdb2pqrPy:
        :param ff:
        :param pqrout:
        :param addhydrogen:
        :param ph:
        :return:
        '''
        if addhydrogen :
            cmd = "python %s -ff %s --with-ph %f --ffout %s --chain %s %s " %\
                  (pdb2pqrPy, ff, ph, ff+"_"+pqrout, pdbin, pqrout)
        else :
            cmd = "python %s -ff %s --ffout %s --chain %s %s " %\
                  (pdb2pqrPy, ff, ff+"_"+pqrout, pdbin, pqrout)
        job = sp.Popen(cmd, shell=True)
        job.communicate()
        job.kill()

        return(ff+"_"+pqrout)

    def addLoopsAutoModel(self, pdbCode, chain1,
                         alignCode, chain2,
                         loopRange,
                         verbose=False):
        '''
        model missing loop region
        :param pdbCode:
        :param chain1:
        :param alignCode:
        :param chain2:
        :param loopRange:
        :param verbose: print more information
        :return: the final model pdb file
        '''
        if MODELLER_EXIST :
            if verbose :
                print("MODELLER exist, perform alignment and loop refinement.")
                log.verbose()

            for pdb in [pdbCode]+[alignCode] :
                if not os.path.exists(pdb+'.pdb') :
                    self.pdbDownloader(pdb, pdb+'.pdb')

            env = environ()
            # directories for input atom files
            env.io.atom_files_directory = ['.', '../atom_files']

            aln = alignment(env)
            mdl_1 = model(env, file=pdbCode, model_segment=('FIRST:%s'%chain1, 'LAST:%s'%chain1))
            aln.append_model(mdl_1, align_codes=pdbCode, atom_files=pdbCode+'.pdb')
            mdl_2 = model(env, file=alignCode, model_segment=('FIRST:%s'%chain2, 'LAST:%s'%chain2))
            aln.append_model(mdl_2, align_codes=alignCode, atom_files=alignCode+'.pdb')
            aln.align2d()
            aln.write(file='alignment.ali', alignment_format='PIR')

            """ select missing loop residues for modeling only"""
            class MyModel(automodel):
                # overide select_atom function, select only some residues
                def select_atoms(self):
                    # Select residues loopRange[0] to loopRange[1] (PDB numbering)
                    return selection(self.residue_range(str(loopRange[0])+":"+str(chain1),
                                                        str(loopRange[1])+":"+str(chain1)
                                                        )
                                     )

            a = MyModel(env,
                        alnfile='alignment.ali',
                        knowns=alignCode,
                        sequence=pdbCode,
                        )
            # a.auto_align()
            # get an automatic loop building
            a.make()

            # obtain successfully modeled pdb file names
            return (self.selectBestModeller(a.outputs))

    def addLoopsSimpleModel(self, pdbIn, chain1,
                            fastaSeq,
                            loopRanges,
                            no_lig=0,
                            verbose=False
                            ):
        '''
        add loop using no template
        :param pdbIn:
        :param chain1:
        :param fastaSeq:
        :param loopRanges:
        :param no_lig:
        :param verbose:
        :return:
        '''
        if pdbIn[-4:] == ".pdb" :
            pdbIn = pdbIn[:-4]

        env = environ()
        env.io.atom_files_directory = ['.', '../atom_files']

        if verbose :
            log.verbose()

        m = model(env, file=pdbIn)
        aln = alignment(env)
        aln.append_model(m, align_codes=pdbIn)
        aln.write(file=pdbIn + '.seq')

        ## add full sequence into the alignment file
        with open("alignment.ali",'wb') as tofile :
            with open(pdbIn + '.seq') as lines :
                tofile.write(lines)
            tofile.write('\nP1;%s_fill \n'%pdbIn)
            tofile.write('sequence::::::::: \n')
            for i in range(len(fastaSeq)) :
                tofile.write(fastaSeq[i])
                if i % 74 == 0 and i != 0 :
                    tofile.write("\n")

            tofile.write("%s"%"."*no_lig)
            tofile.write("*\n")

        class MyModel(automodel):
            def select_atoms(self):
                s = []
                for i in range(len(loopRanges)) :
                    s.append(selection(self.residue_range(
                        str(loopRanges[i][0]),
                        str(loopRanges[i][1])))
                    )
                return s

        a = MyModel(env, alnfile='alignment.ali', knowns = pdbIn, sequence = pdbIn+"_fill")
        a.starting_model = 1
        a.ending_model = 1
        #a.loop.md_level = refine.fast

        a.make()

        # obtain successfully modeled pdb file names
        return (self.selectBestModeller(a.outputs))

    def addLoopRefineModel(self, pdbIn, chain1,
                            fastaSeq,
                            loopRanges,
                            no_lig=0,
                            verbose=False):
        '''
        model the loops and refine them without any templates
        :param pdbIn:
        :param chain1:
        :param fastaSeq:
        :param loopRanges:
        :param no_lig:
        :param verbose:
        :return: the file name of the best model
        '''
        env = environ()
        if verbose :
            log.verbose()

        m = model(env, file=pdbIn)
        aln = alignment(env)
        aln.append_model(m, align_codes=pdbIn)
        aln.write(file=pdbIn + '.seq')

        ## add full sequence into the alignment file
        with open("alignment.ali", 'wb') as tofile:
            with open(pdbIn + '.seq') as lines:
                for s in lines :
                    if s[:-1] == "*" and s[0] != ">" :
                        tofile.write(s[:-1]+"."*no_lig+"*\n")
                    else :
                        tofile.write(s)
            tofile.write('\nP1;%s_fill \n' % pdbIn)
            tofile.write('sequence::::::::: \n')
            for i in range(len(fastaSeq)):
                tofile.write(fastaSeq[i])
                if i % 74 == 0 and i != 0:
                    tofile.write("\n")

            tofile.write("%s" % "." * no_lig)
            tofile.write("*\n")

        a = loopmodel(env, alnfile='alignment.ali', knowns = pdbIn , sequence = pdbIn+'_fill')

        a.starting_model = 1
        a.ending_model = 1

        a.loop.starting_model =  1
        a.loop.ending_model = 10
        a.loop.md_level = refine.fast

        a.make()

        # obtain successfully modeled pdb file names
        return self.selectBestModeller(a.outputs)

    def selectBestModeller(self, modellerOutput, verbose=False):
        # obtain successfully modeled pdb file names
        ok_models = [x for x in modellerOutput if x['failure'] is None]

        # Rank the models by DOPE score
        key = 'DOPE score'
        if sys.version_info[:2] == (2, 3):
            # Python 2.3's sort doesn't have a 'key' argument
            ok_models.sort(lambda a, b: cmp(a[key], b[key]))
        else:
            ok_models.sort(key=lambda a: a[key])

        # Get top model
        finalPDB = ok_models[0]
        if verbose:
            print("Top model: %s (DOPE score %.3f)" % (finalPDB['name'], finalPDB[key]))
        return (finalPDB)