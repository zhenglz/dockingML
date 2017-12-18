# -*- coding: utf-8 -*-

import os, sys
from .gentop import GenerateTop
import numpy as np
from collections import defaultdict
import math
import urllib

class SummaryPDB :

    PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
    aminoLibFile = PROJECT_ROOT + "/data/amino-acid.lib"

    def __init__(self, pdbfile, aminoLibFile=aminoLibFile):
        self.pdbfile = pdbfile

        resShortName = {}
        try :
            with open(aminoLibFile) as lines:
                for s in lines:
                    if "#" not in s:
                        resShortName[s.split()[2]] = s.split()[3]
            self.resShortName = resShortName
        except IOError :
            self.resShortName = {}

    def netCharges(self, inputMol, ligName=None):
        '''
        Deduce the total net charge of a molecule (mostly a small ligand).
        netCharges determine the net charge of a molecule, given a pdb file.
        if the input molecule is a pdbqt (Vina input), or pqr (APBS input type),
           a molecule name (residue code) is need.
        if a mol2 file provide, only the @<TRIPOS>ATOM field will be used, no ligand name required.

        other formats are not supported.

        last column of a pdbqt, pqr and mol2 file generally will be the atomic charge field,
          otherwise, a ValueError exception will be rasied.

        :param inputMol: input file with atomic charges in the last column
        :return int netCharge
        '''
        extension = inputMol.split(".")[-1]
        netCharge = 0.0

        with open(inputMol) as lines :
            if extension in ['pdbqt', 'pqr']:
                for s in lines :

                    if s.split()[0] in ["ATOM", "HETATM"] and len(s.split()) > 5 :
                        try :
                            if ligName and ligName in s :
                                netCharge += float(s.split()[-1])
                            elif not ligName :
                                netCharge += float(s.split()[-1])
                        except ValueError :
                            netCharge += 0.0
                            print("Last column in %s is not a float point charge value."%inputMol)
            elif extension in ['mol2'] :
                condition = 0
                for s in lines :
                    if len(s.split()) and "@" in s :
                        if "<TRIPOS>ATOM" in s :
                            condition += 1
                        elif "@<TRIPOS>BOND" in s :
                            condition = 0
                        else :
                            pass
                    elif condition and len(s.split() ):
                        try:
                            netCharge += float(s.split()[-1])
                        except ValueError :
                            netCharge += 0.0
                            print("Last column in %s is not a float point charge value." % inputMol)
                #print("NET CHARGE " * 10)

            else :
                print("The file extention is not recognized. "
                      "\nPlease provide a pdbqt, pqr, or a mol2 file.")
                netCharge = 0.0

        return(int(netCharge))

    def centerOfMass(self, inputMol, atomNdx, obabelexe='obabel', molBox=False):
        '''
        Given a file (preferable PDB file format),
            if not, the file will be convert into a pdb file,
        and selected atom sequence number,
        determine the COM of the coordinates

        :param inputMol: a file, the input coordinates
        :param atomNdx: a list, atom sequence number
        :return: two lists, center of mass, format [0.0, 0.0, 0.0]
                            boxsize, format [10.0, 10.0, 10.0]
        '''
        pdb = inputMol
        if inputMol.split(".")[-1] not in ['pdb', 'pdbqt'] :
            gpdb = GenerateTop()
            gpdb.runObabel(obabelexe, inputMol, inputMol+".pdb")
            pdb = inputMol + ".pdb"
        coordinates = []

        with open(pdb) as lines :
            for s in lines :
                if "ATOM" == s.split()[0] or "HETATM" == s.split()[0] and s.split()[1] in atomNdx :
                    crd_list = []
                    crd_list.append(float(s[30:38].strip()))
                    crd_list.append(float(s[38:46].strip()))
                    crd_list.append(float(s[46:54].strip()))
                    coordinates.append(crd_list)

        coordinates = np.asarray(coordinates)

        xcenter = np.mean(coordinates[:, 0])
        ycenter = np.mean(coordinates[:, 1])
        zcenter = np.mean(coordinates[:, 2])

        if molBox :
            xsize = 2 * max(np.max(coordinates[:, 0]) - xcenter,
                        np.abs(np.min(coordinates[:,0])-xcenter))
            ysize = 2 * max(np.max(coordinates[:, 1]) - ycenter,
                        np.abs(np.min(coordinates[:, 1]) - ycenter))
            zsize = 2 * max(np.max(coordinates[:, 2]) - zcenter,
                        np.abs(np.min(coordinates[:, 2]) - zcenter))
        else :
            xsize, ysize, zsize = 100, 100, 100
        return([xcenter, ycenter, zcenter], [xsize, ysize, zsize])

    def details(self, verbose=False):
        chains = []
        resNdx = defaultdict(list)
        resName= defaultdict(list)
        resAtom= defaultdict(list)
        resNameNdx = {}
        with open(self.pdbfile) as lines :
            for s in lines :
                if len(s.split()) > 0 and s.split()[0] in ["ATOM","HETATM"] and "TER" not in s :
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

            print( "\nNumber of residues in each chain: ")
            for chain in chains :
                print( "Chain ", chain, " ",len(resNdx[chain]))

            print( "\nResidue Names for each chain are : ")
            for chain in chains :
                print( "For chain ", chain)
                for i in range(10) :
                    print( resNdx[chain][i], "  ", resName[chain][i])
                print( "......")
                for j in range(10) :
                    print( resNdx[chain][-10+j], "  ", resName[chain][-10+j])

        return chains, resNdx, resName, resAtom, resNameNdx

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
            print( "Error in function SummaryPDB.details.\nNumber of index is different with number of residues.")
            print( "Exit Now!")
            sys.exit(0)

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

        ## get protein chain sequence
        startNdx = proResNdx[0]
        finalNdx = proResNdx[-1]
        #print startNdx,finalNdx, chain, " chain starting and end index"

        fullNdx = range(startNdx, finalNdx+1)
        fullResNdx[chain] = fullNdx
        #print fullNdx
        for i in range(len(fullNdx) ):
            residueid = str(fullNdx[i]) + "_" + chain
            #print residueid
            #print i, fullNdx[i]
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
        #print chain, proteinSequence[chain]

        ## print some information
        if verbose :
            print("\nFull sequence of protein chains are: ")
            #for chain in chains :
            print( chain, "  ", proteinSequence[chain])

            print( "\nSome missing protein residues in each chain ")
            #for chain in chains :
            print( chain, "  ", missingResNdx[chain])

            print("\nThe information of the non-protein residues here ")
            #for chain in chains :
            print( chain, "  ", noProteinResName[chain])

            print( "\nThe sequence of the full protein chain is: ")
            print( "Chain  ", chain)
            sections = math.ceil(float(len(proteinSequence[chain])) / 20.0)

            for i in range(int(sections - 1)):
                print( fullResNdx[chain][i * 20],
                       " " * (-len(str(fullResNdx[chain][i * 20])) + 11),
                       fullResNdx[chain][i * 20 + 10])
                print( proteinSequence[chain][i * 20: i * 20 + 10],
                       " ",
                       proteinSequence[chain][i * 20 + 10: i * 20 + 20])
            print( fullResNdx[chain][i * 20 + 20])
            print( proteinSequence[chain][i * 20 + 20:])

        return proteinSequence, missingResNdx, noProteinResNdx, noProteinResName, fullResNdx

    def missingRes(self, chain, fastaSeq='', matchNum=10, pdbCode='', verbose=False,):
        ## find the missing residues sequences in the pdb file
        #chains, resNdx, resName, resAtom, resNameNdx = self.details()
        proteinSequence, missingResNdx, noProteinResNdx, noProteinResName, fullResNdx = self.summary(chain)
        #print fullResNdx[chain]
        trueResName = defaultdict(list)
        ## full fasta sequences here :
        fastaseq = ''
        if os.path.isfile(fastaSeq):
            with open(fastaSeq) as lines:
                for s in lines:
                    if '>' not in s :
                        #print s #strip()
                        fastaseq += s.strip()
        #print fastaseq
        elif len(fastaSeq) > 4 :
            fastaseq = fastaSeq

        else :
            fastaseq = self.downloadFasta(pdbCode, chain=chain)

        # align the fasta sequence with the protein residue sequences
        startMatch, matchedFastaSeq = proteinSequence[chain][: matchNum], fastaseq
        for i in range(len(fastaseq)-matchNum) :
            if fastaseq[i:i+matchNum] == startMatch :
                matchedFastaSeq = fastaseq[i:]
                break

        # print protein missing residue information
        startndx = fullResNdx[chain][0]
        for index in missingResNdx[chain] :
            if chain not in trueResName.keys() :
                trueResName[chain] = []
            trueResName[chain].append(matchedFastaSeq[int(index)-startndx]+'_'+index)

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
                missedRSeq += matchedFastaSeq[i]
                #print missedRSeq
            elif brokenSeq[i] != '-' and len(siteRange) > 0 :
                ranges.append([siteRange[0],siteRange[-1]])
                #missedRSeq.append(seq[i])

                missedSeqsList[str(siteRange[0])+"_"+str(siteRange[-1])] = missedRSeq
                missedRSeq = ''
                siteRange = []
            else :
                pass
        #print missedSeqsList

        if verbose :
            ## print full sequence information
            print( "The fasta sequence of the full protein chain is: ")
            print( "Chain  ", chain)
            sections = math.ceil(float(len(matchedFastaSeq) / 20.0))

            for i in range(int(sections-1)) :
                print( fullResNdx[chain][i*20],
                       " "*(-len(str(fullResNdx[chain][i*20]))+11),
                       fullResNdx[chain][i*20+10])
                print( matchedFastaSeq[i*20 : i*20+10],
                       " ",
                       matchedFastaSeq[i*20+10 : i*20+20])
            print( fullResNdx[chain][i*20+20])
            print( matchedFastaSeq[i*20+20:])

            ## print the true residue name
            print( "\nThe missing residues are :")
            print(chain,"  ", trueResName[chain])

            ## information about the missing sequences
            for ndxrange in missedSeqsList.keys() :
                print(("%12s   %-s ")%(ndxrange, missedSeqsList[ndxrange]))


        return matchedFastaSeq, missedSeqsList, fullResNdx

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
            print("\nExtend %d residues at broken part for chain %s : "%(extend,chain))
            for newrange in extMissSeq.keys() :
                print(("%12s   %-s")%(newrange, extMissSeq[newrange]))

        return extMissSeq

    def downloadFasta(self, pdbCode='', uniportID='', chain='A', writeFasta=True):
        '''
        download fasta sequence from RCSB or UniPort
        :param pdbCode:
        :param uniportID:
        :param chain:
        :param writeFasta:
        :return: the fasta sequence of a protein chain
        '''
        fastaSeq, url = '', ''
        if writeFasta :
            tofile = open(pdbCode+"_"+chain+".fasta","wb")

        if len(pdbCode) == 4 :
            try :
                url = 'http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=%s' % pdbCode
            except :
                print("Fasta Sequence of %s download failed."%pdbCode)
        elif len(uniportID)  :
            try :
                url = 'http://www.uniport.org/uniport/%s.fasta' % uniportID
            except:
                print("Fasta Sequence of %s download failed." % uniportID)
        if len(url) :
            f = urllib.urlopen(url)
            data = f.read()
            condition = False
            if len(pdbCode) == 4 :
                for item in data.split("\n") :
                    if ">" in item and item.split(":")[1][0] == chain :
                        condition = True
                    elif ">" in item and item.split(":")[1][0] != chain :
                        condition = False
                    else :
                        pass

                    if condition and '>' not in item :
                        fastaSeq += item

                    if condition and writeFasta :
                        tofile.write(item+"\n")
            else :
                tofile.write(data)
                for item in data.split("\n") :
                    if ">" not in item :
                        fastaSeq += item

            tofile.close()

        return fastaSeq