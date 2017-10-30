#import pymol

import os
import numpy as np
import pandas as pd
import linecache

class RewritePDB :
    """
    Modify pdb file by changing atom indexing, resname, res sequence number and chain id
    """
    def __init__(self, filename):
        self.pdb = filename

    def pdbRewrite(self, output, chain, atomStartNdx, resStartNdx):
        """
        change atom id, residue id and chain id
        :param output: str, output file name
        :param chain: str, chain id
        :param atomStartNdx:
        :param resStartNdx:
        :return:
        """
        resseq = resStartNdx
        atomseq = int(atomStartNdx)
        chainname = chain

        newfile = open(output,'w')
        resseq_list = []

        try :
            with open(self.pdb) as lines :
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
        except FileExistsError :
            print("File %s not exist" % self.pdb)

        newfile.close()
        return 1

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

class PDB_Parser :
    """
    parse pdb file
    """
    def __init__(self, inPDB='1a28.pdb'):
        self.pdbin = inPDB

    def pdbListInfor(self, pdbList):
        #pdbInfor = defaultdict(list)
        pdbCode, resolution, ligandName, Ktype, Kvalue = [], [], [], [], []
        with open(pdbList) as lines :
            for s in lines :
                if len(s.split()) > 3 and "#" not in s :
                    pdbCode.append(s.split()[0])
                    resolution.append(s.split()[1])
                    ligandName.append(s.split()[-2])
                    Ktype.append(s.split()[-1][:2])
                    Kvalue.append(s.split()[-1][3:])
        pdbInfor = pd.DataFrame({
            "PDBCode" : pdbCode,
            "Resolution" : np.asarray(resolution ),
            "LigandName" : ligandName,
            "AffinityType" : Ktype,
            "AffinityValue": Kvalue,
        })

        return pdbInfor

    def parsePDB(self, pdbin, proteinChains, pdbOut, ligandInfor=["SUB","A"]):
        """
        subset the pdb file and only output part of the file
        :param pdbin:
        :param proteinChains:
        :param pdbOut:
        :param ligandInfor:
        :return:
        """
        pdbout = open(pdbOut, 'w')

        with open(pdbin) as lines :
            for s in lines :
                if "ATOM" in s :
                    if s[21] in proteinChains :
                        pdbout.write(s)
                elif "HETATM" in s :
                    if s[21] == ligandInfor[1] and s[17:20].strip() == ligandInfor[0] :
                        pdbout.write(s)
                else :
                    pass

        pdbout.close()
        return 1

class handlePBC :
    def __init__(self):
        pass

    def getPBCFromPBD(self, inp, ):
        """
        get the xyz pbc information from a standard pdb file,
        for lipid especially
        :param inp: str, a pdb file with CRYST1 information
        :return: list, xyz pbc range
        """
        pbcline = linecache.getline(inp, 1)
        x, y, z = pbcline.split()[1], pbcline.split()[2], pbcline.split()[3]

        return [[0.0, float(x)],
                [0.0, float(y)],
                [0.0, float(z)],
                ]

    def crdPBCRestore(self, crd, pbc):
        """
        restore xyz PBC conditions by add or minus the pbc vectors
        only cubic pbc could be handled
        :param crd: list of floats, original coordinates
        :param pbc: list of lists, pbc vectors
        :return: list of floats, the new coordinates
        """
        newcrd = []
        for i in range(3):
            if len(pbc[i]) == 2:
                if crd[i] <= pbc[i][0]:
                    newcrd.append(crd[i] + pbc[i][1])
                elif crd[i] >= pbc[i][1]:
                    newcrd.append(crd[i] - pbc[i][1])
                else:
                    newcrd.append(crd[i])
            else:
                newcrd.append(crd[i])
        return newcrd

    def checkAtomInPBCBox(self, crd, pbc):
        """
        check whether an atom in pbc box
        :param crd: list of floats, original coordinates
        :param pbc: list of lists, pbc vectors
        :return:
        """

        return all([ True for x in range(len(crd))
                     if (crd[x] <= pbc[x][1] and crd[x] >= pbc[x][0])
                     ])

    def checkResInPBCBox(self, crds, pbc, ratio):
        """
        if the ratio of atoms in box large than the given cutoff, turn True
        :param crds: list of lists, list of lists floats, coordiations
        :param pbc: 2d list, pbc ranges
        :param ratio:
        :return:
        """

        counts = map( lambda crd : all([ True for x in range(len(crd))
                                         if (crd[x] <= pbc[x][1] and crd[x] >= pbc[x][0])]), crds
                      )

        if float(sum(counts)) / float( len(crds) ) >= ratio :
            return True
        else :
            return False

class coordinatesPDB :
    def __init__(self):
        pass

    def getAtomCrdFromLines(self, lines):
        """
        given a list of atom pbd lines, return their coordinates in a 2d list
        :param lines: list of str, a list of pdb lines contains coordinates
        :return: 2d list, list of 3 element list
        """
        atomCrd = map(lambda x: [float(x[30:38].strip()),
                                 float(x[38:46].strip()),
                                 float(x[46:54].strip())],
                      lines)

        return atomCrd

    def getAtomCrdByNdx(self, singleFramePDB, atomNdx):
        """
        input a pdb file and the atom index, return the crd of the atoms
        :param singleFramePDB: file, string
        :param atomNdx: atom index, list of strings
        :return: atom coordinates, list
        """
        atomCrd = []
        with open(singleFramePDB) as lines :
            lines = [s for s in lines if len(s) > 4 and
                     s[:4] in ["ATOM","HETA"] and
                     s.split()[1] in atomNdx]
            atomCrd = map(lambda x: [float(x[30:38].strip()),
                                     float(x[38:46].strip()),
                                     float(x[46:54].strip())],
                          lines)
        return atomCrd

if __name__ == "__main__" :
    os.chdir(os.getcwd())

    pdbin = "npt_M1_control.pdb"

    pdbr = RewritePDB(pdbin)

    pdbr.pdbRewrite("new_"+pdbin, ' ', 1, 683)
