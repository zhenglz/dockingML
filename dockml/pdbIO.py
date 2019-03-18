# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
from collections import defaultdict


class rewritePDB(object):
    """
    Modify pdb file by changing atom indexing, resname, res sequence number and chain id

    Parameters
    ----------

    Attributes
    ----------

    """
    def __init__(self, inpdb):
        self.pdb = inpdb

    def pdbRewrite(self, input, output, chain, atomStartNdx, resStartNdx):
        """
        change atom id, residue id and chain id
        :param input: str, input pdb file
        :param output: str, output file name
        :param chain: str, chain id
        :param atomStartNdx: int,
        :param resStartNdx: int
        :return:
        """
        resseq = int(resStartNdx)
        atomseq = int(atomStartNdx)
        chainname = chain

        newfile = open(output,'w')
        resseq_list = []

        try :
            with open(input) as lines :
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
        resseqstring = " "*(4 - len(str(resseq)))+str(resseq)
        newline = inline[:22] + resseqstring + inline[26:]
        return newline

    def atomSeqChanger(self, inline, atomseq):
        atomseqstring = " " * (5 - len(str(atomseq))) + str(atomseq)
        newline = inline[:6] + atomseqstring + inline[11:]
        return newline

    def resNameChanger(self, inline, resname):
        resnamestr = " " * (4 - len(str(resname))) + str(resname)
        newline = inline[:16] + resnamestr + inline[20:]
        return newline

    def chainIDChanger(self, inline, chainid) :
        newline = inline[:21] + str(chainid) + inline[22:]
        return newline

    def atomNameChanger(self, inline, new_atom_name):
        newline = inline[:12] + "%4s" % new_atom_name + inline[16:]
        return newline

    def combinePDBFromLines(self, output, lines):
        """
        combine a list of lines to a pdb file

        Parameters
        ----------
        output
        lines

        Returns
        -------

        """

        with open(output, "wb") as tofile :
            tmp = map(lambda x: tofile.write(x), lines)

        return 1

    def swampPDB(self, input, atomseq_pdb, out_pdb, chain="B"):
        """
        given a pdb file (with coordinates in a protein pocket), but with wrong atom
        sequence order, try to re-order the pdb for amber topology building

        Parameters
        ----------
        input:str,
            the pdb file with the correct coordinates
        atomseq_pdb:str,
            the pdb file with correct atom sequences
        out_pdb: str,
            output pdb file name
        chain: str, default is B
            the chain identifier of a molecule

        Returns
        -------

        """

        tofile = open("temp.pdb", 'w')

        crd_list = {}

        ln_target, ln_source = 0, 0
        # generate a dict { atomname: pdbline}
        with open(input) as lines :
            for s in [x for x in lines if ("ATOM" in x or "HETATM" in x)]:
                crd_list[s.split()[2]] = s
                ln_source += 1

        # reorder the crd_pdb pdblines, according to the atomseq_pdb lines
        with open(atomseq_pdb) as lines:
            for s in [x for x in lines if ("ATOM" in x or "HETATM" in x)]:
                newline = crd_list[s.split()[2]]
                tofile.write(newline)
                ln_target += 1

        tofile.close()

        if ln_source != ln_target:
            print("Error: Number of lines in source and target pdb files are not equal. (%s %s)"%(input, atomseq_pdb))

        # re-sequence the atom index
        self.pdbRewrite(input="temp.pdb", atomStartNdx=1, chain=chain, output=out_pdb, resStartNdx=1)

        os.remove("temp.pdb")

        return None


class parsePDB(object):
    """
    A pdb file parser.

    Parameters
    ----------
    inPDB : str, default = '1a28.pdb'

    Attributes
    ----------
    resSideChainAtoms : list
        The sidechain atoms of all natural amino acids
    resMainchainAtoms : list
        The mainchain atoms of all natural amino acids
    prores : str
        The amino acid standard names
    nucleic : str
        The nucleic acid standard atom names

    """
    def __init__(self, inPDB='1a28.pdb'):
        self.pdbin = inPDB

        self.resSideChainAtoms = ['CZ2', 'OE2', 'OE1', 'OG1', 'CD1', 'CD2', 'CG2', 'NE', 'NZ', 'OD1',
                                  'ND1', 'ND2', 'OD2', 'CB', 'CZ3', 'CG', 'CZ',
                                  'NH1', 'CE', 'CE1', 'NH2', 'CG1', 'CD', 'OH', 'OG', 'SG', 'CH2',
                                  'NE1', 'CE3', 'SD', 'NE2', 'CE2']
        self.resMainChainAtoms = ['CA', 'N', 'C', 'O']

        self.prores = "amino-acid.lib"
        if not os.path.exists("amino-acid.lib") :
            PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
            self.prores = PROJECT_ROOT + '/../data/amino-acid.lib'

        self.nucleic = "nucleic-acid.lib"
        if not os.path.exists("nucleic-acid.lib") :
            PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
            self.nucleic = PROJECT_ROOT + '/../data/nucleic-acid.lib'

    def readDomainRes(self, filein):
        """Load the domain data file, return a list of domain residue information

        Parameters
        ----------
        filein : str,
            Input file name

        Returns
        -------
        drange: list
            List of list (>3 elements per item),
            domain-residue-range information, example
            [ [ 'HNH', 1, 10], ['RuvC', 11, 92, 150, 190], ['LBD', 93, 149]]

        """

        drange = []

        try:
            with open(filein) as lines:
                for s in lines:
                    if "#" not in s and len(s.split()) >= 3:
                        d = [s.split()[0]]
                        d += [int(x) for x in s.split()[1:]]
                        drange.append(d)
        except FileNotFoundError:
            print("Either the file does not exist or there are empty lines in the file. ")

        return drange

    def getStdProRes(self):
        """Get the standard protein residue list

        Parameters
        ----------
        resname : list
            standard residue name list
        """

        with open(self.prores) as lines :
            lines = [ x for x in lines if "#" not in x ]

            resname = [ x.split()[2] for x in lines ]

        return resname

    def shortRes2LongRes(self):
        """Convert the short single-character residue name to long 3-code name

        Parameters
        ----------
        resmap: dict,
            { shortname: longname}
        """

        resmap = {}
        with open(self.prores) as lines:
            lines = [x for x in lines if "#" not in x]

            for s in lines :
                resmap[s.split()[3]] = s.split()[2]

        return resmap

    def longRes2ShortRes(self):
        """Convert the long 3-code name to short single-character residue name

        Parameters
        ----------
        resmap: dict, { longname: shortname }
        """

        resmap = {}
        with open(self.prores) as lines:
            lines = [x for x in lines if "#" not in x]

            for s in lines:
                resmap[s.split()[2]] = s.split()[3]

        return resmap

    def pdbListInfor(self, pdbList):
        """Get a list of information from a rcsb database pdb file

        Parameters
        ----------
        pdbInfor: pd.DataFrame,
            the pdb information for each atom
        """
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
            "PDBCode": pdbCode,
            "Resolution": np.asarray(resolution),
            "LigandName": ligandName,
            "AffinityType": Ktype,
            "AffinityValue": Kvalue,
        })

        return pdbInfor

    def subsetPDB(self, pdbin, proteinChains, pdbOut, ligandInfor=["SUB","A"]):
        """Subset the pdb file and only output part of the file

        Parameters
        ----------
        pdbin : str
            Input pdb file name
        proteinChains : list
            The selected chain identifiers
        pdbOut : str
            The output file name
        ligandInfor : list
            The resname of the ligand and its chain identifier

        Returns
        -------
        self : an instance of itself

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
                else:
                    pass

        pdbout.close()
        return self

    def withSubGroup(self, isProtein=True):
        """
        check whether a pdb line is protein or nucleic acids
        :param isProtein:
        :return:
        """

        if isProtein:
            subgroup = {}

            for atom in self.resMainChainAtoms :
                subgroup[atom] = "mainchain"
            for atom in self.resSideChainAtoms :
                subgroup[atom] = 'sidechain'

            return subgroup
        else:
            xna = {}
            try:
                with open(self.nucleic) as lines:
                    for s in [ x for x in lines.readlines() if "#" not in x ] :
                        xna[s.split()[-1]] = s.split()[1]

            except FileNotFoundError:
                print("File %s not exist" % self.nucleic)

            return xna

    def atomInformation(self, pdbin):
        """
          # elements in atom infor
          # key: str, atom index
          # value: [atomname, molecule type, is_hydrogen,  resname, resndx, chainid,(mainchian, sidechain,
          #           sugar ring, phosphate group, base ring), element]
          :param pdbin:
          :return: a dictionary of list, key is atom ndx
        """

        atominfor = defaultdict(list)

        if not os.path.exists(self.prores):
            PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
            proteinres = PROJECT_ROOT + '/../data/amino-acid.lib'

        with open(self.prores) as lines:
            protRes = [s.split()[2] for s in lines if "#" not in s]
        DNARes = ['DA', 'DT', 'DC', 'DG']
        RNARes = ['A', 'G', 'C', 'U']

        prosubgroup = self.withSubGroup(True)
        xnasubgroup = self.withSubGroup(False)

        with open(pdbin) as lines:
            for s in lines:
                if len(s.split()) and s[:4] in ["ATOM", "HETA"]:

                    atomndx = s[6:11].strip()
                    atomname = s[12:16].strip()
                    resname = s[17:20].strip()
                    resndx = int(s[22:26].strip())
                    chain = s[21]
                    if len(s) > 76:
                        if s[77] != " ":
                            element = s[77]
                        else:
                            element = s[13]

                        if s[13] == "B" :
                            element = "Br"
                        if "Cl" in s.split()[1] or "CL" in s.split()[1] :
                            element = "Cl"

                    else:
                        element = s[13]

                    hydrogen = {"H": True}
                    is_hydrogen = hydrogen.get(s[13], False)

                    if resname in protRes:
                        moltype = "Protein"
                    elif resname in DNARes:
                        moltype = "DNA"
                    elif resname in RNARes:
                        moltype = "RNA"
                    else:
                        moltype = "Unknown"

                    if moltype == "Protein":
                        subgroup = prosubgroup.get(atomname, "Unknown")
                    elif moltype in ["RNA", "DNA"]:
                        subgroup = xnasubgroup.get(atomname, "Unknown")
                    else:
                        subgroup = "Unknown"

                    atominfor[atomndx] = [atomname, moltype, is_hydrogen, resname, resndx, chain, subgroup, element]
                else:
                    pass

        return atominfor

    def getResNamesList(self, pdbin, chains):
        """get a list of residue names in specific chains

        Parameters
        ----------
        pdbin: str,
            input pdb file name
        chains: list of str,
            chain identifiers

        Returns
        -------
        resname: list
            the list of residue names
        """

        with open(pdbin) as lines:
            lines = [x for x in lines if (x[:4] in ["ATOM", "HETA"] and x[21] in chains)]

            resname = []
            for s in lines:
                if s.split()[3] not in resname:
                    resname.append(s.split()[3])

        return resname

    def getNdxForRes(self, pdbin, chains):
        """get the residue names, as well as seq ndx, and chain id

        :param pdbin:
        :param chains: list or str
        :return: a list of sets, [ (resname, seqid, chainid), ]
        """

        with open(pdbin) as lines:
            lines = [x for x in lines if (x[:4] in ["ATOM", "HETA"]
                                          and x[21] in chains)]

            reslist = []
            for s in lines:
                id = (s.split()[3], s[22:26].strip(), s[21])
                if id not in reslist:
                    reslist.append(id)

        return reslist

    def getNdxForMol(self, pdbin, molres=[], resndx=[]):
        """
        generally, a molecule is composed by a list of atoms
        these atoms' sequence number is a list of intergers (atomndx)
        :param pdbin:
        :param molres:
        :param resndx:
        :return:
        """

        atominf = self.atomInformation(pdbin)
        # ndx in atominf keys are string
        atomndx = [str(x) for x in sorted([x for x in atominf.keys()])]

        return atomndx


class handlePBC(object):
    def __init__(self):
        pass

    def getPBCFromPBD(self, inp, ):
        """
        get the xyz pbc information from a standard pdb file,
        for lipid especially
        :param inp: str, a pdb file with CRYST1 information
        :return: list, xyz pbc range
        """
        with open(inp) as lines :
            pbcline = [ x for x in lines if "CRYST1" in x ][0]

        x, y, z = pbcline.split()[1], pbcline.split()[2], pbcline.split()[3]

        return [
                [0.0, float(x)],
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


class coordinatesPDB(object):
    def __init__(self):
        pass

    def replaceCrdInPdbLine(self, line, newxyz):
        """Input a line of pdb file, and a new vector of xyz values,
        the old xyz values will be replaces, and return a new pdb line

        Parameters
        ----------
        line : str
            A line from pdb file
        newxyz : list, shape = [3, ]
            The new xyz coordinates, in unit nanometer

        Returns
        -------
        new_line : str
            The new PDB line with new xyz coordinates

        """

        if "ATOM" in line or "HETATM" in line:
            head = line[:30]
            tail = line[54:]

            newline = head + "{0:8.3f}{1:8.3f}{2:8.3f}".format(newxyz[0], newxyz[1], newxyz[2]) + tail

        else :
            print("WARNING: %s is not a coordination line"%line)
            newline = ""

        return newline

    def getAtomCrdFromLines(self, lines):
        """Given a list of atom pbd lines, return their coordinates in a 2d list

        Parameters
        ----------
        lines : list
            A list of pdb lines contains coordinates

        Returns
        -------
        coordinates : list, shape = [ N, 3]
            The coordinates of selected pdb lines, N is the number of
            lines.
        """

        atomCrd = list(map(lambda x: [float(x[30:38].strip()),float(x[38:46].strip()),
                                      float(x[46:54].strip())],lines))

        return atomCrd

    def getAtomCrdByNdx(self, singleFramePDB, atomNdx=['1',]):
        """Input a pdb file and the atom index, return the crd of the atoms

        Parameters
        ----------
        singleFramePDB : str
            Input pdb file name.
        atomNdx : list
            A lis of atoms for coordinates extraction.

        Returns
        -------
        coordinates : list, shape = [ N, 3]
            The coordinates for the selected lines. N is the number of
            selected atom indices.

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
        return list(atomCrd)


if __name__ == "__main__" :
    os.chdir(os.getcwd())

    pdbin = "npt_M1_control.pdb"

    pdbr = rewritePDB(pdbin)

    pdbr.pdbRewrite( pdbin, "new_"+pdbin, ' ', 1, 683)

