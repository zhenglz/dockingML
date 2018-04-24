#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Indexing PDB file
Generate Gromacs Format Index File
"""

import sys, os
from collections import defaultdict
from collections import OrderedDict
import argparse
from argparse import RawTextHelpFormatter
import linecache

class PdbIndex :
    '''
    Input the reisude number sequence, then out put the required index atoms
    '''
    def __init__(self, ) :
        self.backbone  = ['CA', 'N']
        self.mainchain = ['CA', 'N', 'C', 'O']
        self.ca        = ['CA']
        self.phi       = ['C', 'N', 'CA', 'C']
        self.psi       = ['N', 'CA', 'C', 'N']

    def res_index(self, inpdb, chain, atomtype, residueNdx, atomList, dihedraltype="None"):
        '''
        Obtain atom index from a reference pdb file
        provide information including: residue indexing, atomtype, atomname

        :param inpdb:
        :param chain:
        :param atomtype: str, options: all-atom, non-hydrogen, dihedral
        :param residueNdx: list, dimension 2*1
        :param atomList: list, explicit atom name list
        :param dihedraltype:
        :param atomName:
        :return: a list, of atom index
        '''

        # atomInfor = {}
        indexlist = []
        resSeqNdx = []
        if len(residueNdx) == 1 :
            resSeqNdx = range(residueNdx[0], residueNdx[0] +1)
        elif len(residueNdx) >= 2 and len(residueNdx) % 2 == 0 :
            for k in range(int(len(residueNdx)/2)) :
                resSeqNdx += range(residueNdx[k*2], residueNdx[k*2+1]+1)
        else :
            print("Error!! No residue index provided. ")
            resSeqNdx = []

        if atomtype == "dihedral" :
            indexlist = []
            for resndx in resSeqNdx :
                phitype = self.phi
                phipair = [-1, -1, -1, -1]
                phindxadd = [-1, 0, 0, 0]
                for i in range(4) :
                    with open(inpdb) as lines :
                        for s in lines :
                            if len(s.split()) > 1 and s.split()[0] == "ATOM" \
                                    and s[21] == chain and int(s[22:26].strip()) == resndx + phindxadd[i] \
                                    and s[12:16].strip() == phitype[i] :
                                phipair[i] = int(s.split()[1])

                psitype = self.psi
                psindxadd = [ 0, 0, 0, 1]
                psipair = [-1, -1, -1, -1]
                for i in range(4):
                    with open(inpdb) as lines:
                        for s in lines:
                            if len(s.split()) > 1 and s.split()[0] == "ATOM" and s[21] == chain and \
                                            int(s[22:26].strip()) == resndx + psindxadd[i] and\
                                            s[12:16].strip() == psitype[i]:
                                psipair[i] = int(s.split()[1])
                                # break

                if "PHI" in dihedraltype :
                    indexlist.append(psipair)
                if "PSI" in dihedraltype :
                    indexlist.append(phipair)
        else :
            with open(inpdb) as lines :
                for s in lines :
                    if len(s.split()) > 1 and s.split()[0] == "ATOM" and s[21] == chain and \
                                    int(s[22:26].strip()) in resSeqNdx :
                        if atomtype == "non-hydrogen" :
                            if s[13] != "H" and s.split()[2][0] != "H" and s.split()[-1] != "H" :
                                indexlist.append(s.split()[1])
                        elif atomtype == "all-atom" :
                            ## all atoms
                            indexlist.append(s.split()[1])
                        elif atomtype == "side-chain-noH" :
                            if s.split()[3] != "GLY" :
                                if s[12:16].strip() not in self.mainchain \
                                        and s[13] != "H" and s.split()[-1] != "H":
                                    indexlist.append(s.split()[1])
                            else :
                                if s.split()[2] == "H" :
                                    indexlist.append(s.split()[1])
                        elif atomtype == "side-chain" :
                            if s[12:16].strip() not in self.mainchain :
                                indexlist.append(s.split()[1])
                        else:
                            if s[12:16].strip() in atomList :
                                indexlist.append(s.split()[1])
        return(indexlist)

    def atomList(self, atomtype, atomname):
        '''
        provide information of atom type and atomname
        :param atomtype:
        :param atomname:
        :return:
        '''
        atomList = []
        if atomname :
            atomList = atomname
        else :
            if "mainchain" in atomtype or "Main" in atomtype or "main" in atomtype :
                atomList = self.mainchain

            elif "CA" in atomtype or "ca" in atomtype or "Ca" in atomtype or "alpha" in atomtype :
                atomList = self.ca

            elif "backbone" in atomtype or "Back" in atomtype or "bone" in atomtype :
                atomList = self.backbone

            elif "all" in atomtype :
                atomtype = "all-atom"

            elif "no hy" in atomtype or "non-hy" in atomtype :
                atomtype = "non-hydrogen"

            elif "side" in atomtype or "Sidechain" in atomtype or "sidechain" in atomtype :
                atomtype = "side-chain"

            elif "PSI" in atomtype or "PHI" in atomtype or "phi" in atomtype or 'psi' in atomtype :
                atomtype = "dihedral"

        return(atomList , atomtype)

    def atomList2File(self, atomNdxList, groupName, outputNdxFile, append=True):
        if append :
            tofile = open(outputNdxFile, 'a')
        else :
            tofile = open(outputNdxFile, 'wb')

        tofile.write('[ %s ] \n'% groupName )

        for i, atom in enumerate(atomNdxList) :
            tofile.write('%6d ' % int(atom))
            if i % 15 == 0:
                tofile.write('  \n')

        tofile.write("\n")
        tofile.close()
        return 1

    def withSubGroup(self, isProtein=True, nucleic="nucleic-acid.lib"):

        if not os.path.exists(nucleic):
            PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
            nucleic = os.path.join(PROJECT_ROOT, '/../data/nucleic-acid.lib')

        if isProtein:
            subgroup = {}
            for atom in ['CA', 'N', 'C', 'O']:
                subgroup[atom] = "mainchain"
            for atom in ['CZ2', 'OE2', 'OE1', 'OG1', 'CD1', 'CD2', 'CG2', 'NE', 'NZ', 'OD1',
                         'ND1', 'ND2', 'OD2', 'CB', 'CZ3', 'CG', 'CZ',
                         'NH1', 'CE', 'CE1', 'NH2', 'CG1', 'CD', 'OH', 'OG', 'SG', 'CH2',
                         'NE1', 'CE3', 'SD', 'NE2', 'CE2']:
                subgroup[atom] = 'sidechain'

            return subgroup
        else:
            xna = {}
            try :
                with open(nucleic) as lines:
                    for s in lines:
                        if "#" not in s:
                            xna[s.split()[-1]] = s.split()[1]

            except FileNotFoundError :
                print("File %s not exist" % nucleic)
                xna = {}
            return xna

    def atomInformation(self, pdbin, proteinres="amino-acid.lib"):
        # elements in atom infor
        # key: atom index
        # value: [atomname, molecule type, is_hydrogen,  resname, resndx, chainid,(mainchian, sidechain,
        #           sugar ring, phosphate group, base ring)]
        atominfor = defaultdict(list)

        if not os.path.exists(proteinres) :
            PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
            proteinres = PROJECT_ROOT + '/../data/amino-acid.lib'

        with open(proteinres) as lines:
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

    def arguements(self) :
        d ='''
        ################################################################
        # Generate GMX Index from a PDB file
        # Generate POSRES file for a PDB File or Gro File
        # Contact  LZHENG002@E.NTU.EDU.SG
        # Version  2.2
        # Update Mar 23, 2017
        ################################################################

        Usage examples

        Generate backbone index for residue 1 to 1000 within Chain B
        python autoMD.py index -res 1 100 -at backbone -chain B -out index.ndx

        Generate all atom index for residue 78 to 100
        python autoMD.py index -inp S_1.pdb -out index.ndx -at allatom -chain ' ' -res 78 100

        Generate dihedral (PHI only) quadroplex index
        python autoMD.py index -inp S_1.pdb -out index.ndx -at dihedral -chain ' ' -res 78 100 -dihedral PHI

        '''

        parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
        # parser.add_argument('-h','--help', help="Show this help information. \n")
        parser.add_argument('-pdb', '--pdbfile',type=str, default='input.pdb',
                            help="Input PDB file for generating Index. ")
        parser.add_argument('-out', '--output',type=str, default='output.ndx',
                            help="Output Index file including atoms sequence number.\n"
                                 "Default name is output.ndx \n")
        parser.add_argument('-at', '--atomtype',type=str, default ='allatom',
                            help="Selected atom type for generating index. \n"
                                 "Options including: allatom, mainchain, \n"
                                 "non-hydrogen, c-alpha, backbone, sidechain, dihedral\n"
                                 "Default choice is: allatom \n")
        parser.add_argument('-an','--atomname', type=str, default=[], nargs='+',
                            help="Select atom by atom names. A list of atom names \n"
                                 "could be supplied. If not given, all atom names are \n"
                                 "considered. Default is [].")
        parser.add_argument('-chain', '--chainId',type=str, default= "A",
                            help="Protein chain identifier. Default chain ID is A. \n")
        parser.add_argument('-res', '--residueRange',type=int, nargs= '+',
                            help="Residue sequence number for index generating. \n"
                                 "Example, -res 1 100, generateing atom index within \n"
                                 "residues 1 to 100. Default is None.")
        parser.add_argument('-posres', '--positionRestraint',default=False,
                            help="Generate a postion restraint file for selected index.\n"
                                 "Default name is posres.itp \n")
        parser.add_argument('-dihe', '--dihedralType',default=None, type =str,
                            help="Whether generate dihedral angles index (quadruplex).\n"
                                 "Phi and Psi are considered. Optional choices are: \n"
                                 "PHI, PSI, PHI_PSI, or NA. Default is NA. \n")
        parser.add_argument('-gn','--groupName', type=str, default=None,
                            help="The name of the group of atoms selected. \n"
                                 "Default is None.")
        parser.add_argument('-append', '--appendFile', default=True, type=bool,
                            help="Append the group of index to the output file. \n"
                                 "Options: True, False. \n"
                                 "Default is True.")

        args, unknown = parser.parse_known_args()

        # decide to print help message
        if len(sys.argv) < 3 :# no enough arguements, exit now
            parser.print_help()
            print("\nYou chose non of the arguement!\nDo nothing and exit now!\n")
            sys.exit(1)

        return(args)

    def genGMXIndex( self):
        '''
        run geneIndex
        initialize argument parser function.
        :return:
        '''
        args = self.arguements()

        if len(args.residueRange) == 2:
            residueNdx = range(args.residueRange[0], args.residueRange[-1] + 1)
        elif len(args.residueRange) == 1:
            residueNdx = range(args.residueRange[0], args.residueRange[0] + 1)
        elif len(args.residueRange) > 2:
            residueNdx = args.residueRange
        else:
            print("\nNumber of residue id is not correct. \nExit Now. ")
            sys.exit(1)

        if args.dihedralType :

            atomlist,atomtype = self.atomList(args.dihedralType, args.atomname)
        else :
            atomlist,atomtype = self.atomList(args.atomtype, args.atomname)

        ## generate index
        atomndx = self.res_index(inpdb=args.pdbfile, chain=args.chainId,
                                 residueNdx= args.residueRange, atomList=atomlist,
                                 atomtype=atomtype, dihedraltype=args.dihedralType
                                 )
        if args.appendFile :
            append = 'a'
        else :
            append = 'wb'
        tofile = open(args.output, append)

        if args.groupName :
            tofile.write('[ %s ] \n' % (args.groupName.strip()))
        else :
            tofile.write('[ %s_%d_%d ] \n' % (args.chainId, args.residueRange[0], args.residueRange[-1]))

        if atomtype == "dihedral":
            for index in atomndx:
                for i in range(4):
                    tofile.write("%5d " % index[i])
                tofile.write("  \n")

        else:
            i = 0
            for atom in atomndx:
                i += 1
                tofile.write('%4s ' % atom)
                if i % 15 == 0:
                    tofile.write('  \n')
            tofile.write(" \n")
            tofile.close()

            if args.positionRestraint:
                with open('posres.itp', 'w') as tofile :
                    tofile.write("[ position_restraints ] \n; ai  funct  fcx    fcy    fcz  \n")
                    for atom in atomndx:
                        tofile.write("%12d  1  1000  1000  1000  \n" % int(atom))

        print("\nGenerating Index File Completed!")
        return(1)

class GmxIndex :

    def __init__(self, index):
        if os.path.exists(index) :
            self.index = index
        else:
            print("File {} not exists!".format(index))
            sys.exit(0)

        self.ndxlines   = open(self.index).readlines()
        self.groups     = [ x.split()[1] for x in self.ndxlines if ("[" in x and "]" in x) ]
        self.totallines = len(self.ndxlines)

    def groupsLineNumber(self):

        groupLN = OrderedDict()

        linecount = 0
        for s in self.ndxlines :

            if "[" in s and "]" in s :
                groupLN[s.split()[1]] = linecount

            linecount += 1

        return groupLN

    def groupContent(self, group):

        gln = self.groupsLineNumber()

        start_ln = gln[group] + 1
        end_ln = 0

        if self.groups.index(group) == len(self.groups) - 1 :
            end_ln = self.totallines - 1
        else :
            end_ln = gln[self.groups[self.groups.index(group) + 1]] - 1

        contents = [ self.ndxlines[x].strip("\n")+" " for x in range(start_ln, end_ln+1) ]

        contents = " ".join(contents)

        return contents.split()

    def writeNdxGroup(self, group, elements, output="output.ndx"):

        PdbIndex().atomList2File(elements, group, output)

        return 1

def main() :

    ndx = PdbIndex()

    ndx.genGMXIndex()

    return 1