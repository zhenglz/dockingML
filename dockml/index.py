#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Indexing PDB file
Generate Gromacs Format Index File
"""

import sys, os
from collections import OrderedDict
from mdanaly import gmxcli
import mdtraj as mt


class PdbIndex(object):
    """Input the residue number sequence, then out put the required index atoms

    Attributes
    ----------
    backbone: list,
        residue backbone atom name list
    mainchain: list,
        residue mainchain atom name list
    ca: list,
        residue alpha carbon atom list
    phi: list,
        residue phi torsion angle atom name list
    psi: list
        residue psi torsion angle atom name list

    Examples
    --------
    >>> from dockml import index
    >>> ndx = index.PdbIndex()
    >>> ndx
    <dockml.index.PdbIndex object at 0x2b704a297860>
    >>> ndx.res_index("./test/input.pdb", "C", "all-atom", residueNdx=[1, 2], atomList=[])
    ['22378', '22379', '22380', '22381', '22382', '22383', '22384', '22385', '22386', '22387', '22388', '22389', '22390',
    '22391', '22392', '22393', '22394', '22395', '22396', '22397', '22398', '22399', '22400', '22401', '22402', '22403',
    '22404', '22406', '22407', '22408', '22409', '22410', '22411', '22412', '22413', '22414', '22415', '22416', '22417',
    '22418', '22419', '22420', '22421', '22422', '22423', '22424', '22425', '22426', '22427', '22428', '22429', '22430',
    '22431', '22432', '22433', '22434', '22435', '22436', '22437']

    See Also
    --------
    GmxIndex
        parse gromacs index file

    """

    def __init__(self, reference="input.pdb",
                 chain=["A"], resSeq=[0, -1],
                 atomtype=["CA"], dihedral=["NA"]):

        self.phi = ['C', 'N', 'CA', 'C']
        self.psi = ['N', 'CA', 'C', 'N']
        self.at_direct = ["all", "none", "sidechain", "mainchain", "protein", "water", ]

        self.dihedral = dihedral
        self.pdb = reference
        self.top = None

        # chain identifiers, list
        self.chain = chain
        self.resSeq = resSeq
        self.at = atomtype
        self.atomndx = None
        #self.molecule = None

    def load_pdb(self):

        if not os.path.exists(self.pdb) :
            print("PDB reference file %s is not existed, exit now!" % self.pdb)
            sys.exit(0)

        if self.pdb[-4:] != ".pdb":
            print("You just loaded a file %s \nPlease provide a reference pdb file."%self.pdb)

        traj = mt.load(self.pdb)
        self.top = traj.topology

        return self

    def res_seq(self):
        if len(self.resSeq) == 1:
            self.resSeq = [self.resSeq[0], self.resSeq[0]]

        if self.resSeq[1] == -1:
            self.resSeq[1] = self.top.n_residues - self.resSeq[0] + 1

        return self

    def res_index(self):
        if self.dihedral[0] != "NA":
            indexlist = []
            for resndx in range(self.resSeq[0], self.resSeq[-1]+1):
                phitype = self.phi
                phipair = [-1, -1, -1, -1]
                phindxadd = [-1, 0, 0, 0]
                for i in range(4):
                    with self.pdb as lines:
                        for s in lines:
                            if len(s.split()) > 1 and s.split()[0] == "ATOM" \
                                    and s[21] in self.chain \
                                    and int(s[22:26].strip()) == resndx + phindxadd[i] \
                                    and s[12:16].strip() == phitype[i]:
                                phipair[i] = int(s.split()[1])

                psitype = self.psi
                psindxadd = [0, 0, 0, 1]
                psipair = [-1, -1, -1, -1]
                for i in range(4):
                    with self.pdb as lines:
                        for s in lines:
                            if len(s.split()) > 1 and s.split()[0] == "ATOM" \
                                    and s[21] in self.chain and \
                                    int(s[22:26].strip()) == resndx + psindxadd[i] \
                                    and s[12:16].strip() == psitype[i]:
                                psipair[i] = int(s.split()[1])
                                # break

                if "PSI" in self.dihedral:
                    psipair = [x for x in psipair if all(x) > 0]
                    indexlist.append(psipair)
                if "PHI" in self.dihedral:
                    phipair = [x for x in phipair if all(x) > 0]
                    indexlist.append(phipair)

        else:

            if self.at in self.at_direct:
                names = self.at + " and"
            elif self.at in ["CA", "ca"]:
                names = "name %s and " % self.at
            else:
                names = ""

            chains = " and ".join(["chain " + x for x in self.chain])
            resndx = " and resid %d to %d" % (self.resSeq[0], self.resSeq[1])

            selections = names + resndx + chains

            indexlist = self.top.select(selections)

        return indexlist

    '''def res_index(self, inpdb, chain, atomtype, residueNdx, atomList, dihedraltype=["None"]):
        """Obtain atom index from a reference pdb file
        provide information including: residue indexing, atomtype, atomname

        Parameters
        ----------
        inpdb : str,
            the input pdb file name
        chain : str,
            the chain identifier in a pdb file
        atomtype : str, options: all-atom, non-hydrogen, dihedral
            the atom type for index generation
        residueNdx: iterable, shape = [2]
            the residue index [ start, end ]
        atomList : list,
            explicit atom name list
        dihedraltype : str, options= [ PHI, PSI, PHI_PSI]
            the dihedral angle type

        Returns
        -------
        indexlist : list,
            a list of elements (integers) in certain group

        """

        # atomInfor = {}
        indexlist = []
        resSeqNdx = []
        if len(residueNdx) == 1 :
            resSeqNdx = range(residueNdx[0], residueNdx[0] +1)
        elif len(residueNdx) >= 2 and len(residueNdx) % 2 == 0:
            for k in range(int(len(residueNdx)/2)):
                resSeqNdx += range(residueNdx[k*2], residueNdx[k*2+1]+1)
        else:
            print("Error!! No residue index provided. ")
            resSeqNdx = []

        if atomtype == "dihedral":
            indexlist = []
            for resndx in resSeqNdx:
                phitype = self.phi
                phipair = [-1, -1, -1, -1]
                phindxadd = [-1, 0, 0, 0]
                for i in range(4):
                    with open(inpdb) as lines :
                        for s in lines :
                            if len(s.split()) > 1 and s.split()[0] == "ATOM" \
                                    and s[21] == chain and int(s[22:26].strip()) == resndx + phindxadd[i] \
                                    and s[12:16].strip() == phitype[i]:
                                phipair[i] = int(s.split()[1])

                psitype = self.psi
                psindxadd = [0, 0, 0, 1]
                psipair = [-1, -1, -1, -1]
                for i in range(4):
                    with open(inpdb) as lines:
                        for s in lines:
                            if len(s.split()) > 1 and s.split()[0] == "ATOM" and s[21] == chain and \
                                            int(s[22:26].strip()) == resndx + psindxadd[i] and\
                                            s[12:16].strip() == psitype[i]:
                                psipair[i] = int(s.split()[1])
                                # break

                if "PHI" in dihedraltype:
                    indexlist.append(psipair)
                if "PSI" in dihedraltype:
                    indexlist.append(phipair)
        else:
            with open(inpdb) as lines:
                for s in lines:
                    if len(s.split()) > 1 and s.split()[0] == "ATOM" and s[21] == chain and \
                                    int(s[22:26].strip()) in resSeqNdx:
                        if atomtype == "non-hydrogen":
                            if s[13] != "H" and s.split()[2][0] != "H" and s.split()[-1] != "H" :
                                indexlist.append(s.split()[1])
                        elif atomtype == "all-atom":
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
                            if s[12:16].strip() in atomList:
                                indexlist.append(s.split()[1])
        return indexlist '''

    '''def atomList(self, atomtype, atomname):
        """Provide information of atom type and atomname

        Parameters
        ----------
        atomtype : str,
            the type of atoms for groupping
        atomname : list, iterable, or array
            a list of explict atom names

        Returns
        -------
        atomList : list, or iterable, or array
            a list of atom names
        atomtype : str,
            the atom types for groupping

        """
        atomList = []
        if atomname :
            atomList = atomname
        else :
            if "mainchain" in atomtype or "Main" in atomtype or "main" in atomtype:
                atomList = self.mainchain

            elif "CA" in atomtype or "ca" in atomtype or "Ca" in atomtype or "alpha" in atomtype:
                atomList = self.ca

            elif "backbone" in atomtype or "Back" in atomtype or "bone" in atomtype:
                atomList = self.backbone

            elif "all" in atomtype:
                atomtype = "all-atom"

            elif "no hy" in atomtype or "non-hy" in atomtype or "heavy":
                atomtype = "non-hydrogen"

            elif "side" in atomtype or "Sidechain" in atomtype or "sidechain" in atomtype:
                atomtype = "side-chain"

            elif "PSI" in atomtype or "PHI" in atomtype or "phi" in atomtype or 'psi' in atomtype:
                atomtype = "dihedral"

        return atomList, atomtype'''

    def atomList2File(self, atom_list, group_name, write_dihe=False,
                      out_filen="output.ndx", append=True):
        """
        Save a group of atom index into a gromacs-format index file

        Parameters
        ----------
        atom_list : list, or array, or iterable, or sequence
            a list of atom index to be written into the output file
            if write_dihe = True, atom_list is a 2d list or np.ndarry
        group_name : str,
            the group name to be written into the output file
        write_dihe : bool, default = False
            whether write dihedrals to output file
        out_filen : str, default = 'output.ndx'
            the output file name
        append : bool, default = True
            append the content into the end of the output file,
            if it exists

        Returns
        -------
        self : returns an instance of self.
        """
        if append:
            tofile = open(out_filen, 'a')
        else:
            tofile = open(out_filen, 'wb')

        tofile.write('[ %s ] \n'% group_name)

        if not write_dihe:
            for i, atom in enumerate(atom_list):
                tofile.write('%6d ' % int(atom))
                if i % 15 == 0:
                    tofile.write('  \n')
            tofile.write("\n")

        else:
            for dihe in atom_list:
                tofile.write("   ".join(dihe) + " \n")

        tofile.close()
        return self

    '''def withSubGroup(self, isProtein=True, nucleic="nucleic-acid.lib"):

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

        return atominfor'''

    def arguements(self) :
        """Argument Parser

        Returns
        -------
        args: ArgParser object
        """

        d ='''
        ################################################################
        # Generate GMX Index from a PDB file                           #
        # Generate POSRES file for a PDB File                          #
        # Contact  Zheng Liangzhen, LZHENG002@E.NTU.EDU.SG             #
        # Version  2.3                                                 #
        # Update Dec 26, 2018                                          #
        ################################################################

        Usage examples:

        Generate backbone index for residue 1 to 1000 within Chain B
        gmx_index -f reference.pdb -res 1 100 -at backbone -chain B -o index.ndx

        Generate all atom index for residue 78 to 100
        gmx_index -f reference.pdb -o index.ndx -at allatom -chain ' ' -res 78 100

        Generate dihedral (PHI only) quadroplex index
        gmx_index -f input.pdb -o index.ndx -at dihedral -chain ' ' -res 78 100 -dihe PHI

        Note:
        The -s option is not used. So if you provide a -s reference.pdb, nothing will 
        be affected. 
        '''
        parser = gmxcli.GromacsCommanLine(d=d)

        parser.arguments()

        parser.parser.add_argument('-at', type=str, default ='allatom',
                                   help="Selected atom type for generating index. \n"
                                   "Options including: allatom, mainchain, \n"
                                   "non-hydrogen, c-alpha, backbone, sidechain, dihedral\n"
                                   "Default choice is: allatom \n")
        '''parser.parser.add_argument('-an', type=str, default=[], nargs='+',
                           help="Select atom by atom names. A list of atom names \n"
                           "could be supplied. If not given, all atom names are \n"
                           "considered. Default is [].")'''
        parser.parser.add_argument('-chain', type=str, default= "A", nargs="+",
                                   help="Protein chain identifier. Default chain ID is [\'A\',]. ")
        parser.parser.add_argument('-res', type=int, nargs= '+',
                                   help="Residue sequence number for index generating. \n"
                                   "Example, -res 1 100, generateing atom index within \n"
                                   "residues 1 to 100. Default is None.")
        parser.parser.add_argument('-posres', default=False,
                                   help="Generate a postion restraint file for selected index.\n"
                                   "Default name is posres.itp \n")
        parser.parser.add_argument('-dihe', default=None, type =str,
                                   help="Whether generate dihedral angles index (quadruplex).\n"
                                   "Phi and Psi are considered. Optional choices are: \n"
                                   "PHI, PSI, PHI_PSI, or NA. Default is NA. \n")
        parser.parser.add_argument('-gn', type=str, default=None,
                                   help="The name of the group of atoms selected. \n"
                                   "Default is None.")
        parser.parser.add_argument('-append', default=True, type=bool,
                                   help="Append the group of index to the output file. \n"
                                   "Options: True, False. \n"
                                   "Default is True.")

        parser.parse_arguments()
        args = parser.args

        # decide to print help message
        if len(sys.argv) < 3:
            # no enough arguments, exit now
            parser.parser.print_help()
            print("\nYou chose non of the arguement!\nDo nothing and exit now!\n")
            sys.exit(1)

        return args

    def genGMXIndex(self):
        """run geneIndex, initialize argument parser function.

        """
        args = self.arguements()

        self.resSeq = args.res
        self.res_seq()

        self.at = args.at
        self.dihedral = args.dihedral

        if os.path.exists(args.f):
            self.pdb = args.f
        elif os.path.exists(args.s):
            self.pdb = args.s
        else:
            print("Reference pdb file does not exist.")
            sys.exit(0)

        # generate index
        self.load_pdb()
        self.atomndx = self.res_index()

        if args.append:
            append = 'a'
        else:
            append = 'wb'
        #tofile = open(args.o, append)

        if args.dihedral[0] != "NA":
            write_dihedral = True
        else:
            write_dihedral = False

        self.atomList2File(self.atomndx, group_name=args.gn, write_dihe=write_dihedral,
                           out_filen=args.o, append=append)

        if args.posres:
            with open('posres.itp', 'w') as tofile:
                tofile.write("[ position_restraints ] \n; ai  funct  fcx    fcy    fcz  \n")
                for atom in self.atomndx:
                    tofile.write("%12d  1  1000  1000  1000  \n"
                                 % int(atom))

        print("\nGenerating Index File Completed!")
        return self


class GmxIndex(object):
    """Parse Gromacs Index File

    Parameters
    ------------
    index: str,
          the file name of the input index file

    Attributes
    ------------
    ndxlines: list,
          the lines in the index file
    groups: list,
          the list of groups in the index file
    totallines: int,
          total number of lines in the index file

    Methods
    -------
    groupsLineNumber()
        Get the line number for each of the groups
    groupContent(group)
        Fetch the content in a group given its group name

    """

    def __init__(self, index):
        if os.path.exists(index):
            self.index = index
        else:
            print("File {} not exists!".format(index))
            sys.exit(0)

        self.ndxlines   = open(self.index).readlines()
        self.groups     = [x.split()[1] for x in self.ndxlines if ("[" in x and "]" in x)]
        self.totallines = len(self.ndxlines)

    def groupsLineNumber(self):
        """Get the line number for each of the groups

        Parameters
        ----------

        Returns
        -------
        groupLN: orderedDict, dict
               the group-line_number information
        """

        groupLN = OrderedDict()

        linecount = 0
        for s in self.ndxlines:

            if "[" in s and "]" in s:
                groupLN[s.split()[1]] = linecount

            linecount += 1

        return groupLN

    def groupContent(self, group):
        """Fetch the content in a group given its group name

        Parameters
        ----------
        group: str,
              the name of a group

        Returns
        -------
        elements: list
            The element in each group. Each element in the list is
            in str format.

        """

        gln = self.groupsLineNumber()

        start_ln = gln[group] + 1

        if self.groups.index(group) == len(self.groups) - 1:
            end_ln = self.totallines - 1
        else:
            end_ln = gln[self.groups[self.groups.index(group) + 1]] - 1

        contents = [self.ndxlines[x].strip("\n")+" " for x in range(start_ln, end_ln+1)]
        contents = " ".join(contents)

        return contents.split()

    def writeNdxGroup(self, group, elements, output="output.ndx"):
        """ Write a list of elements for a group into a file

        Parameters
        ----------
        group : str,
               the name of a new group
        elements : list,
               the list containing elements of a group
        output : str, default = "output.ndx"
               the file name of an output index file

        Returns
        -------
        self : returns an instance of self.
        """

        PdbIndex().atomList2File(elements, group, output)
        return  self

def main():
    """Entry_point of the pdb index and gromacs index modules

    """

    ndx = PdbIndex()

    ndx.genGMXIndex()

    return None

