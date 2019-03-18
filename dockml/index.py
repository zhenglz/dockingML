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
import numpy as np


class PdbIndex(object):
    """Input the residue number sequence, then out put the required index atoms

    Parameters
    ----------
    dihedral : list, default = ["NA"]
        options = ["NA", "PSI", "PHI"]
        The type of dihedrals. 'NA' indicates none dihedral angles will be
        selected.
    reference : str, default = 'input.pdb'
        the input pdb file name
    chain : list, default = [ 'A', ]
        the chain identifiers in the pdb file name
    atomtype : list, default = [ 'CA', ]
        the atom names for atoms

    Attributes
    ----------
    phi: list,
        residue phi torsion angle atom name list
    psi: list
        residue psi torsion angle atom name list
    at_direct : list,
        the molecule types for atom selections
    at : str
        options: CA, ca, mainchain, sidechain, water, protein,
        water, all, none
        the atom type for selection.
    atomndx : list
        the atom index of selected molecules and residues
    chain_ndx : list
        The chain index list, starting from 0.
    selections : str
        The DSL selection quary for mt.topology.select

    Examples
    --------
    >>> from dockml import index
    Backend Qt5Agg is interactive backend. Turning interactive mode on.
    >>> ndx = index.PdbIndex("./test/input.pdb", atomtype="all", chain=["C"])
    >>> # prepare necessary parameters for index generation
    >>> ndx.prepare_selection()
    Out[4]: <dockml.index.PdbIndex at 0x7f7cb917a668>
    >>> ndx.selections_
    Out[5]: ['']
    >>> ndx.res_index()
    Out[6]: <dockml.index.PdbIndex at 0x7f7cb917a668>
    >>> ndx.selections_
    Out[7]: 'name  and resid 1 to 40 and chainid 0'
    >>> len(ndx.atomndx_)
    Out[8]: 602
    >>> ndx.chain_ndx
    Out[9]: [0]
    >>> ndx.chain
    Out[10]: ['C']
    >>> ndx.at
    Out[11]: 'all'
    >>> ndx.atomndx_mt_style_[:10]
    Out[12]: [22407, 22408, 22409, 22410, 22411, 22412, 22413, 22414, 22415, 22416]

    """

    def __init__(self, reference="input.pdb",
                 chain=["A"], resSeq=[1, -1],
                 atomtype="CA", dihedral=["NA"]):

        self.phi = ['C', 'N', 'CA', 'C']
        self.psi = ['N', 'CA', 'C', 'N']
        self.at_direct = ["all", "none", "sidechain", "mainchain",
                          "protein", "water", ]

        # diheral is a list
        self.dihedral = dihedral
        self.pdb = reference
        self.top = None
        self.pdb_lines_ = []

        # chain identifiers, list
        self.chain = chain
        self.resSeq = [int(x) for x in resSeq]
        self.at = atomtype
        self.atomndx_ = np.array([])
        self.atomndx_mt_style_ = np.array([])
        #self.molecule = None

        self.resid_mapper_ = {}
        self.resid_mapped_ = False

        self.selections_ = [""]
        self.start_atom = 1
        self.chain_ndx = [0, ]

        self.index_calculated_ = False

    def load_pdb(self):
        """Load pdb file for mdtraj processing

        Returns
        -------
        self : return an instance of self.
        """

        if not os.path.exists(self.pdb):
            print("PDB reference file %s is not existed, exit now!" % self.pdb)
            sys.exit(0)

        if self.pdb[-4:] != ".pdb":
            print("You just loaded a file %s \nPlease provide a reference pdb file."
                  %self.pdb)

        traj = mt.load(self.pdb)
        self.top = traj.topology

        with open(self.pdb) as lines:
            lines = [int(x.split()[1]) for x in lines if
                     len(x.split()) and (x.split()[0] in ["ATOM", "HETATM"])]
            self.start_atom = lines[0]

        return self

    def res_seq(self):
        """Residue ID index processing.

        Returns
        -------
        self : return an instance of self.
        """
        if len(self.resSeq) == 1:
            self.resSeq = [self.resSeq[0], self.resSeq[0]]

        if self.resSeq[1] == -1:
            self.resSeq[1] = self.top.n_residues - self.resSeq[0] + 1

        return self

    def resid_mt_style(self, chain, start_res, end_res):
        """Query the 0-based resid given the chain and residue
        sequence information.

        Parameters
        ----------
        chain
        start_res
        end_res

        Returns
        -------
        resid_start : int
        resid_end : int
        """
        if not self.resid_mapped_:
            self.resid_mapper()

        start = chain + "_" + str(start_res)
        end = chain + "_" + str(end_res)

        return self.resid_mapper_[start], self.resid_mapper_[end]

    def resid_mapper(self):

        chain_seq = []

        with open(self.pdb) as lines:
            lines = [x for x in lines if x.split()[0] in ["ATOM", "HETATM"]]
            self.pdb_lines_ = lines + []
            for s in lines:
                if s[21] + "_" + s[22:26].strip() not in chain_seq:
                    chain_seq.append(s[21] + "_" + s[22:26].strip())

        resid_dict = {}
        for i, id in enumerate(chain_seq):
            resid_dict[id] = i

        self.resid_mapper_ = resid_dict
        self.resid_mapped_ = True

        return self

    def prepare_selection(self):
        """Prepare necessary information for Gromacs index generation.

        Returns
        -------
        self: return a instance of itself
        """

        self.load_pdb()
        self.res_seq()
        #self.get_chains()
        self.resid_mapper()

        return self

    def res_index(self, atom_name_list=[]):
        """Get atom index list from given information:
        Resid, chain, atom_type_names

        Parameters
        ----------
        atom_name_list : list, default is empty.
            The list of atom names provided for index selection.

        Returns
        -------
        self : return an instance of self.
        """

        if self.dihedral[0] != "NA":
            indexlist = []
            for resndx in range(self.resSeq[0], self.resSeq[-1]+1):
                phitype = self.phi
                phipair = [-1, -1, -1, -1]
                phindxadd = [-1, 0, 0, 0]
                for i in range(4):
                    with open(self.pdb) as lines:
                        for s in lines:
                            if len(s.split()) > 1 and s.split()[0] == "ATOM" \
                                    and s[21] in self.chain \
                                    and int(s[22:26].strip()) == resndx + phindxadd[i] \
                                    and s[12:16].strip() == phitype[i]:
                                phipair[i] = int(s.split()[1])
                #indexlist.append(phipair)

                psitype = self.psi
                psindxadd = [0, 0, 0, 1]
                psipair = [-1, -1, -1, -1]
                for i in range(4):
                    with open(self.pdb) as lines:
                        for s in lines:
                            if len(s.split()) > 1 and s.split()[0] == "ATOM" \
                                    and s[21] in self.chain and \
                                    int(s[22:26].strip()) == resndx + psindxadd[i] \
                                    and s[12:16].strip() == psitype[i]:
                                psipair[i] = int(s.split()[1])
                                # break
                #indexlist.append(psipair)

                # supposing atom index larger than 0
                if "PSI" in self.dihedral[0] and all(np.array(phipair)>0) > 0:
                    indexlist.append(psipair)
                if "PHI" in self.dihedral[0] and all(np.array(psipair)>0) > 0:
                    indexlist.append(phipair)

            self.atomndx_ = indexlist
            print(self.atomndx_)
        else:

            if len(atom_name_list) == 0:
                if self.at in self.at_direct:
                    names = self.at
                elif self.at in ['Heavy', 'heavy', 'non-hdyrogen', 'no-H', 'non-H']:
                    names = "symbol != 'H'"
                else:
                    names = "name %s"%self.at
            else:
                names = "name " + "  ".join(atom_name_list)

            start, end = self.resid_mt_style(self.chain[0],
                                             self.resSeq[0],
                                             self.resSeq[1])
            resndx = "resid %d to %d" % (start, end)

            self.selections_ = names + " and " + resndx
            print(self.selections_)
            try:
                self.atomndx_mt_style_ = self.top.select(self.selections_)
            except:
                print("Error: Atom indexing not successful. ")
                #self.atomndx_mt_style_ = np.array([])
            #self.atom_index_original()
            #self.atomndx = [x + self.start_atom for x in self.atomndx_mt_style]

            self.index_calculated_ = True

        return self

    def atom_index_original(self, dihedral=False):

        if not dihedral:
            if not self.resid_mapped_:
                self.resid_mapper()

            if not self.atomndx_mt_style_.shape[0]:
                if not self.index_calculated_:
                    self.res_index()



            try:
                lines_selected = np.array(self.pdb_lines_)[self.atomndx_mt_style_]
                self.atomndx_ = [int(x.split()[1]) for x in lines_selected]
            except:
                print("Check whether you have generated the mt_style index first. ")
                self.atomndx_ = np.array([])

        return self

    def atomList2File(self, atom_list, group_name, write_dihe=False,
                      out_filen="output.ndx", append=True):
        """Save a group of atom index into a gromacs-format index file

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
        self : return an instance of self.
        """
        if append:
            tofile = open(out_filen, 'a')
        else:
            tofile = open(out_filen, 'wb')

        tofile.write('[ %s ] \n'% group_name)

        if not write_dihe:
            for i, atom in enumerate(atom_list):

                if i % 15 == 0 and i != 0:
                    tofile.write('  \n')
                tofile.write('%6d ' % int(atom))
            tofile.write("\n")

        else:
            for dihe in atom_list:
                tofile.write("   ".join([ str(x) for x in dihe]) + " \n")

        tofile.close()
        return self

    def arguements(self):
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
        parser.parser.add_argument('-an', type=str, default=[], nargs='+',
                                   help="Select atom by atom names. A list of atom names \n"
                                   "could be supplied. If not given, all atom names are \n"
                                   "considered. Default is [].")
        parser.parser.add_argument('-chain', type=str, default=["A"], nargs="+",
                                   help="Protein chain identifier. Default chain ID is [\'A\',]. ")
        parser.parser.add_argument('-res', type=int, nargs= '+', default=[1, -1],
                                   help="Residue sequence number for index generating. \n"
                                   "Example, -res 1 100, generateing atom index within \n"
                                   "residues 1 to 100. Default is None.")
        parser.parser.add_argument('-posres', default=False,
                                   help="Generate a postion restraint file for selected index.\n"
                                   "Default name is posres.itp \n")
        parser.parser.add_argument('-dihe', default=["NA"], type =str, nargs="+",
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

        Returns
        -------
        self : return an instance of self.

        """
        args = self.arguements()

        # define residue index
        self.resSeq = args.res

        # define chains
        self.chain = args.chain

        # define molecule name or atom name
        self.at = args.at
        # define dihedral name
        self.dihedral = args.dihe

        if os.path.exists(args.f):
            self.pdb = args.f
        elif os.path.exists(args.s):
            self.pdb = args.s
        else:
            print("Reference pdb file does not exist.")
            sys.exit(0)

        # prepare
        self.prepare_selection()
        self.res_index(atom_name_list=args.an)
        self.atom_index_original(dihedral=args.dihe)
        print(self.atomndx_)
        self.atomList2File(self.atomndx_, group_name=args.gn,
                           write_dihe=(self.dihedral[0] != "NA"),
                           out_filen=args.o, append=args.append)

        if args.posres:
            with open('posres.itp', 'w') as tofile:
                tofile.write("[ position_restraints ] \n"
                             "; ai  funct  fcx    fcy    fcz  \n")
                for atom in self.atomndx_:
                    tofile.write("%12d  1  1000  1000  1000  \n"
                                 % int(atom))

        print("\nGenerating Index File Completed!")
        return self


class GmxIndex(object):
    """Parse Gromacs Index File

    This class is a parser for gromacs index file. The group names and
    their elements were parsed.

    Parameters
    ------------
    index : str
          the file name of the input index file

    Attributes
    ------------
    ndxlines : list
          the lines in the index file
    groups : list
          the list of groups in the index file
    totallines : int
          total number of lines in the index file

    Methods
    -------
    groupsLineNumber()
        Get the line number for each of the groups
    groupContent(group)
        Fetch the content in a group given its group name

    Examples
    --------
    >>> # process index
    >>> ndx = index.GmxIndex("index.ndx")
    >>> sets = ["receptor", "ligand"]
    >>> used_groups = []
    >>> atom_indices = []
    >>> for i in [0, 1]:
    ...     print("Please select a group for %s: " % sets[i])
    ...     for j, gn in enumerate(ndx.groups):
    ...         print("%d : %s" % (j, gn))
    ...
    ...     used_groups.append(ndx.groups[int(input("Your choice: "))])
    >>> rec_ndx = [int(x)-1 for x in ndx.groupContent(used_groups[0])]
    >>> lig_ndx = [int(x)-1 for x in ndx.groupContent(used_groups[1])]

    """

    def __init__(self, index):
        if os.path.exists(index):
            self.index = index
        else:
            print("File {} not exists!".format(index))
            sys.exit(0)

        self.ndxlines = open(self.index).readlines()
        self.groups = [x.split()[1] for x in self.ndxlines
                       if ("[" in x and "]" in x)]
        self.totallines = len(self.ndxlines)

    def groupsLineNumber(self):
        """Get the line number for each of the groups

        Parameters
        ----------

        Returns
        -------
        groupLN: collections.orderedDict, dict
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
        ndx = PdbIndex()
        ndx.atomList2File(atom_list=elements, group_name=group, out_filen=output)
        return self


def gen_atom_index(pdbin, chain, resSeq, atomtype=['CA'], style='mdtraj'):

    ndx = PdbIndex(reference=pdbin, chain=chain, resSeq=resSeq, atomtype=atomtype, )
    ndx.prepare_selection()
    ndx.res_index()

    #print(ndx.selections_)

    if style == "mdtraj":
        return ndx.atomndx_mt_style_
    elif style == "original":
        ndx.atom_index_original()
        return ndx.atomndx_
    else:
        print("Warning: The atom indices style should be either mdtraj or original. Use mdtraj here. ")
        return ndx.atomndx_mt_style_


def main():
    """Entry_point of the pdb index and gromacs index modules
    """

    ndx = PdbIndex()

    ndx.genGMXIndex()

    return None

