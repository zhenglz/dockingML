#!/usr/bin/env python
# -*- coding: utf-8 -*-

#####################################################
# Script for generating contact probability map     #
# Author: ZHENG Liangzhen                           #
# Email: LZHENG002@e.ntu.edu.sg                     #
# Version: V4.1                                     #
# Date: 23 Nov 2017                                 #
#####################################################

from mdanaly import gmxcli, pca
from dockml import pdbIO, index, algorithms

from matplotlib import pyplot as plt
from datetime import datetime

import math
import os
import sys
import numpy as np
import pandas as pd
import mdtraj as mt


class CoordinatesXYZ(object):
    """Extract atom coordinates using mdtraj

    Parameters
    ----------
    traj : mt.Trajectory
        A mdtraj trajectory object as input
    top : str, format pdb
        the reference pdb file name
    atom_selection : str, default is name CA
        the atom selection for pca calculation.

    Attributes
    ----------
    traj_ : mt.Trajectory object
        a mdtraj trajectory, where the coordinates are stored
    ref_ : mt.Trajectory
        The reference structure.
    topology_ : mt.Topology object
        the mdtraj trajectory topology object
    n_atoms_ : int,
        number of atoms in the trajectory
    superimposed_ : bool,
        whether the trajectory has been superimposed
    superpose_atom_indices_ : np.array,
        the atom indices, starting from 0, format is int
    xyz_ : np.ndarray, shape=[N, M]
        the original xyz coordinates of selected atoms
        N is number of samples, or frames
        M is the multiply of atoms * 3

    Methods
    -------
    superimpose()
        superimpose the trajectory to a reference pdb file
    xyz_coordinates(atom_indices)
        extract xyz coordinates given a list of atoms' indices

    See Also
    --------
    ContactMap

    """

    def __init__(self, traj, top, atom_selection="CA"):
        self.traj_ = traj
        if os.path.exists(top):
            self.ref_ = mt.load(top)
        else:
            print("Reference pdb structure %s is not existed or accessible." % top)
            sys.exit(0)

        # topology
        self.topology_ = self.ref_.topology

        self.n_atoms_ = self.traj_.n_atoms

        self.superimposed_ = False
        self.superpose_atom_indices_ = self.topology_.select(
            "name {}".format(atom_selection))

        self.xyz_ = None

    def superimpose(self):
        """
        Superimpose the trajectory to a reference structure.

        Returns
        -------
        self : return an instance of self.

        """
        if not self.superimposed_:
            self.traj_.superpose(self.ref_, frame=0,
                                 atom_indices=self.superpose_atom_indices_)
            self.superimposed_ = True

        return self

    def xyz_coordinates(self, atom_indices=None):
        """
        Extract xyz coordinates for selected atoms for a trajectory object

        Parameters
        ----------
        atom_indices : np.array
            the atom index for selected atoms

        Returns
        -------
        xyz : np.ndarray, shape = [N, M]
            N is number of samples, or frames
            M is the multiply of number of atoms and 3

        """

        if atom_indices is None:
            atom_indices = self.superpose_atom_indices_

        if not self.superimposed_:
            print("Trajectory is not superimposed. "
                  "Superimpose it to the reference now...")
            self.superimpose()

        # subset a trajectory
        traj = self.traj_.atom_slice(atom_indices=atom_indices, inplace=False)

        # extract the xyz coordinates
        xyz = traj.xyz.reshape((traj.xyz.shape[0], traj.xyz.shape[1]*3))
        self.xyz_ = xyz

        return xyz


class ContactMap(object):
    """Construct a contact map with mdtraj distance matrix

    The distance matrix for the contactmap is firstly constructed,
    and followed by a distance cutoff comparision, 1 is given if
    the distance in a bin of the matrix is less than the distance
    cutoff.

    Parameters
    ----------
    traj : mt.Trajectory object,
        the MDTraj trajectory object for distance calculation
    group_a : list, np.array, or pd.Seris
        the list of atom index in cmap x-axis
    group_b : list, np.array, or pd.Seris
        the list of atom index in cmap x-axis
    cutoff : float, default = 0.35
        the distance cutoff, default is 3.5 angstrom

    Attributes
    ----------
    atom_group_a : np.array
        the list of atom index for cmap x-axis
    atom_group_b : np.array
        the list of atom index for cmap y-axis
    atom_pairs_ : np.ndarray, shape=[N, 2]
        the atom pairs for distance calculation, N is the number
        of atom pairs.
    dist_matrix_ : np.ndarray, shape = [M, N]
        the distance matrix, M is the number of atoms for x-axis
        N is the number of atom pairs in y-axis
    cmap_ : np.ndarray, shape = [M, N]
        the contact map matrix, M is the number of atoms for x-axis
        N is the number of atoms in y-axis
    cmap_computed_ : bool
        whether the contact map has been calculated
    coord_number_ : int
        the number of contact number (coordination number)

    """

    def __init__(self, traj, group_a, group_b, cutoff=0.35):

        self.traj = traj
        self.cutoff = cutoff

        self.atom_group_a = group_a
        self.atom_group_b = group_b

        if not isinstance(self.atom_group_a, np.ndarray):
            self.atom_group_a = np.array(self.atom_group_a)
            self.atom_group_b = np.array(self.atom_group_b)

        self.atom_pairs_ = np.array([])
        self.generate_atom_pairs()

        self.dist_matrix_ = None
        self.cmap_ = None

        self.distmtx_computed_ = False
        self.cmap_computed_ = False

        self.coord_number_ = None

    def distance_matrix(self):
        """Calculate the atom distance matrix.

        Returns
        -------
        self: the instance itself

        """
        #if not self.distmtx_computed_:
        distmtx = mt.compute_distances(self.traj,
                                       self.atom_pairs_)

        self.dist_matrix_ = distmtx
        self.distmtx_computed_ = True

        return self

    def generate_cmap(self, shape='array', switch=False):
        """Calculate atom cmap data.

        Parameters
        ----------
        shape: str,
            the type of cmap output, array like or matrix like.
            Default is array. Options are: array, matrix
        switch: bool, default is False
            apply a switch function to the contact map calculation

        Returns
        -------
        self: the instance itself

        """
        if not self.distmtx_computed_:
            self.distance_matrix()

        if not self.cmap_computed_:
            if switch:
                sf = algorithms.BasicAlgorithm()
                vf = np.vectorize(sf.switchFuction)

                cmap = vf(self.dist_matrix_, d0=self.cutoff*2.0)
            else:
                cmap = (self.dist_matrix_ <= self.cutoff) * 1.0

            if shape == "array":
                pass
            elif shape == "matrix":
                cmap = cmap.reshape((cmap.shape[0], self.atom_group_a.shape[0],
                                     self.atom_group_b.shape[0]))
            else:
                pass

            self.cmap_ = cmap
            self.cmap_computed_ = True

        return self

    def generate_atom_pairs(self):
        """Generate atom pairs list.

        Returns
        -------
        self: the instance itself

        """

        if bool(self.atom_group_a.shape[0]) and bool(self.atom_group_b.shape[0]):
            list_a = np.repeat(self.atom_group_a, self.atom_group_b.shape[0])
            list_b = np.repeat(self.atom_group_b, self.atom_group_a.shape[0])
            list_b = list_b.reshape((self.atom_group_a.shape[0], -1)).T.ravel()

            self.atom_pairs_ = np.array(list(zip(list_a, list_b)))
        else:
            self.atom_pairs_ = np.array([])

        return self

    def coord_num(self):
        """
        Calculate the coordinate number using mdtraj

        Returns
        -------
        self: the instance itself

        """
        if not self.cmap_computed_:
            self.generate_cmap(shape="array", switch=False)

        self.coord_number_ = np.sum(self.cmap_, axis=1)

        return self


class CmapNbyN(ContactMap):
    """Generate side-chain contact map

    This class is inherited from ContactMap

    Parameters
    ----------
    traj
    resids_a
    resids_b
    cutoff

    Attributes
    ----------
    traj : mt.Trajectory
        The input trajectory for analysis
    atom_group_a : np.array
        the list of atom index for cmap x-axis
    atom_group_b : np.array
        the list of atom index for cmap y-axis
    atom_pairs_ : np.ndarray, shape=[N, 2]
        the atom pairs for distance calculation, N is the number
        of atom pairs.
    dist_matrix_ : np.ndarray, shape = [M, N]
        the distance matrix, M is the number of atoms for x-axis
        N is the number of atom pairs in y-axis
    cmap_ : np.ndarray, shape = [M, N]
        the contact map matrix, M is the number of atoms for x-axis
        N is the number of atoms in y-axis
    cmap_computed_ : bool
        whether the contact map has been calculated
    coord_num_
    resids_a
    resids_b
    top_
    contacts_by_res_

    See Also
    --------
    ContactMap

    TODO: Bugs to be fixed in this module

    """

    def __init__(self, traj, resids_a, resids_b, cutoff=0.35):

        super().__init__(traj=traj, group_a=[], group_b=[], cutoff=cutoff)

        self.resids_a_ = resids_a
        self.resids_b_ = resids_b

        self.top_ = self.traj.topology
        self.contacts_by_res_ = np.array([])

    def single_pair_dist(self, resid_a, resid_b, atomtype=["all"]):
        """

        Parameters
        ----------
        resid_a : int
            The residue id for x-axis
        resid_b : int
            The residue id for y-axis

        Returns
        -------

        """

        atoms_list_a = self.top_.select("resid %d and %s"
                                        % (resid_a, atomtype[0]))
        atoms_list_b = self.top_.select("resid %d and %s"
                                        % (resid_b, atomtype[-1]))

        self.atom_group_a = atoms_list_a
        self.atom_group_b = atoms_list_b

        #print(self.atom_group_a, self.atom_group_b)

        self.generate_atom_pairs()
        #self.distance_matrix()
        self.generate_cmap()
        self.coord_num()

        nbyn = self.coord_number_ / self.atom_pairs_.shape[0]

        return nbyn.reshape((-1, 1)), self.coord_number_.reshape((-1, 1))

    def cmap_nbyn(self, atomtype=['all', 'all']):
        """The main function for residue-residue sidechain contact
        normalized by atom number products.

        Returns
        -------
        self : an instance of itself

        """

        cmap = np.array([])

        for res_a in self.resids_a_:
            for res_b in self.resids_b_:
                nbyn = self.single_pair_dist(res_a, res_b, atomtype)[0]
                if cmap.shape[0] == 0:
                    cmap = nbyn
                else:
                    cmap = np.concatenate((cmap, nbyn), axis=1)

        self.cmap_ = cmap

        return self

    def contact_nbyn(self, atomtype=['heavy', 'heavy']):
        contacts = np.array([])
        for res_a in self.resids_a_:
            for res_b in self.resids_b_:
                #print(res_a, res_b)
                nbyn = self.single_pair_dist(res_a, res_b)[1]
                #print(nbyn)
                if contacts.shape[0] == 0:
                    contacts = nbyn
                else:
                    contacts = np.concatenate((contacts, nbyn), axis=1)

        self.contacts_by_res_ = contacts

        return self


class DrawCMap(object):

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
        if len(refpdb) == 3:

            ppdb = pdbIO.parsePDB("")
            fullreslist = ppdb.getNdxForRes(refpdb[0], [refpdb[1]])

            shortresmap = ppdb.longRes2ShortRes()
            fullreslist = [x for x in fullreslist if x[2] in refpdb[1]]

            for resk in key_res:
                if fullreslist[resk][0] in shortresmap.keys():
                    resseq = str(resk + refpdb[2])
                    resname= shortresmap[fullreslist[resk][0]]
                    chainid= fullreslist[resk][2]

                    if yticks_showchainid:
                        id = resname + resseq + chainid
                    else :
                        id = resname + resseq

                    res_labels.append(id)

        # only keep the important residue cmap
        keyres_cmap = np.asarray(cm_sorted)[:, 1:][:, list(key_res)]

        # get the length of x and y axis
        shapex = len(key_res)
        shapey = cmapdata.shape[0] + 1

        print("Protein residue numbers: %d" % shapex)
        print("Time point numbers: %d"% shapey)

        z = np.transpose(keyres_cmap[:, 1:]).T

        plt.pcolormesh(z.T, cmap=plt.get_cmap(cmaptype))

        if colorbar_show :
            plt.colorbar(label=colorbar_label)

        plt.xlabel(xlabel, fontsize=fsize)
        plt.ylabel(ylabel, fontsize=fsize)

        if len(yticks_loc) and len(yticks_labels):
            plt.yticks(yticks_loc, yticks_labels)
        else:
            if len(refpdb) == 3:
                plt.yticks(np.array(range(shapex))+0.5, res_labels)

        if len(xlim):
            plt.xlim(xlim)

        if len(ylim):
            plt.ylim(ylim)
        else:
            plt.ylim([0, shapex+0.5])

        if len(xticks_labels) and len(xticks_loc):
            plt.xticks(xticks_loc, xticks_labels)

        if len(savefig):
            plt.savefig(savefig, dpi=2000)

        plt.show()

        return 1


class CommunityCmap(object):
    """Calculate and parse side-chain contact map for community
    network analysis.

    Parameters
    ----------
    cmap : np.ndarray, or a cmap file name
        str: the filename of a contact map
        np.ndarray: a contactmap matrix np.ndarray object
    sep : str, default = ','
        the delimiter used in the cmap file.

    Attributes
    ----------
    cmap_ : np.ndarray

    """

    def __init__(self, cmap, sep=",", is_file=False):

        if is_file and os.path.exists(cmap):
            try:
                self.cmap_ = np.loadtxt(cmap,
                                        delimiter=sep,
                                        comments="#")
            except FileExistsError:
                print("File %s not exists! " % cmap)
        else:
            self.cmap_ = cmap

        self.final_map_ = None

    def icriticalMap(self, icritical, dat, cutoff=0.48):
        """Calculate whether the icritical map

        Parameters
        ----------
        icritical : float
            The
        dat: np.ndarray, shape = [N, M]
            The contact map matrix, N is number of residues

        Returns
        -------

        """

        self.final_map_ = np.sum(np.greater(dat, icritical), axis=0) / dat.shape[0]
        self.final_map_ = (self.final_map_ >= cutoff) * 1.0

        return self.final_map_

    def diagonalZeroed(self, cmap, xsize, outfile):
        from mdanaly import matrix

        map = cmap.reshape((xsize, cmap.shape[0] / xsize))

        # matrix 2 xyz
        xyz = matrix.MatrixHandle().matrix2xyz(map)
        mtx = matrix.MatrixHandle().neiborhood2zero(xyz, neiborsize=4, outtype='mtx')

        np.savetxt(outfile, mtx, fmt="%3.1f", delimiter=" ")

        return mtx

    def generateCmap(self):

        Icritical = np.arange(0, 15, 1.0)

        for ic in Icritical:
            imap = self.icriticalMap(ic, self.cmap_)

            self.diagonalZeroed(imap, int(np.sqrt(self.cmap_.shape[1])), "cmap_I_%.2f" % ic)

        return 1

    def calculateNbyN(self, pdbfile, dcutoff, res1, res2, ndxlist):

        """
        Todo: To be implemented.

        Parameters
        ----------
        pdbfile
        dcutoff
        res1
        res2
        ndxlist

        Returns
        -------

        """

        return NotImplementedError

    def scatterFileList(self, ranksize, pdbFileList):
        load4each = int(math.ceil(float(len(pdbFileList)) / float(ranksize)))
        filesList = []

        for i in range(ranksize - 1):
            filesList.append(pdbFileList[i * load4each: load4each * (i + 1)])
        filesList.append(pdbFileList[(ranksize - 1) * load4each:])

        return filesList


def descriptions():
    """
    Generate utility description.

    Returns
    -------
    d: str,
        the description content
    """

    d = '''
        ########################################################################
        #  Generating contact probability map                                  #
        #  Author:  ZHENG Liangzhen & Mu Yuguang                               #
        #  Email:   LZHENG002@e.ntu.edu.sg                                     #
        #  Version: V2.2                                                       #
        #  Date:    27 Dec 2017                                                #
        ########################################################################

        Generating contact probability Map (Cmap)

        Input a gromacs trajectory file (.xtc) to construct a contact probability map.
        All the structures should stay whole, broken structures will cause inaccurate results.
        All the frames in trajectory file do not consider PBC conditions, you should keep structures
        as a whole.

        If some arguements not given, default values would be used.

        Usage:
        1. Show help information
        gmx_cmap.py -h

        2. Construct a Ca-Ca Cmap for a protein chain
        gmx_cmap -f traj.xtc -out Cmap.dat -rc A 1 250
        -lc A 1 250 -cutoff 0.5 -switch T -atomtype CA

        3. Generate a full-atom Cmap for a poly-peptide chain
        gmx_cmap -f traj.xtc -o Cmap.dat -rc A 1 250
        -lc A 1 250 -cutoff 0.5 -atomtype all all

        4. Construct a Cmap between a small ligand and a protein
        gmx_cmap -f traj.xtc -o Cmap.dat -rc A 1 250
        -lc A 251 251 -cutoff 0.5 -atomtype all all

        5. Construct a Cmap between a small ligand and a protein, Ca-allatom
        gmx_cmap -f traj.xtc -out Cmap.dat -rc A 1 250
        -lc A 251 251 -cutoff 0.5 -atomtype CA all

        '''
    return d


def arguments():

    parser = gmxcli.GromacsCommanLine(d=descriptions())

    parser.arguments()

    parser.parser.add_argument('-rc', type=str, nargs='+', default=['A', '1', '250'],
                               help="Input, optional. \n"
                                    "The receptor chains and residue index for Cmap construction.\n"
                                    "You must enter a chain name, start residue index, and end chain index.\n"
                                    "Default is: A 1 250 \n")
    parser.parser.add_argument('-lc', type=str, nargs='+', default=['A','1','250'],
                               help="Input, optional. \n"
                                    "The ligand chains and residue index for Cmap construction.\n"
                                    "You must enter a chain name, start residue index, and end chain index.\n"
                                    "Default is: A 1 250 \n")
    parser.parser.add_argument('-cutoff',type=float, default=0.35,
                               help="Input, optional. Default is 0.35 (nanometer). \n"
                                    "Distance Cutoff for determining contacts. \n")
    parser.parser.add_argument('-atomtype', type=str, nargs='+', default=[],
                               help="Input, optional. \n"
                                    "Atom types for Receptor and Ligand in Contact Map Calculation. \n"
                                    "Only selected atoms will be considered.\n"
                                    "Options: CA, Backbone, MainChain, All, non-H(All-H), lig-all. \n"
                                    "CA, alpha-carbon atoms. Backbone, backbone atoms in peptides. \n"
                                    "MainChain, including CA and N atoms. All, means all atoms.\n"
                                    "non-H, non-hydrogen atoms, all the heavy atoms. \n"
                                    "lig-all trys to consider all the atoms of a ligand (H atoms not considered). \n"
                                    "Two choices should be provided for receptor and ligand respectively. \n"
                                    "If only one atomtype given, the 2nd will be the same as 1st.\n"
                                    "Default is: [] \n")
    parser.parser.add_argument('-atomname1', type=str, nargs='+', default=[],
                               help="Input, optional. \n"
                                    "Atom names for Recetpor in Contact Map. \n"
                                    "Default is []. ")
    parser.parser.add_argument('-atomname2', type=str, nargs='+', default=[],
                               help="Input, optional. \n"
                                    "Atom names for Ligand in Contact Map. \n"
                                    "Default is []. ")
    parser.parser.add_argument('-eletype', type=str, nargs="+", default=[],
                               help="Input, optional. \n"
                                    "Choose the specific elements for atom indexing to construct the cmap.\n"
                                    "Default is [].\n")
    parser.parser.add_argument('-switch', type=lambda x: (str(x).lower() == "true"), default=False,
                               help="Input, optional. Default is False. \n"
                                    "Apply a switch function for determing Ca-Ca contacts for a smooth transition. \n"
                                    "Only work with atomtype as CA. \n")
    parser.parser.add_argument('-NbyN', type=lambda x: (str(x).lower() == "true"), default=False,
                               help="Input, optional. Default is False\n"
                                    "For community analysis, calculate atom contact number, normalized. \n")
    parser.parser.add_argument('-details', default=None, type=str,
                               help="Provide detail contact information and write out to a file. \n"
                                    "Default is None.")
    parser.parser.add_argument('-opt', default="TS", type=str,
                               help="Optional setting controls. Default is A. \n"
                                    "Average, calculating the average cmap along the simulations.\n"
                                    "Separated, create time series contact map, suitable for ligand \n"
                                    "protein contact information along time.\n"
                                    "Options: S(Separated), A(Average).\n")

    parser.parse_arguments()

    return parser.args


def verbose(verbose=True, s=""):
    if verbose:
        print(s)


def cmap_general(trajs, ref, rc, lc, at, cutoff=0.35, v=True, switch=False):
    """Computate general type of contact maps

    Parameters
    ----------
    trajs : list of mt.trajectory objects, shape = N
        The input chunks of trajectories, N is number of chunks
    ref : str
        The reference pdb file name
    rc : list, shape = 3
        The residue and chain identifier for x-axis
    lc : list, shape = 3
        The residue and chain identifier for y-axis
    at : str
        The atomtype of atoms used for cmap generation
    cutoff : float, default = 0.35
        The distance cutoff, in unit nanometer
    v : bool, default = True
        Whether print detail information during the calculation
    switch : bool, default = False
        Whether apply a switch function to have continuous contact map

    Returns
    -------
    contact_map : np.ndarray, shape = [ N, M ]
        The output contact map dataset
        N is number of frames, M is number of dimensions.
    """

    contact_map = np.array([])

    # receptor (x-axis) atom selection
    group_a = index.gen_atom_index(pdbin=ref, chain=[rc[0], ], resSeq=rc[1:],
                                   atomtype=at[0], style="mdtraj")

    # ligand (y-axis) atom selection
    group_b = index.gen_atom_index(pdbin=ref, chain=[lc[0], ], resSeq=lc[1:],
                                   atomtype=at[-1], style="mdtraj")

    verbose(verbose=v, s="Atom indices have been processed ......")

    for i, traj in enumerate(trajs):
        # calculate cmap information
        verbose(v, "Generate cmap for chunk %5d ......" % i)
        contmap = ContactMap(traj, group_a, group_b, cutoff=cutoff)
        contmap.generate_cmap(shape="array", switch=switch)

        if i == 0:
            contact_map = contmap.cmap_
        else:
            contact_map = np.concatenate((contact_map, contmap.cmap_),
                                         axis=0)

    return contact_map


def cmap_nbyn(trajs, ref, rc, lc, v=True,
              cutoff=0.35, allchains=" ABCDEFGH",
              atomtype=["sidechain", "sidechain"]):
    """Generate sidechain based contact maps.

    Parameters
    ----------
    trajs : list of mt.trajectory objects, shape = N
        The input chunks of trajectories, N is number of chunks
    ref : str
        The reference pdb file name
    rc : list, shape = 3
        The residue and chain identifier for x-axis
    lc : list, shape = 3
        The residue and chain identifier for y-axis
    cutoff : float, default = 0.35
        The distance cutoff, in unit nanometer
    v : bool, default = True
        Whether print detail information during the calculation
    allchains : str, default = 'ABCDEFGH'
        All available chain identifiers in the reference pdb files
    Returns
    -------

    """

    pdb = pdbIO.parsePDB(inPDB=ref)
    all_resids = pdb.getNdxForRes(ref, chains=allchains)

    print(all_resids)

    # for chain_a
    resids_a = []
    for i, item in enumerate(all_resids):
        if item[2] in rc[0] and \
                int(item[1]) in np.arange(int(rc[1]), int(rc[2])+1):
            resids_a.append(i)
    print(resids_a)

    # for chain_b
    resids_b = []
    for i, item in enumerate(all_resids):
        if item[2] in lc[0] and \
                int(item[1]) in np.arange(int(lc[1]), int(lc[2])+1):
            resids_b.append(i)
    print(resids_b)

    contact_map = np.array([])
    for i, traj in enumerate(trajs):
        # calculate cmap information
        verbose(v, "Generate cmap for chunk %5d ......" % i)
        contmap = CmapNbyN(traj, resids_a=resids_a,
                           resids_b=resids_b, cutoff=cutoff)

        contmap.cmap_nbyn(atomtype=atomtype)
        if i == 0:
            contact_map = contmap.cmap_
        else:
            contact_map = np.concatenate((contact_map, contmap.cmap_),
                                         axis=0)

    return contact_map


def iterload_cmap():
    """Load large trajectory iteratively using mdtraj.iterload function,
    then generate contact map from the trajectories.

    Returns
    -------

    """

    # for calculation time counting
    startTime = datetime.now()

    # argument options
    args = arguments()

    verbose(args.v, "Atom selecting ......")
    # TODO: atom selection method required
    if os.path.exists(args.s) and args.s[-4:] == ".pdb":
        inp = args.s
    elif os.path.exists(args.f) and args.f[-4:] == ".pdb":
        inp = args.f
    else:
        inp = None
        print("Reference pdb file is not existed. Exit now!")
        sys.exit(0)

    rec_index = int(args.rc[2]) - int(args.rc[1]) + 1
    lig_index = int(args.lc[2]) - int(args.lc[1]) + 1

    verbose(args.v, "Loading trajectory ......")
    # read gromacs trajectory
    trajs = gmxcli.read_xtc(args.f, args.s, chunk=1000,
                            stride=int(args.dt/args.ps))

    n_frames = sum([x.n_frames for x in trajs])
    verbose(args.v, "Total number of frames: %d " % n_frames)

    verbose(args.v, "Start calculating contact map ......")
    if args.NbyN:
        contact_map = cmap_nbyn(trajs, inp, args.rc, args.lc,
                                args.v, args.cutoff,
                                " ABCDEFGHIJK", args.atomtype)
    else:
        contact_map = cmap_general(trajs, inp, args.rc, args.lc,
                                   args.atomtype, args.cutoff,
                                   v=args.v, switch=args.switch)

    # subset the results
    verbose(args.v, "Preparing output file ......")
    contact_map = pd.DataFrame(contact_map)
    contact_map.index = np.arange(contact_map.shape[0]) * args.dt
    contact_map = pca.datset_subset(contact_map, args.b, args.e)

    # get mean cmap data
    if args.opt in ['A', 'a', 'average', 'Average']:
        results = np.mean(contact_map, axis=0).reshape((rec_index, lig_index))
        results = pd.DataFrame(results)
        results.index = range(int(args.rc[1]), int(args.rc[2]) + 1)
        results.columns = range(int(args.lc[1]), int(args.lc[2]) + 1)
    else:
        results = contact_map
        results = pd.DataFrame(results)
        results.index = np.arange(results.shape[0]) * args.dt
        results.columns = [str(x) for x in np.arange(results.shape[1])]

    # save results to an output file
    verbose(args.v, "Writing output now ...... ")
    results.to_csv(args.o, sep=",", header=True, index=True, float_format="%.3f")

    print("Total Time Usage: ")
    print(datetime.now() - startTime)

