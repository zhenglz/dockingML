#!/usr/bin/env python
# -*- coding: utf-8 -*-

#####################################################
# Script for generating contact probability map     #
# Author: ZHENG Liangzhen                           #
# Email: LZHENG002@e.ntu.edu.sg                     #
# Version: V4.1                                     #
# Date: 23 Nov 2017                                 #
#####################################################

import glob, math, sys, os
import numpy as np
import pandas as pd
import mdanaly

from mdanaly import gmxcli
from mdanaly import pca
import dockml.pdbIO as pio
from dockml import index as ndx

from matplotlib import pyplot as plt
from collections import defaultdict
from datetime import datetime
from mpi4py import MPI
import mdtraj as mt


class CoordinatesXYZ(object):
    """
    Perform trajectory coordinates PCA analysis using mdtraj

    Parameters
    ----------
    traj: mdtraj.Trajectory object,
        a mdtraj trajectory, where the coordinates are stored
    top: str, format pdb
        the reference pdb file name
    atom_selection: str, default is name CA
        the atom selection for pca calculation.

    Attributes
    ----------
    topology: mdtraj.Trajectory.topology object
        the mdtraj trajectory topology object
    n_atoms_: int,
        number of atoms in the trajectory
    superimposed_: bool,
        whether the trajectory has been superimposed
    superpose_atom_indices: numpy array,
        the atom indices, starting from 0, format is int
    xyz: numpy ndarray, shape=[N, M]
        the original xyz coordinates of selected atoms
        N is number of samples, or frames
        M is the multiply of atoms * 3

    Methods
    -------
    superimpose
    xyz_coordinates

    See Also
    --------
    ContactMap

    """

    def __init__(self, traj, top, atom_selection="CA"):
        self.traj = traj
        self.ref = mt.load(top)

        # topology
        self.topology = self.ref.topology

        self.n_atoms_ = self.traj.n_atoms

        self.superimposed_ = False
        self.superpose_atom_indices_ = self.topology.select("name {}".format(atom_selection))

        self.xyz = None

    def superimpose(self):
        """
        Superimpose the trajectory to a reference structure.

        Returns
        -------
        self: the object itself

        """

        self.traj.superpose(self.ref, frame=0, atom_indices=self.superpose_atom_indices_)
        self.superimposed_ = True

        return self

    def xyz_coordinates(self, atom_indices=None):
        """
        Extract xyz coordinates for selected atoms for a trajectory object

        Parameters
        ----------
        atom_indices: numpy array,
            the atom index for selected atoms

        Returns
        -------
        xyz: numpy ndarray, shape = [N, M]
            N is number of samples, or frames
            M is the multiply of number of atoms and 3

        """

        if atom_indices is None:
            atom_indices = self.superpose_atom_indices_

        if not self.superimposed_:
            print("Trajectory is not superimposed. Superimpose it to reference now...")
            self.superimpose()

        # subset a trajectory
        traj = self.traj.atom_slice(atom_indices=atom_indices, inplace=False)

        # extract the xyz coordinates
        #xyz = traj.xyz
        xyz = traj.xyz.reshape((traj.xyz.shape[0], traj.xyz.shape[1]*3))

        self.xyz = xyz

        return xyz


class ContactMap(object):
    """
    Construct a contact map with mdtraj distance matrix

    Parameters
    ----------
    traj: mdtraj.Trajectory object,
        the MDTraj trajectory object
    group_a: list,
        the list of atom index in cmap x-axis
    group_b: list,
        the list of atom index in cmap x-axis
    cutoff: float,
        the distance cutoff, default is 3.5 angstrom

    Attributes
    ----------
    atom_group_a: list,
        the list of atom index for cmap x-axis
    atom_group_b: list,
        the list of atom index for cmap y-axis
    atom_pairs_: ndarray, shape=[N, 2]
        the atom pairs for distance calculation, N is the number
        of atom pairs.

    Methods
    -------
    generate_atom_pairs
    distance_matrix
    generate_cmap

    """

    def __init__(self, traj, group_a, group_b, cutoff=0.35):

        self.traj = traj
        self.cutoff = cutoff

        self.atom_group_a = group_a
        self.atom_group_b = group_b

        self.atom_pairs_ = None
        self.generate_atom_pairs()

        self.dist_matrix_ = None
        self.cmap_ = None

        self.distmtx_computed_ = False
        self.cmap_computed_ = False

        self.coord_number_ = None

    def distance_matrix(self):
        """
        Calculate the atom distance matrix.

        Returns
        -------
        self.dist_matrix_: np.ndarray, shape=[A, B]
            A is number of atoms in rec, B is number of atoms in lig

        """
        if not self.distmtx_computed_:
            distmtx = mt.compute_distances(self.traj, atom_pairs=self.atom_pairs_)

            self.dist_matrix_ = distmtx
            self.distmtx_computed_ = True

        return self

    def generate_cmap(self, shape='array'):
        """
        Calculate atom cmap data.

        Parameters
        ----------
        shape: str,
            the type of cmap output, array like or matrix like.
            Default is array. Options are: array, matrix

        Returns
        -------
        self.cmap_: np.ndarray, shape=[N, A*B]
            the output cmap data, with N samples, and A*B elements per sample

        """
        if not self.distmtx_computed_:
            self.distance_matrix()

        if not self.cmap_computed_:
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
        """
        Generate atom pairs list.

        Returns
        -------
        self: the class itself

        """
        atom_pairs = []

        for a in self.atom_group_a:
            for b in self.atom_group_b:
                atom_pairs.append([a, b])

        self.atom_pairs_ = pd.DataFrame(atom_pairs).values

        return self

    def coord_num(self):
        # TODO: calculation coordination numbers
        if not self.distmtx_computed_:
            self.generate_cmap(shape="array")

        coord_number = np.sum(self.cmap_, axis=1)

        self.coord_number_ = coord_number

        return self


class DistanceMatrixMap(ContactMap):

    def __init__(self, traj, group_a, group_b, cutoff):
        ContactMap(traj, group_a, group_b, cutoff)

    def distance_matrix(self):
        if not self.distmtx_computed_:
            self.distance_matrix()

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


class CommunityCmap(object):

    def __init__(self, cmap, sep=","):

        if os.path.isfile(cmap):
            try:
                self.cmap = np.loadtxt(cmap, delimiter=sep, comments="#")
            except FileExistsError:
                print("File %s not exists! " % cmap)
        else:
            self.cmap = cmap

    def icriticalMap(self, icritical, dat):
        """
        Calculate whether the icritical number is suitable

        Parameters
        ----------
        icritical
        dat: np.ndarray, shape = [N, M]

        Returns
        -------
        True, or False
        """

        map = np.sum(np.greater(dat, icritical), axis=0) / dat.shape[0]

        return (map >= 0.48) * 1.0

    def diagonalZeroed(self, cmap, xsize, outfile):
        from mdanaly import matrix

        map = cmap.reshape((xsize, cmap.shape[0] / xsize))

        # matrix 2 xyz
        xyz = matrix.MatrixHandle().matrix2xyz(map)
        mtx = matrix.MatrixHandle().neiborhood2zero(xyz, neiborsize=4, outtype='mtx')

        np.savetxt(outfile, mtx, fmt="%3.1f", delimiter=" ")

        return mtx

    def generateCmap(self):

        #dat = np.loadtxt("all_cmapnbyn.csv", delimiter=",")
        Icritical = np.arange(0, 15, 1.0)

        for ic in Icritical:
            imap = self.icriticalMap(ic, self.cmap)

            self.diagonalZeroed(imap, int(np.sqrt(self.cmap.shape[1])), "cmap_I_%.2f" % ic)

        return 1

    def calculateNbyN(self, pdbfile, dcutoff, res1, res2, ndxlist):

        cutoff = dcutoff ** 2

        cmap = mdanaly.ContactMap(pdbfile)

        crd1 = pio.coordinatesPDB().getAtomCrdByNdx(pdbfile, ndxlist[res1])
        crd2 = pio.coordinatesPDB().getAtomCrdByNdx(pdbfile, ndxlist[res2])

        t = cmap.residueContacts(resCrd1=crd1,
                                 resCrd2=crd2,
                                 distcutoff=cutoff,
                                 verbose=False,
                                 rank=0,
                                 NbyN=True
                                 )
        # print(t)
        return t

    def scatterFileList(self, ranksize, pdbFileList):
        load4each = int(math.ceil(float(len(pdbFileList)) / float(ranksize)))
        filesList = []

        for i in range(ranksize - 1):
            filesList.append(pdbFileList[i * load4each: load4each * (i + 1)])
        filesList.append(pdbFileList[(ranksize - 1) * load4each:])

        return filesList

    def mpiCmapNbyN(self, MAX, inpdb_prefix="S_%s.pdb",
                    chain=" ", resindex =  [str(x) for x in range(1, 251)],
                    cutoff=5.0, atomtype = "side-chain-noH"):

        os.chdir(os.getcwd())

        startTime = datetime.now()

        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        rank = comm.rank

        fileList = [ inpdb_prefix % str(x) for x in range(1, MAX)]

        ndxdict = defaultdict(list)
        for res in resindex:
            ndxdict[res] = ndx.PdbIndex().res_index(fileList[0], chain,
                                                    residueNdx=[int(res)],
                                                    atomtype=atomtype,
                                                    atomList=[],
                                                    )
        if rank == 0:
            pdbFileList = self.scatterFileList(size, fileList)
        else:
            pdbFileList = None

        filesList = comm.scatter(pdbFileList, root=0)

        results = []

        for fn in filesList:
            count = 0
            print("progress file name {}, number {} out of ".format(fn, count, len(filesList)))

            nbyn = [ self.calculateNbyN(fn, cutoff, x, y, ndxdict) for x in resindex for y in resindex]
            # pair = [ [x, y] for x in resindex for y in resindex ]

            results.append(nbyn)
            count += 1

        np.savetxt(str(rank) + "_res_sidechain_cmap_nbyn.csv", np.array(results), delimiter=',', fmt="%5.3f")

        overallValuesList = comm.gather(results, root=0)
        if rank == 0:
            print("Total Time Usage: ")
            print(datetime.now() - startTime)


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
                                    "Default is: B 1 250 \n")
    parser.parser.add_argument('-cutoff',type=float,default=0.35,
                               help="Distance Cutoff for determining contacts. \n"
                                    "Default is 3.5 (angstrom). \n")
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
                                    "Choose the specific elements for atom indexing to construct the cmap."
                                    "Default is empty.")
    parser.parser.add_argument('-switch', type=str, default='True',
                               help="Input, optional. \n"
                                    "Apply a switch function for determing Ca-Ca contacts for a smooth transition. \n"
                                    "Only work with atomtype as CA. Options: True(T, t. TRUE), False(F, f, FALSE) \n"
                                    "Default is False. \n")
    parser.parser.add_argument('-np', default=0, type=int,
                               help='Number of Processers for MPI. Interger value expected. \n'
                                    'If 4 is given, means using 4 cores or processers.\n'
                                    'If 1 is given, means not using MPI, using only 1 Core.\n'
                                    'Default is 1. ')
    parser.parser.add_argument('-test', default=0, type=int,
                               help="Do a test with only a number of frames. For example, 4 frames. \n"
                                    "Default value is 0. ")
    parser.parser.add_argument('-NbyN', type=bool, default=False,
                               help="For community analysis, calculate atom contact number, normalized. \n"
                                    "Default is False.")
    parser.parser.add_argument('-verbose', default=False , type=bool,
                               help="Verbose. Default is False.")
    parser.parser.add_argument('-details', default=None, type=str,
                               help="Provide detail contact information and write out to a file. \n"
                                    "Default is None.")
    parser.parser.add_argument('-opt', default="TimeSeries", type=str,
                               help="Optional setting controls. Default is TimeSeries. \n"
                                    "TimeSeries, using splited files and get average cmap.\n"
                                    "Separated, create time series contact map, suitable for ligand "
                                    "protein contact information.\n"
                                    "Options: S(Separated), TS (TimeSeries).\n")

    parser.parse_arguments()

    return parser.args


def iterload_cmap():

    # for calculation time counting
    startTime = datetime.now()

    # argument options
    args = arguments()

    contact_map = np.array([])

    # TODO: atom selection method required
    indx = ndx.PdbIndex()
    atomList, atomType = indx.atomList(args.atomtype[0], atomname=args.atomname1)
    group_a = indx.res_index(args.s, args.rc[0], atomType, [int(args.rc[1]), int(args.rc[2])], atomList,)
    #group_a = np.array(group_a)
    atomList, atomType = indx.atomList(args.atomtype[1], atomname=args.atomname2)
    group_b = indx.res_index(args.s, args.lc[0], atomType, [int(args.lc[1]), int(args.lc[2])], atomList, )

    group_a = [int(x) - 1 for x in group_a]
    group_b = [int(x) - 1 for x in group_b]

    rec_index = int(args.rc[2]) - int(args.rc[1]) + 1
    lig_index = int(args.lc[2]) - int(args.lc[1]) + 1

    # read gromacs trajectory
    trajs = gmxcli.read_xtc(args.f, args.s, chunk=100, stride=int(args.dt/args.ps))

    for i, traj in enumerate(trajs):
        # calculate cmap information
        contmap = ContactMap(traj, group_a, group_b, cutoff=args.cutoff)
        contmap.generate_cmap(shape="array")

        if i == 0:
            contact_map = contmap.cmap_
        else:
            contact_map = np.concatenate((contact_map, contmap.cmap_), axis=0)

    # subset the results
    contact_map = pd.DataFrame(contact_map)
    contact_map.index = np.arange(contact_map.shape[0]) * args.dt
    contact_map = pca.datset_subset(contact_map, args.b, args.e)

    # get mean cmap data
    results = np.mean(contact_map, axis=0).reshape((rec_index, lig_index))

    # save results to an output file
    results = pd.DataFrame(results)
    results.index = range(int(args.rc[1]), int(args.rc[2]) + 1)
    results.columns = range(int(args.lc[1]), int(args.lc[2]) + 1)

    results.to_csv(args.o, sep=",", header=True, index=True, float_format="%.3f")

    print("Total Time Usage: ")
    print(datetime.now() - startTime)

