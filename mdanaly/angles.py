# -*- coding: utf-8 -*-

import mdtraj as mt
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
from dockml import index
import sys
import os


class ComputeAngles(object):
    """
    Calculate angles or dihedral angles of a trajectory

    Parameters
    ----------
    traj: mdtraj.trajectory object
        the mdtraj trajectory object, including coordinates
        and topology information

    Attributes
    ----------

    """

    def __init__(self, traj):
        self.traj = traj

    def get_angles(self, angle_index):
        """
        calculate angles between three atoms for a trajectory

        Parameters
        ----------
        angle_index: list, shape=[n, 3]
            the atom indices for angle calculations

        Returns
        -------
        angles: ndarray, shape=[N, M]
            angles, N is number of frames, M is the number
            of the angles per frame

        """

        angles = mt.compute_angles(self.traj, angle_indices=angle_index)

        return angles

    def get_dihedral_angles(self, angle_index):
        """
        calculate dihedral angles between four atoms for a trajectory

        Parameters
        ----------
        angle_index: list, shape=[n, 4]
            the atom indices for angle calculations

        Returns
        -------
        angles: numpy ndarray, shape=[N, M]
            angles, N is number of frames, M is the number
            of the angles per frame

        """

        angles = mt.compute_dihedrals(self.traj, indices=angle_index)

        return angles


def read_xtc(xtc, top, chunk=100, stride=1):
    """

    Parameters
    ----------
    xtc: str,
        input xtc file name
    top: str,
        input topology information file, a pdb
    chunk: int,
        number of frame per chunk
    stride: int,
        dt, save a frame every N number of frames

    Returns
    -------
    trajs: list,
        a list of mdtraj trajectory object
    """

    trajs = []

    for chunk in mt.iterload(xtc, chunk=chunk, top=top, stride=stride):
        trajs.append(chunk)

    print("Number of chunks: ", len(trajs))

    return trajs


def read_index(ndx, angle_type):
    """
    read gromacs index file and get atom indices

    Parameters
    ----------
    ndx: str,
        input gromacs index file
    angle_type: str,
        input angle type parameter, options: angle, dihedral

    Returns
    -------
    elements: ndarray, shape=[n,4] or [n, 3]
        the atom indices for angle caculation
    """
    indexer = index.GmxIndex(index=ndx)

    print("Please select a group by entering the index of the group: ")
    for i, g in enumerate(indexer.groups):
        print("%3d    %s " % (i, g))
    #print("Your choice: ")
    n = int(input("Your choice:  "))
    print("You have selected: %s" % indexer.groups[n])

    elements = indexer.groupContent(indexer.groups[n])

    elements = [int(x)-1 for x in elements]
    # if dihedral, four atoms index should be defined
    if angle_type == "dihedral":
        elements = np.reshape(np.array(elements), (-1, 4))
    else:
        elements = np.reshape(np.array(elements), (-1, 3))

    return elements

def arguments():
    """

    Returns
    -------
    args: Argparser object,
        the argparse object holding the arguement information
    """

    d = """
    Calculate angles of a xtc trajectory. This function is simply designed to simulation
    gmx angle module, but provide a direct way to save the result.
    
    Examples:
    Print help information
    gmx_angles.py -h
    
    Calculate angles
    gmx_angles.py -f traj.xtc -n index.ndx -s reference.pdb -o angles.csv -type angle -dt 10 -cos False
    
    Calculate dihedral angles
    gmx_angles.py -f traj.xtc -n index.ndx -s reference.pdb -o dihedral_angles.csv -type dihedral -dt 10 -cos False
    
    """

    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)

    parser.add_argument("-f", type=str, default="md.xtc",
                        help="Input. The xtc trajectory file name. ")
    parser.add_argument("-s", type=str, default="reference.pdb",
                        help="Input. Reference pdb file, where topology information holds. ")
    parser.add_argument("-n", type=str, default="index.ndx",
                        help="Input. Gromacs type index file, where atom indices information "
                             "holds for angle calculation.")
    parser.add_argument("-o", type=str, default="angles.csv",
                        help="Output. The output angle file name. Default is angle.csv. ")
    parser.add_argument("-type", type=str, default="angle",
                        help="Input, optional. The angle type for calculation. Options are "
                             "angle, dihedral. Default is angle. ")
    parser.add_argument("-cos", type=bool, default=False,
                        help="Input, optional. Calculate the cosine values of the angles. "
                             "Options are True, False. Default is False. ")
    parser.add_argument("-dt", type=int, default=2,
                        help="Input, optional. Skip frame with a gap of dt frames. "
                             "Default is 0. ")
    parser.add_argument("-ps", default=2, type=int,
                        help="Input, optional. How many picoseconds the frames are stored in"
                             "trajectory file. Default is 2. ")
    parser.add_argument("-v", default=False, type=bool,
                        help="Input, optional. Whether print detail information. "
                             "Default is False. ")

    args = parser.parse_args()

    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(0)

    return args


def write_angles(angles, fout, cosine=True):
    """
    write angles into a file

    Parameters
    ----------
    angles: ndarray, shape=[N * M ]
        the angles, N is number of frames, M is the number of angles
    fout: str,
        output file
    cosine: bool,
        whether save the cosine of the angles

    Returns
    -------

    """
    pi = 3.1415926

    if cosine:
        angles = np.cos(angles)
        #np.savetxt(fout, angles, fmt="%.3f", delimiter=",")
        #return angles
    else:
        angles = (angles / pi) * 180

    np.savetxt(fout, angles, fmt="%.3f", delimiter=",")

    return angles


def gmxangle(args):
    """
    A gromacs g_angle simulator which works the same way as the gromacs tool

    Parameters
    ----------
    args: argparse object,
        the arguements for input and output as well as the parameters

    Returns
    -------
    angles: ndarray, shape=[N * M ]
        the angles, N is number of frames, M is the number of angles per frame

    """

    if os.path.exists(args.f) and os.path.exists(args.n) and os.path.exists(args.s):

        # prepare index atom slices
        ndx = read_index(args.n, args.type)

        if args.v:
            print("Atom indices: ")
            print(ndx)

        # load trajectories
        trajs = read_xtc(xtc=args.f, top=args.s, chunk=1000, stride=int(args.dt / args.ps))
        if args.v:
            print("Frame information: ")
            for i, traj in enumerate(trajs):
                print("Trajectory %3d: %12d frames" % (i, traj.n_frames))

        angles = np.array([])

        for i, traj in enumerate(trajs):
            cangle = ComputeAngles(traj)
            if args.v:
                print("Progress: %12d " % (i * traj.n_frames))
                print(angles.shape)

            if angles.shape[0] == 0:
                angles = cangle.get_dihedral_angles(ndx)
            else:
                angles = np.concatenate((angles, cangle.get_dihedral_angles(ndx)), axis=0)
            #print("Progress: %12d " % (i * traj.n_frames))

        # write angles to an output file
        angles = write_angles(angles, args.o, cosine=args.cos)

        return angles

    else:
        print("Some of the input file is not existed. Input again.")
        return np.array([])

