# -*- coding: utf-8 -*-

import mdtraj as mt
import numpy as np
import pandas as pd
from dockml import index
import sys
import os
from mdanaly import gmxcli


class ComputeAngles(object):
    """Calculate angles or dihedral angles of a trajectory

    Parameters
    ----------
    traj: mdtraj.trajectory object
        the mdtraj trajectory object, including coordinates
        and topology information

    Attributes
    ----------
    traj: mt.Trajectory,
        the input trajectory object (from mdtraj.Trajectory object)

    Methods
    -------
    get_angles(angle_index)
        Calculate angles between three atoms for a trajectory
    get_dihedral_angles(angle_index)
        Calculate dihedral angles between four atoms for a trajectory

    Examples
    --------
    Calculate the dihedral angle between atom 0, 1, 2, 3
    >>> from mdanaly import angles
    >>> import mdtraj as mt
    >>> angl = angles.ComputeAngles(traj=traj)
    >>> angl.get_dihedral_angles([[0, 1, 2, 3],])
    array([[1.8847224]], dtype=float32)
    Calculate a list of angles
    >>> angl.get_angles([[0, 1, 3], [0, 2, 4], [3, 4, 5]])
    array([[0.89176244, 2.0307815 , 1.8160943 ]], dtype=float32)
    """

    def __init__(self, traj):
        self.traj = traj

    def get_angles(self, angle_index):
        """calculate angles between three atoms for a trajectory

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
        """calculate dihedral angles between four atoms for a trajectory

        Parameters
        ----------
        angle_index: list, shape=[n, 4]
            the atom indices for angle calculations

        Returns
        -------
        angles: np.ndarray, shape=[N, M]
            angles, N is number of frames, M is the number
            of the angles per frame

        """

        angles = mt.compute_dihedrals(self.traj, indices=angle_index)

        return angles


def read_index(ndx="index.ndx", angle_type="dihedral"):
    """read gromacs index file and get atom indices

    Parameters
    ----------
    ndx: str, default = "index.ndx"
        Input gromacs index file name
    angle_type: str, default = "dihedral"
        Input angle type parameter, options: angle, dihedral

    Returns
    -------
    elements: np.ndarray, shape=[n,4] or [n, 3]
        The atom indices for angle calculation. It is a 2-D array.
    """
    indexer = index.GmxIndex(index=ndx)

    print("Please select a group by entering the index of the group: ")
    for i, g in enumerate(indexer.groups):
        print("%3d    %s " % (i, g))
    #print("Your choice: ")
    n = int(input("Your choice:  "))
    print("You have selected: %s " % indexer.groups[n])

    elements = indexer.groupContent(indexer.groups[n])

    # gromacs index starts from 1, while python-style index starts from 0
    elements = [int(x)-1 for x in elements]

    # if dihedral, four atoms index should be defined
    if angle_type == "dihedral":
        elements = np.reshape(np.array(elements), (-1, 4))
    else:
        elements = np.reshape(np.array(elements), (-1, 3))

    return elements


def arguments():
    """Parse the gmx-style command arguments for angle calculations.

    Returns
    -------
    args: Argparser object,
        the argparse object holding the arguement information
    """

    d = """
    Calculate time-series angles using a xtc trajectory. This function is simply designed to simulation
    gmx angle module, but provide a direct and easy way to save the result.

    Examples:
    Print help information
    gmx_angles.py -h

    Calculate angles
    gmx_angles.py -f traj.xtc -n index.ndx -s reference.pdb -o angles.csv -type angle -dt 10 -cos 0

    Calculate dihedral angles
    gmx_angles.py -f traj.xtc -n index.ndx -s reference.pdb -o dihedral_angles.csv -type dihedral -dt 10 -cos 0

    """

    parser = gmxcli.GromacsCommanLine(d)
    parser.arguments()

    parser.parser.add_argument("-type", type=str, default="angle",
                               help="Input, optional. \n"
                               "The angle type for calculation. Options are \n"
                               "angle, dihedral. Default is angle. ")
    parser.parser.add_argument("-cos", type=lambda x: (str(x).lower() == "true"), default=False,
                               help="Input, optional. \n"
                               "Calculate the cosine values of the angles.\n"
                               "Options are True, False. Default is False. ")
    parser.parser.add_argument("-sin", type=lambda x: (str(x).lower() == "true"), default=False,
                               help="Input, optional. \n"
                               "Calculate the sine values of the angles.\n"
                               "Options are True, False. Default is False. ")

    parser.parse_arguments()
    args = parser.args

    return args


def write_angles(angles, fout, cosine, sine, dt=2, begin=0, end=-1):
    """write angles (unit: degree) into a file

    Parameters
    ----------
    angles : ndarray, shape=[N * M ]
        the angles, N is number of frames, M is the number of angles
    fout : str,
        output file name, comma separated
    cosine: int,
        whether save the cosine of the angles, default is 0
    sine : int,
        whether save the sine of the angles, default is 0
    dt : int, default=2
        the stride of the frames were saved or angles were calculated
    begin : int, default = 0,
        the beginning frame index
    end : int, default = -1
        the ending frame index

    """
    pi = 3.141592654

    if cosine and sine:
        angles_1 = np.cos(angles)
        angles_2 = np.sin(angles)
        new_angles = np.concatenate((angles_1, angles_2), axis=1)
    elif cosine and not sine:
        new_angles = np.cos(angles)
    elif not cosine and sine:
        new_angles = np.sin(angles)
    else:
        new_angles = (angles / pi) * 180

    dat = pd.DataFrame(new_angles)
    dat.index = np.arange(new_angles.shape[0]) * dt
    dat.columns = ["Angles_"+str(x) for x in np.arange(new_angles.shape[1])]

    if begin > 0:
        dat = dat[dat.index >= begin]
    if end > 0 and end > begin:
        dat = dat[dat.index <= end]

    dat.to_csv(fout, sep=",", float_format="%.3f", header=True, index=True)
    # np.savetxt(fout, angles, fmt="%.3f", delimiter=",")

    return None


def gmxangle():
    """A gromacs g_angle simulator which works the same way as the gromacs tool

    Returns
    -------
    angles: np.ndarray, shape=[N * M ]
        the angles, N is number of frames, M is the number of angles per frame

    """

    args = arguments()

    if os.path.exists(args.f) and os.path.exists(args.n) \
            and os.path.exists(args.s):

        # prepare index atom slices
        ndx = read_index(args.n, args.type)

        if args.v:
            print("Compute cosine of the angles: ", bool(args.cos))
            print("Atom indices: ")
            print(ndx)

        # load trajectories
        trajs = gmxcli.read_xtc(xtc=args.f, top=args.s, chunk=1000,
                                stride=int(args.dt / args.ps))
        if args.v:
            print("Frame information: ")
            for i, traj in enumerate(trajs):
                print("Trajectory %3d: %12d frames" % (i, traj.n_frames))

        angles = np.array([])

        # for each traj chunk, calculate angles, and cat them together
        for i, traj in enumerate(trajs):
            cangle = ComputeAngles(traj)
            if args.v:
                print("Progress: %12d " % (i * traj.n_frames))
                print(angles.shape)

            if angles.shape[0] == 0:
                angles = cangle.get_dihedral_angles(ndx)
            else:
                angles = np.concatenate((angles, cangle.get_dihedral_angles(ndx)), axis=0)

        if args.v:
            print("Write angles to output file: ", args.o)

        # write angles to an output file
        write_angles(angles, args.o, cosine=args.cos, sine=args.sin, dt=args.dt, begin=args.b, end=args.e)

    else:
        print("Some of the input files are not existed. Input again.")
        if not os.path.exists(args.n):
            print("Index file is not provided. Run again.")
        if not os.path.exists(args.f):
            print("Xtc trajectory file is not provided. Run again.")
        if not os.path.exists(args.s):
            print("Reference file is not provided. Run again.")

    sys.exit(1)

