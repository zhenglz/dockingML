# -*- coding: utf-8 -*-

from dockml import index
from dockml import pdbIO
import mdtraj as mt
import numpy as np
from mdanaly import gmxcli
from mdanaly import cmap
import pandas as pd


def arguments():
    d = """
    Calculate coordination number between molecules. 
    
    Examples:
    
    
    """

    parser = gmxcli.GromacsCommanLine(d=d)
    parser.arguments()

    parser.parser.add_argument("-rec", type=str, default=[],
                               help="Input, optional. ")
    parser.parser.add_argument("-lig", type=str, default=[],
                               help="Input, optional. ")
    parser.parser.add_argument("-cutoff", type=float, default=0.5,
                               help="Input, optional. \n"
                                    "The distance cutoff for coordination number calculation."
                                    "Unit is nanometer, default is 0.5 ")

    args = parser.parse_arguments()

    return args


def run_coord_number():

    args = arguments()

    # TODO: define a way to select atom slices
    group_a = []
    group_b = []

    results = np.array([])

    trajs = gmxcli.read_xtc(args.o, args.s, chunk=1000, stride=int(args.dt/args.ps))

    for traj in trajs:

        coord_num = cmap.CoordinationNumber(traj, group_a, group_b, args.cutoff)

        if results.shape[0] == 0:
            results = coord_num.coord_number_
        else:
            results = np.concatenate((results, coord_num.coord_number_), axis=1)

    results = pd.DataFrame(results)
    results.index = np.array(results.shape[0]) * args.dt

    results.to_csv(args.o, sep=",", header=False, index=True)

