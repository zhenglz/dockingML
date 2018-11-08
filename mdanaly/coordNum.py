# -*- coding: utf-8 -*-

from dockml import index
from mdanaly import gmxcli
from mdanaly import cmap
import pandas as pd
import numpy as np


def arguments():
    d = """
    Calculate coordination number between molecules. 
    
    Examples:
    Caculate the coordination number between two parts in a system only considering non hydrogen atoms
    gmx_coordnum.py -f test_traj_dt1ns.xtc -s reference.pdb -o coordnum -rc " " 1 41 
                    -lc " " 42 137 -atomtype heavy heavy -dt 4
    
    """

    parser = gmxcli.GromacsCommanLine(d=d)
    parser.arguments()

    parser.parser.add_argument("-rc", type=str, default=[], nargs="+",
                               help="Input, optional. \n"
                                    "The chain identifier, start res index, end res index. ")
    parser.parser.add_argument("-lc", type=str, default=[], nargs="+",
                               help="Input, optional. \n"
                                    "The chain identifier, start res index, end res index. ")
    parser.parser.add_argument("-atomtype", type=str, default=[], nargs="+",
                               help="Input, optional. \n"
                                    "The atomtype used for contact calculation. Options: mainchain, sidechain, \n"
                                    "heavy, CA. ")
    parser.parser.add_argument("-cutoff", type=float, default=0.5,
                               help="Input, optional. \n"
                                    "The distance cutoff for coordination number calculation."
                                    "Unit is nanometer, default is 0.5 ")

    parser.parse_arguments()
    args = parser.args

    return args


def run_coord_number():

    args = arguments()

    # TODO: define a way to select atom slices
    # define the atom indices
    indx = index.PdbIndex()
    atomList, atomType = indx.atomList(args.atomtype[0], atomname=[])
    group_a = indx.res_index(args.s, args.rc[0], atomType, [int(args.rc[1]), int(args.rc[2])], atomList,)

    atomList, atomType = indx.atomList(args.atomtype[1], atomname=[])
    group_b = indx.res_index(args.s, args.lc[0], atomType, [int(args.lc[1]), int(args.lc[2])], atomList, )

    group_a = [int(x) - 1 for x in group_a]
    group_b = [int(x) - 1 for x in group_b]

    results = np.array([])

    print("Loading trajectory xtc file ...... ")
    trajs = gmxcli.read_xtc(args.f, args.s, chunk=100, stride=int(args.dt/args.ps))

    print("Performing calculations ...... ")
    for i, traj in enumerate(trajs):
        coord_num = cmap.ContactMap(traj, group_a, group_b, args.cutoff)
        coord_num.coord_num()

        if i == 0:
            results = coord_num.coord_number_
        else:
            results = np.concatenate((results, coord_num.coord_number_), axis=0)

    results = pd.DataFrame(results)
    results.index = np.arange(results.shape[0]) * args.dt

    print("Saving results to output file ...... ")
    results.to_csv(args.o, sep=",", header=False, index=True, float_format="%.1f")

