# -*- coding: utf-8 -*-

from dockml import index
from mdanaly import gmxcli, cmap

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
                               help="Input, optional. Default = 0.5 \n"
                                    "The distance cutoff for coordination number calculation. \n"
                                    "Unit is nanometer.")
    parser.parser.add_argument("-byres", type=lambda x: (str(x).lower() == "true"), default=False,
                               help="Input, optional. Default = False. \n"
                                    "Computate contact number per residue. ")

    parser.parse_arguments()
    args = parser.args

    return args


def run_coord_number():

    args = arguments()

    resides_a, resides_b = [], []
    group_a, group_b = [], []

    if args.byres:
        ndx = index.PdbIndex(args.s, [args.rc[0]], args.rc[1:])
        s, e = ndx.resid_mt_style(args.rc[0], args.rc[1], args.rc[2])
        resides_a = np.arange(s, e+1)

        ndx = index.PdbIndex(args.s, [args.lc[0]], args.lc[1:])
        s, e = ndx.resid_mt_style(args.lc[0], args.lc[1], args.lc[2])
        resides_b = np.arange(s, e+1)

    else:
        # TODO: define a way to select atom slices
        # define the atom indices for receptor
        ndx = index.PdbIndex(reference=args.s, atomtype=args.atomtype[0],
                             resSeq=args.rc[1:],
                             chain=[args.rc[0]])
        ndx.prepare_selection()
        ndx.res_index()
        group_a = ndx.atomndx_mt_style_

        # define the atom indices for ligand
        ndx = index.PdbIndex(reference=args.s, atomtype=args.atomtype[1],
                             resSeq=args.lc[1:],
                             chain=[args.lc[0]])
        ndx.prepare_selection()
        ndx.res_index()
        group_b = ndx.atomndx_mt_style_

    results = np.array([])

    print("Loading trajectory xtc file ...... ")
    trajs = gmxcli.read_xtc(args.f, args.s, chunk=100,
                            stride=int(args.dt/args.ps))

    print("Performing calculations ...... ")
    for i, traj in enumerate(trajs):

        if args.byres:
            coord = cmap.CmapNbyN(traj, resids_a=resides_a,
                                  resids_b=resides_b, cutoff=args.cutoff)
            coord.contact_nbyn()
            if i == 0:
                results = coord.contacts_by_res_
            else:
                results = np.concatenate((results, coord.contacts_by_res_), axis=0)

        else:
            coord = cmap.ContactMap(traj, group_a, group_b, args.cutoff)
            coord.coord_num()

            if i == 0:
                results = coord.coord_number_
            else:
                results = np.concatenate((results, coord.coord_number_), axis=0)

    results = pd.DataFrame(results)
    results.index = np.arange(results.shape[0]) * args.dt

    print("Saving results to output file ...... ")
    results.to_csv(args.o, sep=",", header=False, index=True, float_format="%.1f")

