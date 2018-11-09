# -*- coding: utf-8 -*-

from mdanaly import gmxcli
import numpy as np
import os
import pandas as pd
import math


class DsspParser(object):

    def __init__(self, dssp_xpm, dt=2, ps=2, res_ndx=[0, 1, 2, 3]):
        self.input = dssp_xpm
        self.df = np.array([])
        self.dt = dt
        self.ps = 2

        self.res_range = res_ndx

    def read_dssp(self):
        dat = []

        if os.path.exists(self.input):
            with open(self.input) as lines:
                for s in lines:
                    if "/*" not in s and len(s.split()) <= 2:
                        dat.append(list(s.split(",")[0].strip("\"")))
        self.df = pd.DataFrame(dat).values

        return self

    def dssp_part(self, b=0, e=-1):

        if not self.df.shape[0]:
            self.read_dssp()
        try:
            dat = self.df[self.res_range][:, b:e:int(self.dt/self.ps)]
        except IndexError:
            dat = self.df[self.res_range][:, b::int(self.dt/self.ps)]

        return dat


def arguments():
    d = """
    """
    parser = gmxcli.GromacsCommanLine(d=d)
    parser.arguments()

    parser.parser.add_argument("-res", default=[1, 5], type=int, nargs="+",
                               help="Input, optional. \n"
                                    "The input residue start and end. \n")

    parser.parse_arguments()

    return parser.args


def parse_dssp():

    args = arguments()

    dssp = DsspParser(args.f, args.dt, args.ps, list(range(args.res[0]-1, args.res[1])))
    dssp.read_dssp()
    dat = dssp.dssp_part(b=int(args.b/args.ps), e=math.floor(args.e/args.ps))

    np.savetxt(args.o, dat, delimiter=",", fmt="%1s")

