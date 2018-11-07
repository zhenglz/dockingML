# -*- coding: utf-8 -*-

from mdanaly import gmxcli
import numpy as np
import os


class DsspParser(object):

    def __init__(self, dssp_xpm, dt=2, ps=2, res_ndx=[0, 1, 2, 3]):
        self.input = dssp_xpm
        self.df = None
        self.dt = dt
        self.ps = 2

        self.res_range = res_ndx

    def read_dssp(self):
        dat = []
        #self.df = np.loadtxt(self.input, comments=["/", "#", "*"])
        if os.path.exists(self.input):
            with open(self.input) as lines:
                for s in lines:
                    if "/*" not in s:
                        dat.append(list(s))
        self.df = np.array(dat)

        return self

    def dssp_part(self, ):

        if not self.df:
            self.read_dssp()

        dat = self.df[self.res_range][:, ::self.dt/self.ps]

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
    dat = dssp.dssp_part()

    np.savetxt(args.o, dat, delimiter=",",)

