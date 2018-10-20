# -*- coding: utf-8 -*-

import sys, os
import mdtraj as mt
import pandas as pd
import numpy as np
import argparse
from argparse import RawTextHelpFormatter


class GromacsCommanLine(object):
    """
    Gromacs style commond line argument API

    Parameters
    ----------
    d: str,
        the description and example information

    Attributes
    ----------
    args: argparse arguments object
    parser
    """

    def __init__(self, d=""):

        self.description = d

        self.args = None

        self.parser = None

    def arguments(self):
        """
        Gromacs style arguments object

        Returns
        -------
        args: Argparser object,
            the argparse object holding the arguement information
        """

        parser = argparse.ArgumentParser(description=self.description,
                                         formatter_class=RawTextHelpFormatter)

        parser.add_argument("-f", type=str, default="md.xtc",
                            help="Input. Options: .xtc, .pdb \n"
                                 "The xtc trajectory file name. ")
        parser.add_argument("-s", type=str, default="reference.pdb",
                            help="Input. Options: .pdb \n"
                                 "Reference pdb file, where topology information holds. ")
        parser.add_argument("-n", type=str, default="index.ndx",
                            help="Input, optional. Options: .ndx \n"
                                 "Gromacs type index file, where atom indices information \n"
                                 "holds for angle calculation. Default is index.ndx. ")
        parser.add_argument("-o", type=str, default="output.csv",
                            help="Output. Options: .dat, .xvg, .csv \n"
                                 "The output file name. \n"
                                 "Default is output.csv, with header and index, \n"
                                 "as well as the comma separator. ")
        parser.add_argument("-dt", type=int, default=2,
                            help="Input, optional. \n"
                                 "Skip frames with a gap of dt picoseconds. \n"
                                 "Default is 2. ")
        parser.add_argument("-ps", default=2, type=int,
                            help="Input, optional. \n"
                                 "How many picoseconds the frames are stored in\n"
                                 "trajectory file. Default is 2. ")
        parser.add_argument("-v", default=False, type=lambda x: (str(x).lower() == "true"),
                            help="Input, optional. \n"
                                 "Whether print detail information. Default is False. ")
        parser.add_argument("-b", default=0, type=int,
                            help="Input, optional. \n"
                                 "The beginning frame of the calculation. Default is 0. ")
        parser.add_argument("-e", default=-1, type=int,
                            help="Input, optional. \n"
                                 "The ending frame of the calculation. Default is -1. ")

        self.parser = parser

        return self

    def parse_arguments(self):
        """
        parse the arguments

        Returns
        -------
        self: object
            the class object

        """
        self.args = self.parser.parse_args()

        # print help information
        if len(sys.argv) < 3:
            self.parser.print_help()
            sys.exit(0)

        return self


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

