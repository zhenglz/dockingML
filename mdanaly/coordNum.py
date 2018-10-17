# -*- coding: utf-8 -*-

from dockml import index
from dockml import pdbIO
import mdtraj as mt
import numpy as np
from mdanaly import gmxcli
import pandas as pd


class CoordinationNumber(object):

    def __init__(self, xtc, top="reference.pdb", dt=2, atom_pairs=None):
        self.xtc = xtc
        self.top = top
        self.dt = dt
        self.atom_pairs = atom_pairs

        # attributes
        self.xtc_loaded_ = False
        self.trajs = None
        self.distances_ = np.array([])

        self.coord_number_ = np.array([])

    def atom_slices(self):
        # TODO: get atom slices
        return NotImplementedError

    def iterload_xtc(self, chunk=1000):

        if not self.xtc_loaded_:
            self.trajs = []

            for chunk in mt.iterload(self.xtc, top=self.top, chunk=chunk, stride=self.dt):
                self.trajs.append(chunk)

            self.xtc_loaded_ = False
        else:
            pass

        return self

    def compute_pair_distances(self):

        for traj in self.trajs:

            distances = mt.compute_distances(traj, self.atom_pairs, True)

            if self.distances_.shape[0] == 0:
                self.distances_ = distances
            else:
                self.distances_ = np.concatenate((self.distances_, distances), axis=1)

        return self

    def compute_coord_number(self, cutoff=0.35):

        if self.distances_.shape[0] == 0:
            self.compute_pair_distances()

        self.coord_number_ = np.sum(self.distances_ <= cutoff, axis=1)

        return self


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

    parser.parse_arguments()


def parse_atom_indices():
    # TODO: parse the atom indices given the atom information

    return NotImplementedError


def coordination_number():
    # TODO: main entry point of this coordination number module

    return NotImplementedError