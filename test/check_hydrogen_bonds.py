#!/usr/bin/env python

import sander
import parmed as pmd
#from .modeller import Modeller
import argparse
import subprocess as sp
import numpy as np
import os, sys
import itertools
import pandas as pd
import mdtraj as mt


class PDBExtractor(object):

    def __init__(self, traj_fn):
        self.traj_fn = traj_fn
        self.traj = None
        self.n_frames = 0

    def load_pdb(self, stride=4):
        traj = None
        for i, t in enumerate(mt.iterload(self.traj_fn, stride=stride, chunck=100)):
            if i == 0:
                traj = t
            else:
                traj = mt.join([traj, t])

        self.traj = traj.atom_slice(traj.topology.select("not water"))
        self.n_frames = traj.n_frames

        return self.traj

    def get_frame(self, frame_index, outpdb):
        if self.traj is None:
            self.load_pdb()

        self.traj[frame_index].save_pdb(outpdb)

        return self

    def extract_all_frames(self, base_name):
        if self.traj is None:
            self.load_pdb()

        for i in range(self.n_frames):
            self.get_frame(i, "%s_%d.pdb" % (base_name, i))

        return self


if __name__ == "__main__":

    fn = sys.argv[1]

    cutoff = 0.25

    traj = PDBExtractor(traj_fn=fn).load_pdb(1)

    acceptors = traj.topology.select('residue 6 and name N7 O6')

    donars1 = traj.topology.select('residue 26 and name H13 H14')
    donars2 = traj.topology.select('residue 26 and name H15 H16')

    pairs_1 = list(itertools.product(acceptors, donars1))
    pairs_2 = list(itertools.product(acceptors, donars2))

    distances1 = mt.compute_distances(traj, pairs_1)
    distances2 = mt.compute_distances(traj, pairs_2)

    counts1 = (distances1 <= cutoff) * 1.0
    counts2 = (distances2 <= cutoff) * 1.0

    type1, type2 = 0, 0
    for i in range(counts1.shape[0]):
        if counts1[i].sum() == 1 and counts2[i].sum() == 1:
            type1 += 1
        elif counts1[i].sum() == 2 or counts2[i].sum() == 2:
            type2 += 1

    print(type1, type2, traj.n_frames)


