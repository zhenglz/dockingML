#!/usr/bin/env python

import os

import numpy as np
import pandas as pd

from automd import fixpdb

if __name__ == "__main__" :

    MAX = 2001

    cmap = []

    if not os.path.exists("TimeSeries_Cmap_heavy.csv") :

        for i in range(1, MAX+1) :

            filen = "all_cmap_%d.xyz" % i

            c = np.loadtxt(filen, comments="#", usecols=([2]))

            cmap.append(c)

        cmap = np.array(cmap).T

        np.savetxt("TimeSeries_Cmap_heavy.csv", cmap, fmt="%3.1f", delimiter=",")
    else :
        cmap = np.loadtxt("TimeSeries_Cmap_heavy.csv", delimiter=",")

    cmap = pd.DataFrame(cmap)

    res_ndx = []
    with open( "all_cmap_%d.xyz" % 1) as lines :
        res_ndx = sorted(set([ int(x.split()[0].split("_")[0]) for x in lines if "#" not in x ]))

    cmap['Index'] = np.array(res_ndx)

    fpdb = fixpdb.SummaryPDB("all_cmap_%d.xyz" % 1, "")

    x, y, z, m, n, resname_ndx = fpdb.details()

    resname = []
    for item in res_ndx :
        t = str(item) +"_A"
        resname = resname_ndx[t] + str(item)

    cmap['resname'] = np.array(resname)

    cmap.to_csv("TimeSeries_Map.csv", sep=",")
