#!/usr/bin/env python

import numpy as np
import sys

class PMF :
    def __init__(self):
        pass

    def pmf2d(self, data, nbins=20, xcol=0, ycol=1):
        """
        give a 2d array, calculate the 2d distribution
        and then calculate its PMF fes
        :param data:
        :param nbins:
        :param xcol:
        :param ycol:
        :return:
        """
        data = np.array(data)
        X1 = list(data[:, xcol])
        X2 = list(data[:, ycol])

        hist1, edges1 = np.histogram(X1, bins=nbins)
        hist2, edges2 = np.histogram(X2, bins=nbins)

        n = len(edges1)
        pmf = {}
        for j in range(n - 1):
            for k in range(n - 1):
                pmf[str(j) + "_" + str(k)] = 0

        for i in range(len(X1)):
            x = [X1[i], X2[i]]
            for k in range(n - 1):
                for j in range(n - 1):
                    if x[0] >= edges1[k] and x[0] < edges1[k + 1]:
                        if x[1] >= edges2[j] and x[1] < edges2[j + 1]:
                            pmf[str(k) + "_" + str(j)] += 1.0

        pmf_array = []
        for j in range(n - 1):
            for k in range(n - 1):
                z = [edges1[j], edges2[k], pmf[str(k) + "_" + str(j)]]
                pmf_array.append(z)

        return np.asarray(pmf_array)

def main() :
    if len(sys.argv) < 2 :
        print('''
        Transform a 2D data (time series) file into a PMF matrix
        (Only numerical data files are accepted)
        Usage:
        python pmf2d.py your2d.data numberofbins
        python pmf2d.py proj2d.dat
        python pmf2d.py proj2d.dat 40
        ''')
        sys.exit(0)

    dataf = sys.argv[1]
    if len(sys.argv) == 3 :
        nb = int(sys.argv[2])
    else :
        nb = 20
    df = np.loadtxt(dataf, comments=["#", "@"] )#, usecols=[0, 1])

    pmf = PMF()
    dist = pmf.pmf2d(df, nb)

    np.savetxt("pmf_"+dataf, dist, fmt="%.3f")
