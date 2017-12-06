#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
basic math algorithms
"""

import math
import numpy as np

class BasicAlgorithm :

    def __init__(self):
        pass

    def switchFuction( self, x, d0=7.0, m=12, n=6):
        """
        for countting, implement a rational switch function to enable a smooth transition
        the function is lik  s= [1 - (x/d0)^6] / [1 - (x/d0)^12]
        d0 is a cutoff, should be twice the larget than the distance cutoff
        :param x: float
        :param d0: distance cutoff, should be 2 times of normal cutoff
        :param m: int
        :param n: int
        :return: float
        """
        count = 0.0
        try:
            count = (1.0 - math.pow((x / d0), n)) / (1.0 - math.pow((x / d0), m))
        except ZeroDivisionError:
            print("Divide by zero, ", x, d0)

        return count

        #return (1.0 - math.pow((x / d0), n)) / (1.0 - math.pow((x / d0), m))

    def pmf(self, x, minX, kt=2.5, max=1.0):
        """
        calculate PMF of a histogram vector
        :param x:
        :param minX: float, avoid divide by zero problem
        :param kt: float, kt=2.5 when T=300K and unit is kJ/mol
        :param max:
        :return: list of floats
        """
        if x < minX:
            x = minX / 2.0
        return -1.0 * kt * np.log(x / max)

    def entropy1D(self, x, nbins=20, kbt=2.5):
        """
        get entropy from 1d vector
        :param x:
        :param nbins:
        :param kbt:
        :return:
        """

        # calculate histogram of the vector x
        hist, edges = np.histogram(x, bins=nbins)
        # compute the probability distribution of the histogram
        prob = hist / np.asarray(x).shape[0]

        prob = [ x for x in list(prob) if x > 0 ]

        entropy = -1.0 * kbt * sum([ x * math.log(x) for x in prob ])

        return entropy

    def entropy2D(self, x, y, nbins=20, kbt=2.5):
        """
        get entropy from 2d vector
        :param x:
        :param y:
        :param nbins:
        :param kbt:
        :return:
        """

        # calculate histogram of the vectors
        hist, xedges, yedges = np.histogram2d(x, y, bins=nbins)
        # compute the probability distribution of the histogram
        prob = hist / np.asarray(x).shape[0]

        prob = [ x for x in list(prob.ravel()) if x > 0 ]

        entropy = -1.0 * kbt * sum([x * math.log(x) for x in prob ])

        return entropy
