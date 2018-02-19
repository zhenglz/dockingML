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

    def entropy1D(self, x, nbins=20, kbt=1.0):
        """
        get entropy from 1d vector
        :param x:
        :param nbins:
        :param kbt:
        :return:
        """
        #entropy = 0
        # calculate histogram of the vector x
        hist, edges = np.histogram(x, bins=nbins)
        # compute the probability distribution of the histogram
        prob = hist / float(np.asarray(x).shape[0])

        #print(prob)

        prob = [ x for x in list(prob) if x > 0 ]

        entropy = -1.0 * kbt * sum([ x * math.log(x) for x in prob ])

        return entropy

    def entropy2D(self, x, y, nbins=20, kbt=1.0):
        """
        get entropy from 2d vector
        :param x:
        :param y:
        :param nbins:
        :param kbt:
        :return:
        """
        entropy = 0
        # calculate histogram of the vectors
        hist, xedges, yedges = np.histogram2d(x, y, bins=nbins)
        # compute the probability distribution of the histogram
        prob = hist / np.asarray(x).shape[0]

        prob = [ x for x in list(prob.ravel()) if x > 0 ]

        entropy = -1.0 * kbt * sum([x * math.log(x) for x in prob ])

        return entropy

class PlaneFit :
    def __init__(self):
        pass

    def fitPlane(self, points):
        """
        fit some points to a plane
        https://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
        plane function: ax + by - z + c = 0
        or: ax + by - z = -c
        Ax + By + Cz + D = 0
        we need to determine [a, b, c]
        A B C D = a  b -1 c
        :param points: ndarray, M*3
        :return: array, [a, b, c]
        """

        # prepare dataset
        xs = np.array(points)[:, 0]
        ys = np.array(points)[:, 1]
        zs = np.array(points)[:, 2]

        tmp_A = []
        tmp_B = []
        for i in range(xs.shape[0]):
            tmp_A.append([xs[i], ys[i], 1])
            tmp_B.append(zs[i])

        B = np.matrix(tmp_B).T
        A = np.matrix(tmp_A)

        # do fit
        fit = (A.T * A).I * A.T * B
        errors = B - A * fit
        residual = np.linalg.norm(errors)

        return fit

    def point_distance(self, params, point):
        """
        https://mathinsight.org/distance_point_plane
        determine the distance between a point to a plane
        the function of a function Ax + By + Cz + D = 0
        P (x0, y0, z0 ) to the plane
        distance (x0*A + y0*B + z0*C + D)/sqrt(A^2 + B^2 + C^2)
        :param params: list, a list of parameters (A, B, D), while C=-1
        :param point: list, x y z coordinates of a point
        :return: float, distance
        """

        distance = math.sqrt(params[0] ** 2 + params[1] ** 2 + 1)

        distance = (params[0] * point[0] +
                    params[1] * point[1] +
                    (-1.0) * point[2] +
                    params[2]
                    ) / distance

        return distance

