#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
basic math algorithms
"""

import math
import numpy as np
import pandas as pd

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

    def exponentialFunction(self, x, exp=2.0, k=1.0, c0=0.0, bias=0.0):
        '''

        :param x: variable
        :param exp:  float, exponential number
        :param k: float, constant
        :param c0: float
        :param bias: float, bias or intercept
        :return: float, result
        '''

        return k * (x + c0) ** exp + bias

    def pmf(self, x, minX, kbt=2.5, max=1.0):
        """
        calculate PMF of a histogram vector
        :param x: float
        :param minX: float, avoid divide by zero problem
        :param kt: float, kt=2.5 when T=300K and unit is kJ/mol
        :param max:
        :return: list of floats
        """
        if x < minX:
            x = minX / 2.0
        return -1.0 * kbt * np.log(x / max)

    def pmf2d(self, X, Y, minX, kbt=2.5, bins=20):

        hist, edges_1, edges_2= np.histogram2d(X, Y, bins=bins)
        hist = hist / np.max(hist)

        hist = pd.DataFrame(hist)
        hist = hist.replace(0.0, minX).values
        pmf = -1.0 * kbt * np.log(hist / np.max(hist))

        return edges_1, edges_2, pmf

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

class LineFit :
    """
    Perform line fitting related calculations.

    Parameters
    ----------
    points: list, ndarray, or pd.Series
        the input coordinates

    Attributes
    ----------

    Methods
    -------

    """
    def __init__(self, points):
        self.points = np.asarray(points)

    def fit_line(self):
        """
        given a list of atoms (with xyz coordinate)
        return their normalized fitting line vector

        :return: list, vector
        """
        datamean = self.points.mean(axis=0)
        uu, dd, vv = np.linalg.svd(self.points - datamean)

        return vv[0] / np.sqrt(sum([x ** 2 for x in vv[0]]))
        #return vv[0]

    def unit_vector(self, vector):
        """ Returns the unit vector of the vector.  """
        return vector / np.linalg.norm(vector)

    def angle_between(self, v1, v2):
        '''Returns the angle in radians between vectors 'v1' and 'v2'::

                >>> angle_between((1, 0, 0), (0, 1, 0))
                1.5707963267948966
                >>> angle_between((1, 0, 0), (1, 0, 0))
                0.0
                >>> angle_between((1, 0, 0), (-1, 0, 0))
                3.141592653589793'''

        v1_u = self.unit_vector(v1)
        v2_u = self.unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

    def dotproduct(self, v1, v2):
        return np.dot(v1, v2)

    def length(self, v):
        return np.sqrt(self.dotproduct(v, v))

    def vector_angle(self, v1, v2):
        radian = np.arccos(self.dotproduct(v1, v2) / (self.length(v1) * self.length(v2)))
        degree = radian / np.pi * 180.0

        return degree