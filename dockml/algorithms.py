#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
basic math algorithms
"""

import math
import numpy as np
import pandas as pd


class BasicAlgorithm(object):

    def __init__(self):
        pass

    def switchFuction( self, x, d0=7.0, m=12, n=6):
        """Implement a rational switch function
        to enable a smooth transition

        Parameters
        ----------
        x : float
            the input distance value
        d0 : float
            the distance cutoff, it is usually 2 times of
            the real distance cutoff
        m : int
            the exponential index of higher order
        n : int
            the exponential index of lower order

        Returns
        -------
        switched_dist : float
            the switched continuous distance

        Notes
        -----
        the function is like:
          s= [1 - (x/d0)^6] / [1 - (x/d0)^12]
        d0 is a cutoff, should be twice larger than the distance cutoff

        """
        count = 0.0
        try:
            count = (1.0 - math.pow((x / d0), n)) / (1.0 - math.pow((x / d0), m))
        except ZeroDivisionError:
            print("Divide by zero, ", x, d0)

        return count

    def exponentialFunction(self, x, exp=2.0, k=1.0, c0=0.0, bias=0.0):
        """The exponential function

        Parameters
        ----------
        x : np.array
            the input variable X
        exp : float
            the power index
        k : float
            the constant
        c0 : float
            the shift of variable X
        bias : float
            the intercept

        Returns
        -------
        X_transformed : np.array

        """

        return k * (x + c0) ** exp + bias

    def pmf(self, x, minX=0.5, kbt=2.5, bins=20):
        """Calculate PMF 1D matrix given X

        Parameters
        ----------
        x : np.ndarray, shape = [ N, 1]
            The input data vector
        minX : float
            A minimium value for the probability calculation,
            to avoid log by zero error.
        kbt : float, default = 2.5
            The Kb * T value. If using KJ/mol as the unit,
            kbt is around 2.5
        bins : int, or array like, default = 20
            The number of bins, or pre-defined bin edges.

        Returns
        -------

        """
        dist, xedge = np.histogram(x, bins=bins)

        dist = dist / np.max(dist)
        dist[dist == 0.0] = minX / np.max(dist)

        return -1.0 * kbt * np.log(dist)

    def pmf2d(self, X, Y, minX, kbt=2.5, bins=20):
        """Calculate PMF 2D matrix given X and y

        Parameters
        ----------
        X : np.ndarray, or array like, shape = [N, 1]
            The first input vector
        Y : np.ndarray, or array like, shape = [N, 1]
            The 2nd input vector
        minX : float
            A minimium value for the probability calculation,
            to avoid log by zero error.
        kbt : float, default = 2.5
            The Kb * T value. If using KJ/mol as the unit,
            kbt is around 2.5
        bins : int, or array like, default = 20
            The number of bins, or pre-defined bin edges.

        Returns
        -------
        Xedge : np.ndarray
            The edges of X-axis
        Yedge : np.ndarray
            The edges of Y-axis
        pmf : np.ndarray, shape = [ nbins, nbins]
            The PMF matrix
        """

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


class PlaneFit(object):
    """

    Parameters
    ----------
    points: ndarray, M*3
        input points, in np.ndarray shape

    Methods
    -------

    """

    def __init__(self, points):
        self.points = points

    def fitPlane(self):
        """
        fit some points to a plane
        https://math.stackexchange.com/questions/99299/best-fitting-plane-given-a-set-of-points
        plane function: ax + by - z + c = 0
        or: ax + by - z = -c
        Ax + By + Cz + D = 0
        we need to determine [a, b, c]
        A B C D = a  b -1 c

        Parameters
        ----------

        Returns
        -------
        fit: array, [a, b, c]
        """

        # prepare dataset
        xs = np.array(self.points)[:, 0]
        ys = np.array(self.points)[:, 1]
        zs = np.array(self.points)[:, 2]

        tmp_A = []
        tmp_B = []
        for i in range(xs.shape[0]):
            tmp_A.append([xs[i], ys[i], 1])
            tmp_B.append(zs[i])

        B = np.matrix(tmp_B).T
        A = np.matrix(tmp_A)

        # do fit, fit points to a plane
        fit = (A.T * A).I * A.T * B
        errors = B - A * fit
        residual = np.linalg.norm(errors)

        return fit, residual

    def point_distance(self, params, point):
        """
        https://mathinsight.org/distance_point_plane
        determine the distance between a point to a plane
        the function of a function Ax + By + Cz + D = 0
        P (x0, y0, z0 ) to the plane
        distance (x0*A + y0*B + z0*C + D)/sqrt(A^2 + B^2 + C^2)

        Parameters
        ----------
        params: np.array,
            a list of parameters (A, B, D), while C=-1
        point: list, or np.array,
            x y z coordinates of a point

        Returns
        -------
        distance: float,
            the distance of a point to a defined plane
        """

        distance = math.sqrt(params[0] ** 2 + params[1] ** 2 + 1)

        distance = (params[0] * point[0] +
                    params[1] * point[1] +
                    (-1.0) * point[2] +
                    params[2]
                    ) / distance

        return distance


class LineFit(object):
    """Perform line fitting related calculations.

    Parameters
    ----------
    points: list, np.ndarray, or pd.Series, shape = [ N, M]
        the input coordinates, N is number of points
        and M is the number of dimensions. For euclidean axis,
        M = 3.

    Attributes
    ----------
    points : np.ndarray
        the input coordinates

    Methods
    -------

    """
    def __init__(self, points):
        self.points = np.asarray(points)

    def fit_line(self):
        """Given a list of atoms (with xyz coordinate)
        return their normalized fitting line vector

        Return
        line_vector: list, vector
        """
        datamean = self.points.mean(axis=0)
        uu, dd, vv = np.linalg.svd(self.points - datamean)

        return vv[0] / np.sqrt(sum([x ** 2 for x in vv[0]]))

    def unit_vector(self, vector):
        """ Returns the unit vector of the vector.

        Parameters
        ----------
        vector: np.array,
            a vector,

        Returns
        -------
        unit_vector: np.array,
            the result unit vector

        """
        return vector / np.linalg.norm(vector)

    def angle_between(self, v1, v2):
        """Returns the angle in radians between vectors 'v1' and 'v2'::

        >>> angle_between((1, 0, 0), (0, 1, 0))
        1.5707963267948966
        >>> angle_between((1, 0, 0), (1, 0, 0))
        0.0
        >>> angle_between((1, 0, 0), (-1, 0, 0))
        3.141592653589793

        Parameters
        ----------
        v1: list,
        v2: list,

        Returns
        -------
        angle: float,
            the cross angle of two lines, in radian unit

        """

        v1_u = self.unit_vector(v1)
        v2_u = self.unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
