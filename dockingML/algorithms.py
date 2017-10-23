#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
basic math algorithms
"""

import math

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