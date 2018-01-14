# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from .algorithms import BasicAlgorithm

class DataClean :

    def __init__(self):
        pass

class FeatureSelection :

    def __init__(self):
        pass

    def informationGain(self, plus_data, negative_data):
        '''
        calculate information gain (mutual information)
        :param plus_data: array, a feature from positive dataset
        :param negative_data: array, a feature from negative dataset
        :return:
        '''
        algo = BasicAlgorithm()

        mi = algo.entropy1D(plus_data) + algo.entropy1D(negative_data) - \
             algo.entropy1D(list(plus_data) + list(negative_data))


        return mi