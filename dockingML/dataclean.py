# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd

class DataClean :

    def __init__(self):
        pass

    def loadDataFile(self, dataf, delimiter=",", header=0):
        """
        load data file to pd dataframe
        :param dataf: str, file name
        :param delimiter: str, seperation of columns
        :param header: int, index of header line
        :return: pd dataframe
        """
        data = pd.read_csv(dataf, delimiter=delimiter, header=header)

        return data

    def normalization(self, data):
        """
        normalize the data set
        :param data: ndarray, not including the "class" label
        :return: ndarray
        """

        X = data.apply(lambda x: (x - np.mean(x)) / (np.max(x) - np.min(x)))

        return X

    def removeAllZeroes(self, data):
        """
        remove the columns with all zeroes
        :param data: ndarray
        :return: ndarray
        """

        return data.dropna(axis=1, how="all")

    def featureImportance(self, X, Y, firstNo=100):
        """
        get the importance of each feature, and get the first several important features
        :param X:
        :param Y:
        :param firstNo:
        :return:
        """
        from sklearn.ensemble import ExtraTreesClassifier

        model = ExtraTreesClassifier()
        model.fit(X, Y)

        importance = np.array(model.feature_importances_)

        select = (importance > sorted(importance, reverse=True)[firstNo])
        dataset_sel = X.iloc[:, select]

        return importance, dataset_sel

    def mutualInformation(self, x, y):
        """
        calculate mutual information
        :param x:
        :param y:
        :return:
        """

        from dockingML import algorithms

        algo = algorithms.BasicAlgorithm()

        HX = algo.entropy1D(x)
        HY = algo.entropy1D(y)
        HXY= algo.entropy2D(x, y)

        return HX + HY - HXY
