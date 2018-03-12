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
        :param X: panda dataframe
        :param Y: panda dataframe
        :param firstNo:
        :return: (dataset importance, selected dataset)
        """
        from sklearn.ensemble import ExtraTreesClassifier

        model = ExtraTreesClassifier()
        model.fit(X, Y)

        importance = np.array(model.feature_importances_)

        select = (importance > sorted(importance, reverse=True)[firstNo])
        dataset_sel = X.iloc[:, select]

        return (importance, dataset_sel)

    def mutualInformation(self, x, y):
        """
        calculate mutual information
        :param x:
        :param y:
        :return:
        """

        from dockml import algorithms

        algo = algorithms.BasicAlgorithm()

        HX = algo.entropy1D(x)
        HY = algo.entropy1D(y)
        HXY= algo.entropy2D(x, y)

        return HX + HY - HXY

    def correlations(self, X):
        '''
        calculate pair-wise correlations between features
        return their correlation matrix
        :param X: pandas dataframe
        :return: N*N matrix, correlation matrix
        '''

        correlations = []
        features = X.columns

        for f1 in features:
            for f2 in features:
                corr = X[f1].corr(X[f2])
                correlations.append(corr)

        corr = np.reshape(correlations, (len(features), len(features)))

        return corr

    def removeCorrelated(self, X, corr_cutoff = 0.85):
        """
        remove the highly correlated features,
        and only keep the more important uncorrelated features
        :param corr: N*N correlation matrix
        :param X: pandas dataframe
        :param corr_cutoff: float, the cutoff value for correlation determination
        :return: pandas dataframe with only uncorrelated features
        """
        import math

        # get feature correlation matrix
        corr = self.correlations(X)

        features = X.columns
        features = list(features)
        key_fe = list(X.columns.values)

        for i in range(len(features)):
            for j in range(len(features)):
                if math.fabs(corr[i, j]) > corr_cutoff and i != j and \
                                features[i] in key_fe and \
                                features[j] in key_fe:
                        key_fe.remove(features[i])
        return X[key_fe]

    def PCA(self, X):
        """
        perform PCA analysis on dataset X
        :param X: ndarray or pandas dataframe
        :return: transformed X (same dimensionality) and pca object for transforming other dataset
        """

        from sklearn import decomposition

        pca = decomposition.PCA()
        # preform PCA analysis
        pca.fit(X)
        # predict data set
        X_trans = pca.transform(X)

        return (X_trans, pca)