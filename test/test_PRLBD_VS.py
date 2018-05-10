#!/usr/bin/env python
import pandas as pd
import numpy as np
from dockml import mlearn

class BindingFeatureClean :

    def __init__(self, x, y):

        self.X = x
        self.Y = y

        # remove all zero
        self.X = self.removeAllZeroFeatures(self.X)

        # normalized X
        self.X = mlearn.DataClean().normalization(self.X)

        # data importance
        self.importance, self.X = \
            mlearn.FeatureSelection().featureImportance(self.X, self.Y, firstNo=int(self.X.shape[0] / 2))

        # remove highly correlated features
        self.X, self.key_features = mlearn.FeatureSelection().removeCorrelated(self.X, 0.85)

        # perform PCA projection
        self.X, self.pca = mlearn.FeatureSelection().PCA(self.X)

    def loadDataSet(self, fn='positive.csv'):

        return pd.read_csv(fn, header=0, sep=',')

    def combineDataSet(self, df_1, df_2):
        '''
        combine two dataframe into one dataset
        :param df_1:
        :param df_2:
        :return:
        '''
        return pd.concat([df_1, df_2])

    def removeAllZeroFeatures(self, data):

        return data.dropna(axis=0, how='all')

    def splitDataSet(self, data, ratio=0.8):
        '''
        split a dataset based on a ratio
        return two dataset
        :param data: panda dataframe
        :param ratio: float
        :return:
        '''

        msk = np.random.rand(data.shape[0]) < ratio
        train_ = data[msk]
        test_  = data[~msk]

        return train_, test_

print(BindingFeatureClean('', ''))
