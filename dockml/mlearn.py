# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from .algorithms import BasicAlgorithm
from sklearn import preprocessing

class FeatureSelection :

    def __init__(self):
        pass

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

        #if firstNo > importance.shape[0] :
        #    firstNo = importance.shape[0]

        select = (importance > sorted(importance, reverse=True)[firstNo])
        dataset_sel = X.iloc[:, select]

        return (importance, dataset_sel)

    def mutualInformation(self, x, y):
        """
        calculate mutual information
        :param x: a feature vector
        :param y: a feature vector
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
        return X[key_fe], key_fe

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

        return X_trans, pca

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

    def bootstrap_resample(X, n=None):
        """ Bootstrap resample an array_like
        source: https://gist.github.com/aflaxman/6871948
        Parameters
        ----------
        X : array_like
          data to resample
        n : int, optional
          length of resampled array, equal to len(X) if n==None
        Results
        -------
        returns X_resamples
        """
        if isinstance(X, pd.Series):
            X = X.copy()
            X.index = range(len(X.index))
        if n == None:
            n = len(X)

        resample_i = np.floor(np.random.rand(n) * len(X)).astype(int)
        X_resample = np.array(X[resample_i])  # TODO: write a test demonstrating why array() is important
        return X_resample

class ModelEvaluation :

    def __init__(self, y, pred_y, threshold=0.5):
        '''

        :param y: original Y label, list (or vector)
        :param pred_y: predicted Y, list or (vector)
        '''

        self.Y      = list(y)
        self.Y_pred = list(pred_y)

        self.threshold = threshold

    def getMatrix(self):

        self.tp, self.tn, self.fp, self.fn = self.confusionMatrix()

    def scores(self):
        self.mccf = self.mcc()
        self.f1   = self.f1_score()
        self.fpr = self.fallout()

    def confusionMatrix(self):
        '''
        Confusion Matrix
        return True Positive, True Negative, False Positive, Negative

        parameter
        ------------------------------------------------
        :return: tuplex
        '''
        tp, tn, fp, fn = 0, 0, 0, 0

        for i in range(len(self.Y)):

            if self.Y[i] > self.threshold and self.Y_pred[i] > self.threshold :
                tp += 1.0
            elif self.Y[i] > self.threshold and self.Y_pred[i] <= self.threshold :
                fn += 1.0
            elif self.Y[i] <= self.threshold and self.Y_pred[i] > self.threshold :
                fp += 1.0
            elif self.Y[i] <= self.threshold and  self.Y_pred[i] <= self.threshold :
                tn += 1.0
            else :
                print("True Value: ", self.Y[i], " Predicted Value ", self.Y_pred[i])

        return tp, tn, fp, fn

    def mcc(self):
        '''
        Matthew’s correlation coefficient
        :return: float
        '''
        mccf = (self.tp * self.tn - self.fp * self.fn) / \
               np.sqrt((self.tp + self.fp) * (self.tp + self.fn) * (self.tn + self.fp) * (self.tn + self.fn))

        return mccf

    def sensitivity(self):
        '''
        sensitivity: True positive rate
        sensitivity  TPR = TP / P = TP / (TP+FN)
        :return: float
        '''

        return self.tp / (self.tp + self.fn)

    def specificity(self):
        '''
        specificity: True negative rate
        specificity (SPC) or true negative rate SPC=TN/T=TN/(TN+FP)
        :return: float
        '''

        return self.tn / (self.tn + self.fp)

    def accuracy(self):
        '''
        accuracy: ACC. the prediction accuracy Q
        accuracy (TP + TN) / (TP + TN + FP + FN)
        :return:
        '''
        return (self.tp + self.tn) / (self.tp + self.tn + self.fn + self.fp)

    def fallout(self):
        '''
        fall-out rate, the false positive rate
        FPR = 1 - SPC
        :return:
        '''

        return 1.0 - self.specificity()

    def f1_score(self):
        '''
        F1 score, the harmonic mean of precision and sensitivity
        F1 score 2*TP / (2*TP + FP + FN)
        :return:
        '''

        return 2.0 * self.tp / ( 2.0 * self.tp + self.fp + self.fn )

    def auc(self):
        '''
        calculate auc
        :return: float, auc score
        '''

        from sklearn.metrics import roc_curve, auc

        # find the TP-FP pairs
        fpr, tpr, thresholds = roc_curve(self.Y, self.Y_pred)

        # calculate Area Under Curve
        auc_score = auc(fpr, tpr)

        return auc_score

    def draw_roc_curve(self, pos_label=1,
                       c='black', figout="",
                       title="", dash_c="orange",
                       show=False,
                       ):
        '''
        draw roc curve
        :param pos_label:
        :param c:
        :param figout:
        :param label:
        :param title:
        :return:
        '''

        from sklearn.metrics import roc_curve, auc
        from matplotlib import pyplot as plt

        # find the TP-FP pairs
        fpr, tpr, thresholds = roc_curve(self.Y, self.Y_pred, pos_label=pos_label)

        # calculate Area Under Curve
        auc_score = auc(fpr, tpr)

        plt.plot(fpr, tpr, label="ROC curve (area = %.02f)"%auc_score, color=c, lw=2.0)
        plt.plot([0, 1], [0, 1], 'k--', lw=2.0, c=dash_c)

        plt.xlabel('False positive rate', fontsize=16)
        plt.ylabel('True positive rate', fontsize=16)

        if len(title) :
            plt.title('ROC curve')

        plt.legend(loc='best')

        # equalize the x-axis and y-axis
        plt.axis('scaled')

        if len(figout) :
            plt.savefig(figout)

        if show :
            plt.show()

        return 1

    def enrichment(self, score_list, percent=0.2):
        '''
        report the enrichment factor
        ef = hit_select / hit_total * ( N_total / N_select )
        :param score_list: list, a list containing unordered ligand:score information
                           [ [ score, is_active], ..., ]
        :param percent: float, ratio of ef
        :return:
        '''

        sorted_list = sorted(score_list, key=lambda x: x[0], reverse=True)

        self.N_total = float(len(sorted_list))
        self.N_select= float(self.N_total * percent)

        select_list = sorted_list[: int(self.N_select)]

        self.hit_select = float(sum([ x[-1] for x in select_list]))
        self.hit_total  = float(sum([ x[-1] for x in score_list ]))

        #print(N_total, N_select, hit_total, hit_select)

        return (self.hit_select / self.hit_total) * ( self.N_total / self.N_select)

class BindingFeatureClean:
    def __init__(self, x, y):
        self.X = x
        self.Y = y

    def processPipeLine(self):
        # remove all zero
        self.X = self.removeAllZeroFeatures(self.X)
        print("Dropping all-zeroes columns ... ")

        # normalized X
        self.X = preprocessing.scale(self.X)
        print("Normalized dataset ... ")
        self.X = pd.DataFrame(self.X)

        # data importance
        self.importance, self.X = \
            FeatureSelection().featureImportance(self.X, self.Y, firstNo=800)
        print("Calculating feature importances ... ")

        # remove highly correlated features
        self.X, self.key_features = FeatureSelection().removeCorrelated(self.X, 0.85)
        print("Remove highly correlated features ... ")

        # perform PCA projection
        self.X, self.pca = FeatureSelection().PCA(self.X)
        print("PCA transform completed ... ")

    def loadDataSet(self, fn='positive.csv'):
        '''
        load data set and return a dataframe
        :param fn:
        :return: pandas dataframe
        '''

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
        test_ = data[~msk]

        return train_, test_

class ModelTraining :
    def __init__(self, X, Y):

        self.train_X = None
        self.train_Y = None
        self.test_X  = None
        self.test_Y  = None

        self.cv_indices = []

        self.X = X
        self.Y = Y

    def cross_validation_split(self, X, k_fold=5):
        from sklearn.model_selection import KFold

        kf = KFold(n_splits=5)

        indices = []

        for train, test in kf.split(X) :
            indices.append([train, test])

        self.cv_indices = indices


    def planar_split(self, X, n=5):
        '''
        given a dataset, split it into n equal size random small sets
        :param X:
        :param n:
        :return:
        '''

        indexes = np.arange(X.shape[0])
        np.random.shuffle(indexes)

        chunk = len(indexes) / n
        groups = []

        for i in range(n-1) :
            groups.append(indexes[i*chunk: (i+1)*chunk ])
        groups.append(indexes[(n-1)*chunk : ])

        return groups
