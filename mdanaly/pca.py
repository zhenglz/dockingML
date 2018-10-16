import numpy as np
import pandas as np
import sklearn

class PCA(object):

    def __init__(self, n_components=2):
        self.n_components = n_components

        # attributes
        self.trained_ = False
        self.X_transformed_ = None

        self.scaled_ = False
        self.scaler_ = None
        self.X_scaled = None

        self.pca_obj = None
        self.eigvalues_ = None
        self.eigvectors_ = None

    def fit(self, X):
        """
        fit a pca object

        Parameters
        ----------
        X: numpy ndarray, shape = [N, M]
            the input data matrix, N is the number of samples
            M is the number of dimensions

        Returns
        -------
        pca_obj: sklearn.decomposition.PCA object
            the pca object from sklearn

        """

        if self.scaled_:
            Xs = X
        else:
            self.scaler_ = sklearn.preprocessing.StandardScaler()
            Xs = self.scaler_.fit_transform(X)

        self.X_scaled = Xs

        pca_obj = sklearn.decomposition.PCA(n_components=self.n_components)
        pca_obj.fit(self.X_scaled)

        self.trained_ = True
        self.X_transformed_ = pca_obj.transform(X)

        self.pca_obj = pca_obj

        return pca_obj

    def transform(self, X):

        if self.trained_:
            return self.pca_obj.transform(X)
        else:
            print("Your pca object is not trained yet. Training it now ...")
            pca_obj = self.fit(X)

            return pca_obj.transform(X)

    def eigvalues(self):

        egival_ = self.pca_obj.explained_ratio

        self.eigvalues_ = egival_

        return egival_