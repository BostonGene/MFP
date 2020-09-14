import copy

import pandas as pd
from sklearn.neighbors import KNeighborsClassifier

from portraits.utils import median_scale


class KNeighborsClusterClassifier:
    def __init__(self, norm=True, algorithm='auto', clip=2, scale=False, k=35):
        """
        Classification using KNN. Fit with signature matrix and labels. Predict on the signatures.
        To use training cohort parameters for scaling set norm=True (If the cohorts are from the same batch)
        Data sets from different batches should be scaled
        :param norm:
        :param algorithm:
        :param clip:
        :param scale:
        :param k:
        """
        self.norm = norm
        self.median = 0
        self.mad = 1
        self.X = None
        self.y = None
        self.algorithm = algorithm
        self.model = None
        self.clip = clip
        self.scale = scale
        self.k = k

    def check_is_fitted(self):
        return (self.X is not None) and (self.y is not None) and (self.model is not None)

    def preprocess_data(self, X):
        x = copy.deepcopy(X)
        if self.scale:
            x = median_scale(x)

        x = (x - self.median) / self.mad
        if self.clip is not None and self.clip > 0:
            x = x.clip(-1 * self.clip, self.clip)
        return x

    def check_columns(self, X):
        if hasattr(self.X, 'columns'):
            try:
                return X[self.X.columns]
            except KeyError:
                raise Exception('Columns do not match')
        return X

    def fit(self, X, y):
        """
        :param X: pd.DataFrame, RNA data, columns - features, index - samples
        :param y: pd.Series, cluster labels
        """
        if X.shape[0] != len(y):
            raise Exception('Shapes do not match')

        if self.norm:
            self.median = X.median()
            self.mad = X.mad()

        self.X = self.preprocess_data(X)
        self.y = copy.deepcopy(y)

        self.model = KNeighborsClassifier(algorithm=self.algorithm, n_neighbors=self.k).fit(self.X, self.y)

        return self

    def predict(self, X):
        """
        Predict - return a pd.Series with the predicted cluster labels

        :param X: pd.DataFrame, RNA data, columns - features, index - samples
        :return: pd.Series, predicted cluster labels
        """
        if X.shape[1] != self.X.shape[1]:
            raise Exception('Shapes do not match')

        x_scaled = self.preprocess_data(self.check_columns(X))
        # Here self.model.predict is used in order to mimic its' way to select the class in case of equal probabilities
        return pd.Series(self.model.predict(x_scaled), index=x_scaled.index)

    def predict_proba(self, X):
        """
        :param X: pd.DataFrame, RNA data, columns - features, index - samples
        :return: pd.DataFrame, probabilities for each cluster. Index - samples, columns - clusters
        """

        if X.shape[1] != self.X.shape[1]:
            raise Exception('Shapes do not match')

        x_scaled = self.preprocess_data(self.check_columns(X))

        return pd.DataFrame(self.model.predict_proba(x_scaled).astype(float), index=x_scaled.index,
                            columns=self.model.classes_)
