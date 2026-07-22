#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Classical ML rescoring models
==============================

Provides wrappers for the most successful tree-ensemble models used in
docking rescoring (XGBoost, LightGBM, Random Forest) plus simple
logistic regression / linear baselines. Models accept fixed-length
feature vectors from ``dockml.modern.features`` and predict
either pose-quality (RMSD classification: native-like vs. decoy) or
binding affinity (regression).

For deep-learning models (3D CNN, GNN), see ``dockml.modern.models``.
"""
from __future__ import annotations

import logging
import pickle
from pathlib import Path
from typing import Dict, List, Optional, Union, Tuple

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.model_selection import cross_validate, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, average_precision_score, mean_squared_error

log = logging.getLogger("dockml.modern.scoring.classical")


class ClassicalRescorer:
    """Unified wrapper for classical ML rescoring.

    Parameters
    ----------
    model : str
        One of ``{"rf", "xgb", "lgbm", "lr", "ridge", "gbt"}``.
    task : str
        ``"classify"`` (active/decoy or good/bad pose) or ``"regress"`` (pKd/ΔG).
    scale : bool
        Standardize features (recommended for lr/ridge).
    model_kwargs : dict
        Extra keyword arguments passed to the underlying constructor.
    """

    def __init__(self, model: str = "xgb", task: str = "classify",
                 scale: bool = False, **model_kwargs):
        self.model_name = model
        self.task = task
        self.scaler = StandardScaler() if scale else None
        self.model = self._build(model, task, **model_kwargs)

    # ------------------------------------------------------------------
    @staticmethod
    def _build(name: str, task: str, **kw):
        if task == "classify":
            if name == "rf":
                return RandomForestClassifier(
                    n_estimators=kw.pop("n_estimators", 500),
                    max_depth=kw.pop("max_depth", None),
                    n_jobs=kw.pop("n_jobs", -1),
                    random_state=42, class_weight="balanced", **kw)
            if name == "gbt":
                return GradientBoostingClassifier(
                    n_estimators=kw.pop("n_estimators", 300),
                    learning_rate=kw.pop("learning_rate", 0.05),
                    max_depth=kw.pop("max_depth", 4), random_state=42, **kw)
            if name == "xgb":
                try:
                    from xgboost import XGBClassifier
                except ImportError as e:
                    raise ImportError("Install xgboost: pip install xgboost") from e
                return XGBClassifier(
                    n_estimators=kw.pop("n_estimators", 500),
                    learning_rate=kw.pop("learning_rate", 0.05),
                    max_depth=kw.pop("max_depth", 6),
                    use_label_encoder=False, eval_metric="logloss",
                    tree_method=kw.pop("tree_method", "hist"),
                    n_jobs=kw.pop("n_jobs", -1), random_state=42, **kw)
            if name == "lgbm":
                try:
                    from lightgbm import LGBMClassifier
                except ImportError as e:
                    raise ImportError("Install lightgbm: pip install lightgbm") from e
                return LGBMClassifier(
                    n_estimators=kw.pop("n_estimators", 500),
                    learning_rate=kw.pop("learning_rate", 0.05),
                    num_leaves=kw.pop("num_leaves", 31),
                    n_jobs=kw.pop("n_jobs", -1), random_state=42, **kw)
            if name == "lr":
                return LogisticRegression(max_iter=2000, class_weight="balanced", **kw)
            raise ValueError(f"Unknown classifier {name}")
        else:  # regress
            if name == "rf":
                return RandomForestRegressor(n_estimators=500, n_jobs=-1, random_state=42, **kw)
            if name == "ridge":
                return Ridge(alpha=1.0, **kw)
            if name == "xgb":
                from xgboost import XGBRegressor
                return XGBRegressor(n_estimators=500, learning_rate=0.05,
                                    max_depth=6, tree_method="hist",
                                    n_jobs=-1, random_state=42, **kw)
            if name == "lgbm":
                from lightgbm import LGBMRegressor
                return LGBMRegressor(n_estimators=500, learning_rate=0.05,
                                     num_leaves=31, n_jobs=-1, random_state=42, **kw)
            raise ValueError(f"Unknown regressor {name}")

    # ------------------------------------------------------------------
    def fit(self, X: Union[pd.DataFrame, np.ndarray],
            y: Union[pd.Series, np.ndarray]) -> "ClassicalRescorer":
        X = np.asarray(X, dtype=float)
        if self.scaler:
            X = self.scaler.fit_transform(X)
        self.model.fit(X, np.asarray(y))
        return self

    def predict_proba(self, X) -> np.ndarray:
        X = np.asarray(X, dtype=float)
        if self.scaler:
            X = self.scaler.transform(X)
        if hasattr(self.model, "predict_proba"):
            return self.model.predict_proba(X)[:, 1]
        return self.model.predict(X)

    def predict(self, X) -> np.ndarray:
        X = np.asarray(X, dtype=float)
        if self.scaler:
            X = self.scaler.transform(X)
        return self.model.predict(X)

    # ------------------------------------------------------------------
    def cv_evaluate(self, X, y, n_fold: int = 5,
                    metrics=("roc_auc", "average_precision")) -> Dict[str, float]:
        X = np.asarray(X, dtype=float)
        y = np.asarray(y)
        if self.scaler:
            X = StandardScaler().fit_transform(X)
        scoring = list(metrics) if self.task == "classify" else ["neg_root_mean_squared_error", "r2"]
        skf = StratifiedKFold(n_splits=n_fold, shuffle=True, random_state=42) \
            if self.task == "classify" else n_fold
        res = cross_validate(self.model, X, y, cv=skf, scoring=scoring, n_jobs=1)
        return {k: float(np.mean(v)) for k, v in res.items() if k.startswith("test_")}

    # ------------------------------------------------------------------
    def save(self, path: str):
        with open(path, "wb") as f:
            pickle.dump({"model": self.model, "scaler": self.scaler,
                         "model_name": self.model_name, "task": self.task}, f)

    @classmethod
    def load(cls, path: str) -> "ClassicalRescorer":
        import pickle as _pk
        with open(path, "rb") as f:
            d = _pk.load(f)
        obj = cls.__new__(cls)
        obj.model = d["model"]; obj.scaler = d["scaler"]
        obj.model_name = d["model_name"]; obj.task = d["task"]
        return obj


class ConsensusRescorer:
    """Simple arithmetic / rank-based consensus of multiple models/scores.

    Consensus scoring consistently improves VS hit rates over any single SF.
    """
    def __init__(self, models: Optional[List] = None,
                 weights: Optional[List[float]] = None):
        self.models = models or []
        self.weights = np.asarray(weights or [1.0] * len(self.models))
        self.weights = self.weights / self.weights.sum()

    def add(self, model, weight: float = 1.0):
        self.models.append(model)
        self.weights = np.append(self.weights, weight)
        self.weights = self.weights / self.weights.sum()

    def predict(self, X) -> np.ndarray:
        scores = []
        for m in self.models:
            if hasattr(m, "predict_proba"):
                s = m.predict_proba(X)
            else:
                s = m.predict(X)
            scores.append(np.asarray(s))
        scores = np.vstack(scores).T  # (N, K)
        # z-score normalization per model, then weighted mean
        z = (scores - scores.mean(axis=0)) / (scores.std(axis=0) + 1e-9)
        return z @ self.weights

    def rank_average(self, X, ascending: bool = False) -> np.ndarray:
        scores = []
        for m in self.models:
            if hasattr(m, "predict_proba"):
                s = m.predict_proba(X)
            else:
                s = m.predict(X)
            order = np.argsort(np.argsort(s))  # ranks
            if ascending:
                order = len(s) - order
            scores.append(order)
        return np.mean(scores, axis=0)
