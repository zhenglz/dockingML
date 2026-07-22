"""
Reference 3D-CNN rescoring model for voxelised protein-ligand grids.
==================================================================

This is a minimal reimplementation of the idea behind KDEEP/DeepAtomicCharge
(Jimenez et al., J Cheminf 2018; Stepniewska-Dziubinska et al., Bioinformatics
2018): represent the binding site as a multi-channel voxel grid (each channel
= one atom type / property) and train a small 3D CNN to predict pose
probability or pKd.

Intended for learning / small-scale experiments; for heavy production use
gnina/DeepChem/OpenDrugDiscovery toolkits directly.
"""
from __future__ import annotations

import logging
from typing import Tuple, Optional

import numpy as np

log = logging.getLogger("dockml.modern.models.cnn")


def voxelise(lig_xyz: np.ndarray, rec_xyz: np.ndarray,
             lig_types: np.ndarray, rec_types: np.ndarray,
             grid_size: int = 24, resolution: float = 1.0,
             n_channels: int = 8) -> np.ndarray:
    """Create a (C, G, G, G) voxel grid centred at the ligand centroid.

    Channels (default 8) correspond to:
        0 receptor hydrophobic (C), 1 rec aromatic, 2 rec HBD, 3 rec HBA/negative,
        4 rec positive, 5 ligand C/halogen, 6 ligand polar (N/O/S), 7 ligand metal/other.

    This is a coarse Gaussian-density voxelisation; for full-feature grids
    use ``oddt`` or ``gnina``'s grid maker.
    """
    centre = lig_xyz.mean(axis=0)
    half = grid_size * resolution / 2.0
    edges = np.linspace(-half, half, grid_size + 1)
    gridc = 0.5 * (edges[:-1] + edges[1:]) + centre[:, None, None, None]
    # gridc shape (3, G, G, G)
    G = grid_size
    grid = np.zeros((n_channels, G, G, G), dtype=np.float32)
    sigma = resolution * 0.8

    def _fill(xyz, types, ch_offset):
        d = np.sqrt(((gridc[0] - xyz[:, 0].reshape(-1, 1, 1, 1)) ** 2 +
                     (gridc[1] - xyz[:, 1].reshape(-1, 1, 1, 1)) ** 2 +
                     (gridc[2] - xyz[:, 2].reshape(-1, 1, 1, 1)) ** 2))
        dens = np.exp(-d**2 / (2 * sigma**2))  # (N, G, G, G)
        for ch in range(4):
            m = (types == ch + ch_offset)
            if m.any():
                grid[ch + ch_offset] += dens[m].sum(0)

    _fill(rec_xyz, rec_types, 0)
    _fill(lig_xyz, lig_types, 4)
    return grid


class Simple3DCNN:
    """Minimal 3D-CNN classifier for rescoring voxelised poses.

    PyTorch is imported lazily so the rest of the pipeline works without it.
    """

    def __init__(self, channels: int = 8, grid_size: int = 24, lr: float = 1e-3,
                 device: str = "cpu"):
        self.channels = channels
        self.grid_size = grid_size
        self.lr = lr
        self.device = device
        self.net = None
        self._built = False

    def _build(self):
        try:
            import torch
            import torch.nn as nn
        except ImportError as e:
            raise ImportError("Install torch: pip install torch") from e

        class Net(nn.Module):
            def __init__(self, c, g):
                super().__init__()
                self.features = nn.Sequential(
                    nn.Conv3d(c, 32, 3, padding=1), nn.ReLU(), nn.MaxPool3d(2),
                    nn.Conv3d(32, 64, 3, padding=1), nn.ReLU(), nn.MaxPool3d(2),
                    nn.Conv3d(64, 128, 3, padding=1), nn.ReLU(),
                    nn.AdaptiveAvgPool3d(1),
                )
                self.head = nn.Sequential(
                    nn.Flatten(), nn.Linear(128, 64), nn.ReLU(),
                    nn.Dropout(0.3), nn.Linear(64, 2),
                )
            def forward(self, x):
                return self.head(self.features(x))

        self.net = Net(self.channels, self.grid_size).to(self.device)
        self.opt = __import__("torch").optim.Adam(self.net.parameters(), lr=self.lr)
        self.loss = __import__("torch.nn").CrossEntropyLoss()
        self._built = True

    def fit(self, X: np.ndarray, y: np.ndarray, epochs: int = 10, batch: int = 32,
            val_split: float = 0.1):
        if not self._built:
            self._build()
        import torch
        from torch.utils.data import DataLoader, TensorDataset
        X = torch.as_tensor(X, dtype=torch.float32)
        y = torch.as_tensor(y, dtype=torch.long)
        n = X.shape[0]
        nv = int(n * val_split)
        idx = torch.randperm(n)
        trX, trY = X[idx[nv:]], y[idx[nv:]]
        vaX, vaY = X[idx[:nv]], y[idx[:nv]]
        loader = DataLoader(TensorDataset(trX, trY), batch_size=batch, shuffle=True)
        for ep in range(epochs):
            self.net.train()
            for xb, yb in loader:
                xb = xb.to(self.device); yb = yb.to(self.device)
                self.opt.zero_grad()
                out = self.net(xb)
                loss = self.loss(out, yb)
                loss.backward(); self.opt.step()
            if nv:
                self.net.eval()
                with torch.no_grad():
                    vout = self.net(vaX.to(self.device))
                    vacc = (vout.argmax(1) == vaY.to(self.device)).float().mean().item()
                log.info("epoch %d  train_loss=%.3f  val_acc=%.3f", ep+1, loss.item(), vacc)

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        import torch
        self.net.eval()
        X = torch.as_tensor(X, dtype=torch.float32).to(self.device)
        with torch.no_grad():
            p = torch.softmax(self.net(X), 1)[:, 1].cpu().numpy()
        return p

    def save(self, path: str):
        import torch
        torch.save({"state": self.net.state_dict(),
                    "channels": self.channels, "grid_size": self.grid_size}, path)

    @classmethod
    def load(cls, path: str, device: str = "cpu") -> "Simple3DCNN":
        import torch
        d = torch.load(path, map_location=device)
        obj = cls(channels=d["channels"], grid_size=d["grid_size"], device=device)
        obj._build()
        obj.net.load_state_dict(d["state"])
        return obj
