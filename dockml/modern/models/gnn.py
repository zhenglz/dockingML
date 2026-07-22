"""
Reference Graph Neural Network (GNN) rescorer.
==============================================

This is a minimal GNN built on top of ``torch_geometric``. Nodes are ligand +
binding-pocket heavy atoms (within e.g. 6 Å of the ligand). Edges are added
between pairs within 4.5 Å. Each node carries a small feature vector
(element one-hot, residue-type one-hot, partial charge, etc.).

For real work prefer GraphDelta, PIGNet, or POTENTIALNET; this module
provides an easy to hack starting point.
"""
from __future__ import annotations

import logging
from typing import List, Tuple

import numpy as np

log = logging.getLogger("dockml.modern.models.gnn")

try:
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    from torch_geometric.data import Data, DataLoader
    from torch_geometric.nn import GCNConv, global_mean_pool
    _HAVE_PYG = True
except Exception:
    _HAVE_PYG = False


ELEMENTS = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "Du"]
RESIDUES = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
            "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
            "TYR", "VAL", "LIG"]


def complex_to_graph(rec_xyz, rec_ele, rec_resname,
                     lig_xyz, lig_ele,
                     cutoff: float = 4.5):
    """Build a ``torch_geometric.Data`` graph from a protein-ligand complex.

    Node features: one-hot element + one-hot residue (LIG for ligand atoms)
    + is-ligand flag.
    """
    if not _HAVE_PYG:
        raise ImportError("Install torch-geometric: pip install torch-geometric")
    nr, nl = len(rec_ele), len(lig_ele)
    xyz = np.vstack([rec_xyz, lig_xyz])
    ele = list(rec_ele) + list(lig_ele)
    res = list(rec_resname) + ["LIG"] * nl
    is_lig = np.array([0] * nr + [1] * nl, dtype=np.float32)

    def _onehot(v, vocab):
        v = [0] * len(vocab)
        return v  # filled below
    feat = np.zeros((len(xyz), len(ELEMENTS) + len(RESIDUES) + 1), dtype=np.float32)
    for i, (e, r) in enumerate(zip(ele, res)):
        try:
            feat[i, ELEMENTS.index(e)] = 1.0
        except ValueError:
            feat[i, ELEMENTS.index("Du")] = 1.0
        try:
            feat[i, len(ELEMENTS) + RESIDUES.index(r)] = 1.0
        except ValueError:
            pass
        feat[i, -1] = is_lig[i]

    # edges based on distance cutoff (k-nearest for ligand to avoid huge graphs)
    d2 = np.sum((xyz[:, None, :] - xyz[None, :, :]) ** 2, axis=2)
    src, dst = np.where((d2 < cutoff ** 2) & (d2 > 1e-6))
    edge_index = torch.as_tensor(np.vstack([src, dst]), dtype=torch.long)
    x = torch.as_tensor(feat, dtype=torch.float32)
    pos = torch.as_tensor(xyz, dtype=torch.float32)
    batch = torch.zeros(x.shape[0], dtype=torch.long)
    return Data(x=x, edge_index=edge_index, pos=pos, batch=batch)


class GNNRescorer:
    """Simple GCN classifier for protein-ligand graphs."""

    def __init__(self, hidden: int = 64, lr: float = 1e-3, device: str = "cpu"):
        if not _HAVE_PYG:
            raise ImportError("Install torch-geometric to use the GNN rescorer")
        self.hidden = hidden
        self.lr = lr
        self.device = device
        self._build()

    def _build(self):
        in_dim = len(ELEMENTS) + len(RESIDUES) + 1

        class Net(nn.Module):
            def __init__(self):
                super().__init__()
                self.conv1 = GCNConv(in_dim, self.hidden)
                self.conv2 = GCNConv(self.hidden, self.hidden)
                self.conv3 = GCNConv(self.hidden, self.hidden)
                self.lin = nn.Linear(self.hidden, 2)
            def forward(self, x, edge_index, batch):
                x = F.relu(self.conv1(x, edge_index))
                x = F.relu(self.conv2(x, edge_index))
                x = F.relu(self.conv3(x, edge_index))
                x = global_mean_pool(x, batch)
                return self.lin(x)
        self.net = Net().to(self.device)
        self.opt = torch.optim.Adam(self.net.parameters(), lr=self.lr)
        self.loss = nn.CrossEntropyLoss()

    def fit(self, graphs: List, labels: np.ndarray, epochs: int = 30,
            batch: int = 16):
        loader = DataLoader(list(zip(graphs, labels.astype(int))),
                            batch_size=batch, shuffle=True)
        for ep in range(epochs):
            self.net.train()
            total = 0
            for data, y in loader:
                data = data.to(self.device); y = y.to(self.device)
                self.opt.zero_grad()
                out = self.net(data.x, data.edge_index, data.batch)
                loss = self.loss(out, y)
                loss.backward(); self.opt.step()
                total += loss.item()
            log.info("epoch %d  loss=%.3f", ep+1, total / max(len(loader),1))

    def predict_proba(self, graphs: List) -> np.ndarray:
        self.net.eval()
        loader = DataLoader(graphs, batch_size=32, shuffle=False)
        ps = []
        with torch.no_grad():
            for data in loader:
                data = data.to(self.device)
                p = torch.softmax(self.net(data.x, data.edge_index, data.batch), 1)[:, 1]
                ps.append(p.cpu().numpy())
        return np.concatenate(ps)
