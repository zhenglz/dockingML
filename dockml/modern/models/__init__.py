"""Optional deep-learning rescoring models.

These require PyTorch (installed via ``pip install dockingML[dl]``) and
provide two reference architectures:

    * :class:`Simple3DCNN`    - a KDEEP/DeepAtomicCharge-style 3D CNN that
      operates on voxelised atom-type grids of the binding pocket.
    * :class:`GNNRescorer`    - a DGL/PyG graph neural network that uses
      inter-atomic contact graphs of the pocket-ligand complex.

Both classes expose a ``sklearn``-like ``fit``/``predict_proba`` interface.
For production models, prefer dedicated libraries (gnina, PIGNet2, OnionNet-2).
"""
from .cnn import Simple3DCNN  # noqa: F401
from .gnn import GNNRescorer  # noqa: F401
