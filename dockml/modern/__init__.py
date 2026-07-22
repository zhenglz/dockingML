"""
Modernized docking + ML pipeline.
==================================

Subpackages
-----------
docking    - wrappers around Vina, Vina 1.2+, Smina, QVina2, GNINA, Vina-CUDA
features   - modern protein-ligand interaction fingerprints
scoring    - ML/DL rescoring functions
models     - ML/DL model wrappers
pipeline   - end-to-end pipeline (prepare -> dock -> featurise -> rescore)
cli        - command-line entry points
"""

__version__ = "2.0.0"
