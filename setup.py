from setuptools import setup, find_packages

"""
Modernized dockingML package.

Supports:
  - Modern docking engines: AutoDock Vina 1.2+, GNINA, QVina2, Smina, Vina-CUDA
  - Modern MD: GROMACS 2024+ recommended force fields (ff19SB/OPC, CHARMM36m/TIP3P-CHARMM, OpenFF 2.x)
  - Topology building: ParmEd, OpenFF Toolkit 2.x (besides legacy ACPYPE/AmberTools)
  - Features: vectorised residue-atom contacts, OnionNet-style multi-shell features,
             hydrogen-bond, hydrophobic, pi-pi, salt-bridge fingerprints
  - ML: scikit-learn, XGBoost, LightGBM, optional PyTorch/TorchDrug/Geometric deep learning rescoring
"""


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='dockingML',
    version='2.0.0',
    long_description=readme(),
    long_description_content_type='text/markdown',
    description='Modernized Molecular Docking, MD and Machine Learning pipeline for '
                'Structure-Based Drug Discovery',
    url='https://github.com/zhenglz/dockingML',
    author='zhenglz & contributors',
    author_email='zhenglz@outlook.com',
    license='GPL-3.0',
    packages=find_packages(exclude=['test', 'tests', '*.test', '*.tests']),
    include_package_data=True,
    package_data={
        '': ['data/*.dat', 'data/*.mdp', 'data/*.lib', 'data/*.gro',
             'automd/data/*.mdp', 'automd/data/*.gro', 'automd/data/*.lib'],
    },
    install_requires=[
        # numerical / data
        'numpy>=1.22',
        'scipy>=1.8',
        'pandas>=1.4',
        # ML classical
        'scikit-learn>=1.1',
        'xgboost>=1.7',
        'lightgbm>=3.3',
        # plotting
        'matplotlib>=3.5',
        'seaborn>=0.12',
        # network / community analysis
        'networkx>=2.8',
        # structure / trajectory handling
        'mdtraj>=1.9.7',
        'MDAnalysis>=2.3',
        'biotite>=0.36',
        'parmed>=3.4',
        # cheminformatics
        'rdkit-pypi; python_version < "3.12"',
        'openbabel-wheel; python_version < "3.12"',
        # parallel
        'joblib>=1.2',
        'tqdm>=4.65',
        # config
        'pyyaml>=6.0',
    ],
    extras_require={
        'dl': [
            'torch>=2.0',
            'torch-geometric>=2.3',
            'dgllife>=0.3',
        ],
        'docs': ['sphinx', 'sphinx-rtd-theme'],
        'openff': [
            'openff-toolkit>=0.14',
            'openff-interchange>=0.3',
        ],
        'dev': ['pytest>=7.0', 'pytest-cov'],
    },
    entry_points={
        'console_scripts': [
            # Legacy (compatible with original CLI)
            'genfeatures=dockml.bin.genfeatures:main',
            'gmx_cmap=bin.gmx_cmap:main',
            # Modern CLI
            'dml-prepare=dockml.modern.cli:prepare_main',
            'dml-dock=dockml.modern.cli:dock_main',
            'dml-features=dockml.modern.cli:features_main',
            'dml-train=dockml.modern.cli:train_main',
            'dml-rescore=dockml.modern.cli:rescore_main',
            'dml-mdtop=automd.cli:gentop_main',
            'dml-md=automd.cli:md_main',
        ],
    },
    zip_safe=False,
    python_requires='>=3.9',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
