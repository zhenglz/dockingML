# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import sklearn
from sklearn import decomposition
import mdtraj as mt
from mdanaly import gmxcli


class PCA(object):
    """
    A PCA module for data analysis.

    Parameters
    ----------
    n_components: int, default is 20.
        number of remain dimensions after PCA analysis

    Attributes
    ----------
    trained_: bool,
        whether the pca object has been trained
    X_transformed_: numpy ndarray, shape=[N, M]
        the transformed X dataset,
        N is number of samples, M is the reduced dimensions,
        M = n_components
    scaled_: bool,
        whether the X dataset has been scaled. It is required for
        PCA calculation.
    scaler_: sklearn preprocessing StandardScaler object
        the scaler object
    X_scaled_: numpy ndarray, shape=[N, M]
        the scaled X dataset,
        N is number of samples, M is the number of dimensions
    pca_obj: sklearn decomposition PCA object,
        the pca object
    eigvalues_: numpy array,
        the eigvalues of the new axis
    eigvalues_ratio: numpy array,
        the percentage ratio of the new axis variances
    eigvectors_: numpy ndarray,
        # TODO: to to completed here

    """

    def __init__(self, n_components=20):
        self.n_components = n_components

        # attributes
        self.trained_ = False
        self.X_transformed_ = None

        self.scaled_ = False
        self.scaler_ = None
        self.X_scaled = None

        self.pca_obj = None
        self.eigvalues_ = None
        self.eigvalues_ratio_ = None
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
            print("Dataset not scaled. Scaling the dataset now ...... ")
            self.scaler_ = sklearn.preprocessing.StandardScaler()
            Xs = self.scaler_.fit_transform(X)

            self.scaled_ = True
        self.X_scaled = Xs

        # using sklearn, perform PCA analysis based on scaled dataset
        print("Perform PCA decompostion now ...... ")
        pca_obj = decomposition.PCA(n_components=self.n_components)
        pca_obj.fit(self.X_scaled)

        # train and transform the dataset
        print("Transform dataset now ...... ")
        self.X_transformed_ = pca_obj.transform(X)
        self.trained_ = True

        self.pca_obj = pca_obj

        # get eigval and eigvect
        self.eigvalues()
        self.eigvectors()

        return pca_obj

    def transform(self, X):
        """

        Parameters
        ----------
        X: numpy, ndarray, shape = [ N, M]
            the input dataset, N is number of samples,
            M is number of dimensions

        Returns
        -------
        X_transformed: numpy, ndarray, shape=[N, M_1]
            transformed dataset, N is number of samples
            M_1 is number of reduced dimensions, M_1 == n_components

        """

        if self.trained_:
            return self.pca_obj.transform(X)
        else:
            print("Your pca object is not trained yet. Training it now ...")
            pca_obj = self.fit(X)

            return pca_obj.transform(X)

    def eigvalues(self):
        """
        Eigen values, the variance of the eigenvectors

        Returns
        -------
        eigval_: numpy array, shape=[1, M]
            the eigvalues, M is number of reduced dimensions.
            M == n_components

        """

        eigval_ = self.pca_obj.explained_variance

        self.eigvalues_ = eigval_
        self.eigvalues_ratio_ = self.pca_obj.explained_variance_ratio_

        return None

    def eigvectors(self):
        """
        Eigenvectors. This vectors could be applied for essential dynamics
        analysis.

        Returns
        -------
        eigvect_: numpy ndarray, shape=[M_1, M]
            # TODO: add further explanations here

        """

        eigvect_ = self.pca_obj.components_

        self.eigvectors_ = eigvect_

        return None


class CoordinationPCA(object):
    """
    Perform trajectory coordinates PCA analysis using mdtraj

    Parameters
    ----------
    traj: mdtraj.Trajectory object,
        a mdtraj trajectory, where the coordinates are stored
    top: str, format pdb
        the reference pdb file name
    atom_selection: str, default is name CA
        the atom selection for pca calculation.

    Attributes
    ----------
    topology: mdtraj.Trajectory.topology object
        the mdtraj trajectory topology object
    n_atoms_: int,
        number of atoms in the trajectory
    superimposed_: bool,
        whether the trajectory has been superimposed
    superpose_atom_indices: numpy array,
        the atom indices, starting from 0, format is int
    xyz: numpy ndarray, shape=[N, M]
        the original xyz coordinates of selected atoms
        N is number of samples, or frames
        M is the multiply of atoms * 3

    """

    def __init__(self, traj, top, atom_selection="name CA"):
        self.traj = traj
        self.ref = top

        self.topology = mt.load(self.ref).toology

        self.n_atoms_ = self.traj.n_atoms

        self.superimposed_ = False
        self.superpose_atom_indices_ = self.topology.select(atom_selection)

        self.xyz = None

    def superimpose(self):
        """
        Superimpose the trajectory to a reference structure.

        Returns
        -------
        self: the object itself

        """

        self.traj.superpose(self.ref, frame=0, atom_indices=self.superpose_atom_indices_)
        self.superimposed_ = True

        return self

    def xyz_coordinates(self, atom_indices=None):
        """
        Extract xyz coordinates for selected atoms for a trajectory object

        Parameters
        ----------
        atom_indices: numpy array,
            the atom index for selected atoms

        Returns
        -------
        xyz: numpy ndarray, shape = [N, M]
            N is number of samples, or frames
            M is the multiply of number of atoms and 3

        """

        if atom_indices is None:
            atom_indices = self.superpose_atom_indices_

        if not self.superimposed_:
            print("Trajectory is not superimposed. Superimpose it to reference now...")
            self.superimpose()

        traj = self.traj.atom_slice(atom_indices=atom_indices, inplace=False)

        xyz = traj.xyz

        xyz = np.reshape(xyz, (xyz.shape[0], xyz.shape[1]*3))

        self.xyz = xyz

        return xyz


def iterload_xyz_coordinates(xtcfile, top, chunk, stride, atom_selection="name CA"):
    """
    Load a large gromacs xtc trajectory file using mdtraj.iterload method, and extract
    the atom xyz coordinates.

    Parameters
    ----------
    xtcfile: str, format gromacs xtc
        the gromacs xtc trajectory file
    top: str, format pdb
        the reference pdb file
    chunk: int,
        number of frames to load per time
    stride: int,
        skip n frames when loading trajectory file
    atom_selection: str,
        the atom selection language

    Returns
    -------
    xyz: numpy array, shape = [N, M]
        the xyz coordinates dataset,
        N is number of samples, or frames
        M is n_atoms * 3

    """

    trajs = gmxcli.read_xtc(xtc=xtcfile, top=top, chunk=chunk, stride=stride)

    xyz = np.array([])

    for traj in trajs:
        copca = CoordinationPCA(traj, top, atom_selection)
        # copca.superimpose()

        if xyz.shape[0] == 0:
            xyz = copca.xyz_coordinates()
        else:
            xyz = np.concatenate((xyz, copca.xyz_coordinates()), axis=0)

    return xyz


def xyz_pca(args):
    """
    Perform xyz coordinate based PCA calculation based on Gromacs XTC files.


    Parameters
    ----------
    args: argparse object,
        the arguments holder

    Returns
    -------

    """

    # obtain all xyz coordinates in a long trajectory file
    xyz = iterload_xyz_coordinates(xtcfile=args.f, top=args.s,
                                   chunk=1000, stride=int(args.dt/args.ps),
                                   atom_selection=args.select)

    # perform PCA calculation using xyz coordinates
    pca = PCA(n_components=args.proj)
    pca.fit(xyz)
    xyz_transformed = pca.X_transformed_

    write_results(X_transformed=xyz_transformed, variance_ratio=pca.eigvalues_ratio_,
                  X_out=args.o, variance_out=args.var_ratio, dt=args.dt)

    return None


def write_results(X_transformed, variance_ratio, X_out, variance_out, dt):
    """
    Write PCA results into files:
    1. transformed dataset
    2. explained variance ratio

    Parameters
    ----------
    X_transformed: numpy ndarray
    variance_ratio: numpy array
    X_out: str
    variance_out: str
    dt: int

    Returns
    -------

    """
    # save the data into a file
    dat = pd.DataFrame(X_transformed)
    dat.index = np.arange(X_transformed.shape[0]) * dt
    dat.columns = ["PC_%d" % x for x in np.arange(X_transformed.shape[1])]
    dat.to_csv(X_out, sep=",", header=True, index=True,
               index_label="time(ps)", float_format="%.3f")

    # save variance ratio into a file
    variance = list(variance_ratio)
    eigval = pd.DataFrame()
    eigval["PC"] = np.arange(len(variance))
    eigval["eigval_ratio"] = variance
    eigval.to_csv(variance_out, sep=",", header=True, index=False, float_format="%.3f")

    return None


def general_pca(args):
    """
    Perform PCA calculation given a clean dataset file.

    Parameters
    ----------
    args: argparse object
        the arguments holder

    Returns
    -------

    """

    X = pd.read_csv(args.f, sep=",", header=0)

    pca = PCA(n_components=args.proj)

    pca.fit(X)

    X_transformed = pca.X_transformed_
    variance_ratio = pca.eigvalues_ratio_

    write_results(X_transformed=X_transformed, variance_ratio=variance_ratio,
                  X_out=args.o, variance_out=args.var_ratio, dt=args.dt)
    return None


def arguments():
    # prepare gmx style argument for pca calculation
    d = """
    Perform PCA analysis of the xyz coordinates of selected atoms

    Example: 

    """
    parser = gmxcli.GromacsCommanLine(d=d)

    parser.arguments()

    parser.parser.add_argument("-mode", type=str, default="general",
                               help="Input, optional. \n"
                                    "The PCA calculation mode. \n"
                                    "Options: general, xyz, cmap \n"
                                    "general: perform general PCA using a well formated dataset file. \n"
                                    "xyz: perform xyz coordinate PCA analysis using a trajectory xtc file. \n"
                                    "cmap: perform contact map PCA analysis using a trajectory xtc file. ")
    parser.parser.add_argument("-proj", type=int, default=2,
                               help="Input, optional. \n"
                                    "How many number of dimensions to output. \n"
                                    "Default is 2.")
    parser.parser.add_argument("-select", type=str, default="name CA",
                               help="Input, optional, it works with mode = xyz \n"
                                    "Atom selected for PCA calculation. \n"
                                    "Default is name CA.")
    parser.parser.add_argument("-var_ratio", type=str, default="explained_variance_ratio.dat",
                               help="Output, optional. \n"
                                    "Output file name containing explained eigen values variance ratio. \n"
                                    "Default is explained_variance_ratio.dat. ")

    parser.parse_arguments()
    args = parser.args

    return args


def main():

    args = arguments()

    if args.mode == "xyz":
        xyz_pca(args=args)

        return None

    elif args.mode == "general":
        general_pca(args=args)

        return None

    elif args.mode == "cmap":
        # TODO: to be completed.

        return NotImplementedError

