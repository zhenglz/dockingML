import numpy as np
import pandas as pd
import sklearn
import mdtraj as mt
from mdanaly import gmxcli

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
        pca_obj = sklearn.decomposition.PCA(n_components=self.n_components)
        pca_obj.fit(self.X_scaled)

        # train the dataset
        self.trained_ = True
        self.X_transformed_ = pca_obj.transform(X)

        self.pca_obj = pca_obj

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

        return eigval_

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

        return eigvect_


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

        self.topology = mt.load(self.ref).topology

        self.n_atoms_ = self.traj.n_atoms

        self.superimposed_ = False
        self.superpose_atom_indices_ = self.top.select(atom_selection)

        self.xyz = None

    def superimpose(self):

        self.traj.superpose(self.top, frame=0, atom_indices=self.superpose_atom_indices_)
        self.superimposed_ = True

        return self

    def xyz_coordinates(self, atom_indices=None):

        if atom_indices == None:
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

    trajs = gmxcli.read_xtc(xtc=xtcfile, top=top, chunk=chunk, stride=stride)

    xyz = np.array([])

    for traj in trajs:
        copca = CoordinationPCA(traj, top, atom_selection)
        #copca.superimpose()

        if xyz.shape[0] == 0:
            xyz = copca.xyz_coordinates()
        else:
            xyz = np.concatenate((xyz, copca.xyz_coordinates()), axis=0)

    return xyz


def xyz_pca():

    # prepare gmx style argument for pca calculation
    d = """
    Perform PCA analysis of the xyz coordinates of selected atoms
    
    Example: 
    
    """
    parser = gmxcli.GromacsCommanLine(d=d)

    parser.arguments()

    parser.parser.add_argument("-proj", type=int, default=2,
                               help="Input, optional. \n"
                                    "How many number of dimensions to output. \n"
                                    "Default is 2.")
    parser.parser.add_argument("-select", type=str, default="name CA",
                               help="Input, optional. \n"
                                    "Atom selected for PCA calculation. \n"
                                    "Default is name CA.")
    parser.parser.add_argument("-var_ratio", type=str, default="explained_variance_ratio.dat",
                               help="Output, optional. \n"
                                    "Output file name containing explained eigen values variance ratio. \n"
                                    "Default is explained_variance_ratio.dat. ")

    parser.parse_arguments()
    args = parser.args

    # obtain all xyz coordinates in a long trajectory file
    xyz = iterload_xyz_coordinates(xtcfile=args.f, top=args.s,
                                   chunk=1000, stride=int(args.dt/args.ps),
                                   atom_selection=args.select)

    # perform PCA calculation using xyz coordinates
    pca = PCA(n_components=args.proj)
    pca.fit(xyz)
    xyz_transformed = pca.X_transformed_

    # save the data into a file
    dat = pd.DataFrame(xyz_transformed)
    dat.index = np.arange(xyz_transformed.shape[0]) * args.dt
    dat.columns = ["PC_%d" % x for x in np.array(xyz_transformed.shape[1])]
    dat.to_csv(args.o, sep=",", header=True, index=True,
               index_label="time(ps)", float_format="%.3f")

    # save variance ratio into a file
    variance = list(pca.eigvalues_ratio_)
    eigval = pd.DataFrame()
    eigval["PC"] = np.arange(len(variance))
    eigval["eigval_ratio"] = variance
    eigval.to_csv(args.var_ratio, sep=",", header=True, index=False, float_format="%.3f")

    return None




