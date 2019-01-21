from mdanaly import gmxcli, pca, cmap
from dockml import index
import sys
import numpy as np

if __name__ == "__main__":

    ref = sys.argv[1]
    xtc = sys.argv[2]
    ind = sys.argv[3]

    # load xtc first
    trajs = gmxcli.read_xtc(xtc, ref, 1000, 50)

    # process index
    ndx = index.GmxIndex(ind)
    sets = ["receptor", "ligand"]
    used_groups = []

    atom_indices = []
    for i in [0, 1]:
        print("Please select a group for %s: " % sets[i])
        for j, gn in enumerate(ndx.groups):
            print("%d : %s" % (j, gn))

        used_groups.append(ndx.groups[int(input("Your choice: "))])

    rec_ndx = [int(x)-1 for x in ndx.groupContent(used_groups[0])]
    lig_ndx = [int(x)-1 for x in ndx.groupContent(used_groups[1])]

    cmaps = np.array([])
    # generate cmap now
    for traj in trajs:
        CMap = cmap.ContactMap(traj, rec_ndx, lig_ndx, cutoff=0.5)
        CMap.generate_atom_pairs()
        CMap.generate_cmap(switch=True)

        c = CMap.cmap_

        if cmaps.shape[0] == 0:
            cmaps = c
        else:
            cmaps = np.concatenate((cmaps, c), axis=0)

    # run PCA now
    pca_obj = pca.PCA()
    pca_obj.fit_transform(cmaps[:, 1:])
    print(pca_obj.eigvalues_ratio_)
    X_transformed = pca_obj.X_transformed_
    #X_transformed = np.concatenate()
    np.savetxt("transformed_X.csv", X_transformed, delimiter=",", fmt="%8.3f")
    print("PCA calculation completed. ")

