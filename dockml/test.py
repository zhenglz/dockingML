#from dockml import pdbIO, index
#from mdanaly import cmap

import sys


def rewritePDB(pdb, chain="B") :

    rew = pdbIO.rewritePDB(pdb)

    rew.pdbRewrite(chain+pdb, chain, 1, 1)

def contacts() :
    inp = sys.argv[1]
    #out = sys.argv[2]

    cutoff = 10.0

    rec_res = 126
    lig_res = 7

    # split file now
    condition = True
    rec = open("rec.pdb", 'w')
    lig = open("lig.pdb", 'w')
    with open(inp) as lines :
        for s in lines :
            if "ACE" in s :
                condition = False

            if condition :
                rec.write(s)
            else :
                lig.write(s)
    rec.close()
    lig.close()

    coord = pdbIO.coordinatesPDB()
    index = index.PdbIndex()

    rec_ndx = index.res_index("rec.pdb","A", "non-hydrogen", [rec_res], [], "none")
    rec_coord = coord.getAtomCrdByNdx("rec.pdb", rec_ndx)

    lig_ndx = index.res_index("lig.pdb","A", "non-hydrogen", [lig_res], [], "none")
    lig_coord = coord.getAtomCrdByNdx("lig.pdb", lig_ndx)

    CMAP = cmap.ContactMap(inp)
    nbyn = CMAP.residueContacts(rec_coord, lig_coord, cutoff, 1, rank=0, NbyN=True)

    print(nbyn)

def getMMPBSA() :

    import numpy as np
    import os

    dataset = []
    for i in range(1, 501) :
        f = "complex_%d/summary_energy.dat" % i
        if os.path.exists(f) :
            energy = []
            energy.append(i)
            with open(f) as lines :
                for s in lines :
                    if "Waal" in s :
                        energy.append(float(s.split()[5]))
                    elif "Electrostattic" in s :
                        energy.append(float(s.strip()[3]))
                    elif "Polar" in s :
                        energy.append(float(s.split()[4]))
                    elif "SASA" in s :
                        energy.append(float(s.split()[3]))
                    elif "Binding" in s :
                        energy.append(float(s.split()[3]))

            dataset.append(energy)

    dataset = np.array(dataset)

    np.savetxt("mmpbsa_summary.dat", dataset, delimiter=" ", fmt="%8.3f")

def findAtomNdx(pdbfile, resList, chain, atomType, verbose=False):
    '''
    give a pdb file, return the atom ndx needed
    :param pdbfile:
    :param resList: a list of index of residues
    :param chain: chains
    :param atomType: atom names
    :param verbose:
    :return: a list of atom ndx, string
    '''
    if verbose :
        print( pdbfile, resList, chain, atomType)

    atomndx = []
    for key in resList.keys() :
        for ndx in resList[key] :
            with open(pdbfile) as lines :
                atomndx +=     [ s.split()[1]
                                 for s in lines
                                 if len(s.split()) and
                                 s[:6].strip() in ["ATOM", "HETATM"] and
                                 s.split()[2] in atomType and
                                 s[21] in list(chain) and
                                 int(s[22:26].strip()) in ndx and
                                 s.split()[1] not in atomndx ]
            '''
            for s in lines :
                if "ATOM" in s or "HETATM" in s :
                    if s.split()[2] in atomType and s[21] in list(chain) and
                        int(s[22:26].strip()) in resList[key] :
                        if s.split()[1] not in atomndx :
                            append(s.split()[1])
                        else :
                            pass '''
    return atomndx



def summpbsa():
    import numpy as np
    import subprocess as sp
    from collections import defaultdict
    from mdanaly import cmap

    ref = "complex_1/em-c1.pdb"
    sp.check_output("echo \"1 0 \" | trjconv -f complex_1/em-c1.gro -s complex_1/em-c1.tpr -o complex_1/em-c1.pdb",
                    shell=True)
    CMAP = cmap.ContactMap("./complex_1/em-c1.pdb")

    at = CMAP.findAtomType('non-hydrogen', ref)
    rd = defaultdict(list)
    rd["A"] = range(1, 300)
    ld = defaultdict(list)
    ld["B"] = range(1, 10)
    rndx = CMAP.findAtomNdx(ref, rd, "A", at, False)
    lndx = CMAP.findAtomNdx(ref, ld, "B", at, False)

    energy = np.loadtxt("./mmpbsa_summary.dat")

    f_idx = [int(x) for x in energy[:, 0]]

    result = np.array([])

    for mmpbsa in energy:
        i = int(mmpbsa[0])
        sp.check_output(
            "echo \"1 0 \" | trjconv -f complex_%d/em-c%d.gro -s complex_%d/em-c%d.tpr -o complex_%d/em-c%d.pdb"
            % (i, i, i, i, i, i), shell=True)
        fn = ["complex_%d/em-c%d.pdb" % (i, i)]

        c = CMAP.cmap_ca(fn, 0.35, False, [rndx, lndx], 0, False,
                         len(rndx) * len(rndx), perAtom=[False, False],
                         ccutoff=2.0, NbyN=False)

        e = mmpbsa[-1]

        if not result.shape[0]:
            result = c * (2.71828 ** (-1.0 * e / 2.5))
        else:
            result += c * (2.71828 ** (-1.0 * e / 2.5))

    MAX = np.max(result)
    cmap = np.reshape(result, (len(lndx), len(rndx)))

    np.savetxt("cmap_energy_normalized.dat", cmap / MAX, delimiter=" ", fmt="%.4f")


def test() :
    import numpy as np

    # Generate some data that lies along a line

    x = np.mgrid[-2:5:120j]
    y = np.mgrid[1:9:120j]
    z = np.mgrid[-5:3:120j]

    print(x, y, z)

    data = np.concatenate((x[:, np.newaxis],
                           y[:, np.newaxis],
                           z[:, np.newaxis]),
                          axis=1)

    print("DATA " * 10)
    print(data.shape)
    # Perturb with some Gaussian noise
    data += np.random.normal(size=data.shape) * 0.4

    # Calculate the mean of the points, i.e. the 'center' of the cloud
    datamean = data.mean(axis=0)

    # Do an SVD on the mean-centered data.
    uu, dd, vv = np.linalg.svd(data - datamean)

    print("VVV " * 50)
    print(vv)
    print(vv.shape)

    linepts = vv[0] * np.mgrid[-7:7:2j][:, np.newaxis]

    # shift by the mean to get the line in the right place
    linepts += datamean

    # Verify that everything looks right.

    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d as m3d

    ax = m3d.Axes3D(plt.figure())
    ax.scatter3D(*data.T)
    ax.plot3D(*linepts.T)
    plt.show()


def shiftDNA(input) :
    import os, sys
    from dockml import pdbIO
    from dockml import algorithms
    import numpy as np
    from automd import shiftpdb

    dt = 10 # angstrom
    out = "output.pdb"

    if not os.path.exists(input) :
        sys.exit(0)

    with open(input) as lines :
        lines = [ x for x in lines if "ATOM" in x ]

        crd = pdbIO.coordinatesPDB()

        p_crds = crd.getAtomCrdFromLines(lines)

    algo = algorithms.LineFit(np.array(p_crds))
    vectors = algo.fit_line()

    vectors = vectors * dt

    tofiles = open(out, 'w')

    shift = shiftpdb.shiftPDB(input)

    with open(input) as lines :
        for atom in lines :
            if "ATOM" in atom :
                newline = shift.xyzChanger(atom, vectors)
                tofiles.write(newline)
            else :
                tofiles.write(atom)

    tofiles.close()

shiftDNA("input.pdb")
