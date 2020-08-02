import os, sys
from dockml import pdbIO
from dockml import algorithms
import numpy as np


class shiftPDB(pdbIO.rewritePDB):
    """
    Shift or rotate molecules in a pdb file

    Parameters
    ----------

    Methods
    -------

    Attributes
    ----------

    """

    def xyzReplace(self, pdbline, XYZ):
        """
        Replace the original coordinates in a pdb file

        Parameters
        ----------
        pdbline: str,
            the atom/hetatm line in a pdb file
        XYZ: list of float,
            the new coordinates vector

        Returns
        -------
        newline: str,
            the new pdb atom line

        """

        if len(XYZ) != 3:
            print("Number of values for XYZ vector not equals to 3. Exit Now!")
            return ""

        head = pdbline[:30]
        tail = pdbline[54:]

        x = XYZ[0]
        y = XYZ[1]
        z = XYZ[2]

        newline = head + '{:8.3f}'.format(x) + '{:8.3f}'.format(y) + '{:8.3f}'.format(z) + tail

        return newline

    def xyzChanger(self, pdbline, XYZ):
        """
        shift x y z by a vector

        Parameters
        ----------
        pdbline: str,
            a line in original pdb file
        XYZ: list of float,
            the new coordinates vector

        Returns
        -------
        newline: str,
            a new pdb line

        """

        # in PDB
        # X 31 - 38
        # Y 39 - 46
        # Z 47 - 54
        if len(XYZ) != 3:
            print("Number of values for XYZ vector not equals to 3. Exit Now!")
            return ""

        head = pdbline[:30]
        tail = pdbline[54:]

        x = XYZ[0] + float(pdbline[30:38].strip())
        y = XYZ[1] + float(pdbline[38:46].strip())
        z = XYZ[2] + float(pdbline[46:54].strip())

        newline = head + '{:8.3f}'.format(x) + '{:8.3f}'.format(y) + '{:8.3f}'.format(z) + tail

        return newline

    def xyzReverser(self, pdbline, XYZ):
        """
        shift x y z by a vector

        Parameters
        ----------
        pdbline: str,
            a line in original pdb file
        XYZ: list of float,
            the new coordinates vector

        Returns
        -------
        newline: str,
            a new pdb line

        """

        # in PDB
        # X 31 - 38
        # Y 39 - 46
        # Z 47 - 54
        if len(XYZ) != 3:
            print("Number of values for XYZ vector not equals to 3. Exit Now!")
            return ""

        head = pdbline[:30]
        tail = pdbline[54:]

        x = XYZ[0] * float(pdbline[30:38].strip())
        y = XYZ[1] * float(pdbline[38:46].strip())
        z = XYZ[2] * float(pdbline[46:54].strip())

        newline = head + '{:8.3f}'.format(x) + '{:8.3f}'.format(y) + '{:8.3f}'.format(z) + tail
        return newline


class simulationModeller(object):
    """
    Manipulate the molecules in a pdb file

    Parameters
    ----------
    pdb_in: str,
        input pdb file for manipulation

    Methods
    -------
    shiftXYZ


    Attributes
    ----------

    """

    def __init__(self, pdb_in):
        self.pdb = pdb_in

    def shift_xyz(self, vectors, output="new.pdb",
                  inplace=False, newchain="",
                  chainID=['A'], residueList=[],
                  res_seq="1"):
        """
        Shift molecules translational

        Parameters
        ----------
        vectors: list of float,
            the xyz shift vector
        output: str,
            output file name
        inplace: bool,
            whether modify the molecule in place
        chainID: list of str,
            the chain identifier for molecule manipulation
        residueList: list of integer
            the residue sequence numbers

        Returns
        -------

        """

        # duplicate the system
        tofile = open(output,'w')
        rwPBD = shiftPDB(self.pdb)

        with open(self.pdb) as lines:
            for s in lines:
                if len(s.split()) > 5 and \
                        s[:6] in ["ATOM  ", "HETATM"] and \
                        s[21] in chainID and \
                        int(s[22:26].strip()) in residueList:
                    if newchain == "":
                        newchain = s[21]

                    line = rwPBD.xyzChanger(s, vectors)
                    line = rwPBD.chainIDChanger(line, newchain)
                    line = rwPBD.resSeqChanger(line, res_seq)
                    tofile.write(line)
                else:
                    if inplace:
                        tofile.write(s)
                    else:
                        pass
        tofile.close()

    def reverse_xyz(self, vectors, output="new.pdb",
                    inplace=False,
                    chainID=['A'], residueList=[]):
        """
        Shift molecules translational

        Parameters
        ----------
        vectors: list of float,
            the xyz shift vector
        output: str,
            output file name
        inplace: bool,
            whether modify the molecule in place
        chainID: list of str,
            the chain identifier for molecule manipulation
        residueList: list of integer
            the residue sequence numbers

        Returns
        -------

        """

        # duplicate the system
        tofile = open(output,'w')
        rwPBD = shiftPDB(self.pdb)

        with open(self.pdb) as lines:
            for s in lines:
                if len(s.split()) > 5 and \
                        s[:6] in ["ATOM  ", "HETATM"] and \
                        s[21] in chainID and \
                        int(s[22:26].strip()) in residueList:
                    line = rwPBD.xyzReverser(s, vectors)
                    line = rwPBD.chainIDChanger(line, s[21])
                    tofile.write(line)
                else:
                    if inplace:
                        tofile.write(s)
                    else:
                        pass
        tofile.close()

    def trunctSimulationBox(self, pbcVector, output="box.pdb",
                            headAtom='P',
                            lipidRes=['DOE', 'DOP', 'LPS']):
        """
        remove the resides out of PBC box

        Parameters
        ----------
        pbcVector: list of float,
            the pbc box dimensions
        output: str,
            output pdb file name
        headAtom: str,
            the head atom for lipids
        lipidRes: list of str,
            the lipid molecule residue names

        Returns
        -------

        """

        tofile = open(output,'wb')
        removeResList = []
        removeAtomList = []
        with open(self.pdb) as lines:
            for s in lines:
                if len(s.split()) > 5 and \
                                s[:4] in ["ATOM", "HETA"] and \
                                s.split()[2] == headAtom and \
                                s.split()[3] in lipidRes:
                    if float(s[30:38].strip()) > pbcVector[0] or \
                            float(s[38:46].strip()) > pbcVector[1]:
                        removeResList.append(s.split()[3]+"_"+s.split()[4]+"_"+s.split()[5])
                        removeAtomList.append(s.split()[1])

                    if float(s[30:38].strip()) < 0 or \
                            float(s[38:46].strip()) < 0:
                        removeResList.append(s.split()[3] + "_" + s.split()[4] + "_" + s.split()[5])
                        removeAtomList.append(s.split()[1])

        with open(self.pdb) as lines:
            for s in lines:
                if len(s.split()) > 5 and \
                                s[:4] in ["ATOM", "HETA"] and \
                                s.split()[3] in lipidRes :
                    if s.split()[3]+"_"+s.split()[4]+"_"+s.split()[5] not in removeResList and \
                                    s.split()[1] not in removeAtomList:
                        tofile.write(s)
        tofile.close()


class shiftDNA(object):
    """
    Shift DNA or RNA moleulces along long axis

    Parameters
    ----------
    pdb_in: str,
        input pdb file name

    Methods
    -------

    Attributes
    ----------

    """

    def __init__(self, pdb_in):

        self.pdb = pdb_in

    def shiftDNA(self,dt=10.0, out = "output.pdb"):
        """
        shift DNA along its helical axis
        using all atoms in the DNA, and calculate the fitting line of the atoms
        then get the vector of the fitting line, and then shift the DNA along
        the vector with movement step size

        Parameters
        ----------
        dt: float,
            could be negative, unit Angstrom, step size of the movement
        out: str,
            output pdb file name

        Returns
        -------

        """

        if not os.path.exists(self.pdb):
            sys.exit(0)

        with open(self.pdb) as lines:
            lines = [x for x in lines if "ATOM" in x]

            crd = pdbIO.coordinatesPDB()

            p_crds = crd.getAtomCrdFromLines(lines)

        # calculate the vector of the DNA helix axis
        algo = algorithms.LineFit(np.array(p_crds))
        vectors = algo.fit_line()

        vectors = vectors * dt

        tofiles = open(out, 'w')
        shift = shiftPDB(self.pdb)

        with open(self.pdb) as lines:
            for atom in lines:
                if "ATOM" in atom:
                    newline = shift.xyzChanger(atom, vectors)
                    tofiles.write(newline)
                else:
                    tofiles.write(atom)

        tofiles.close()

        return 1


if __name__ == "__main__" :

    os.chdir(os.getcwd())

    #sm = SimulationModeller(sys.argv[1])
    #sm.reversXYZ([1,1,-1], "temp1.pdb", atomGroup=range(1, 40000))
    #sm = SimulationModeller("temp2.pdb")
#    for height in range(135, 155, 2) :
#        for y
    #sm.shiftXYZ([10,-35,20],output='temp7.pdb', atomGroup=range(1, 40000), residueList=['DOE','DOG', 'LPS'])
    '''
    sm = SimulationModeller("temp7.pdb")
    sm.shiftXYZ([-95.522, 0, 0], output='toplevel_1.pdb', atomGroup=range(1, 40000), residueList=['DOE', 'DOG', 'LPS'], chainID='B')
    sm.shiftXYZ([95.522, 0, 0], output='toplevel_2.pdb', atomGroup=range(1, 40000), residueList=['DOE', 'DOG', 'LPS'], chainID='C')
    sm.shiftXYZ([0, -95.522, 0], output='toplevel_3.pdb', atomGroup=range(1, 40000), residueList=['DOE', 'DOG', 'LPS'], chainID='D')
    sm.shiftXYZ([0, 95.522, 0], output='toplevel_4.pdb', atomGroup=range(1, 40000), residueList=['DOE', 'DOG', 'LPS'], chainID='E')
    sm.shiftXYZ([-95.522, 95.522, 0], output='toplevel_5.pdb', atomGroup=range(1, 40000), residueList=['DOE', 'DOG', 'LPS'], chainID='F')
    sm.shiftXYZ([95.522, 95.522, 0], output='toplevel_6.pdb', atomGroup=range(1, 40000),
                residueList=['DOE', 'DOG', 'LPS'], chainID='G')
    sm.shiftXYZ([95.522,-95.522, 0], output='toplevel_7.pdb', atomGroup=range(1, 40000),
                residueList=['DOE', 'DOG', 'LPS'], chainID='H')
    sm.shiftXYZ([-95.522,-95.522, 0], output='toplevel_8.pdb', atomGroup=range(1, 40000),
                residueList=['DOE', 'DOG', 'LPS'], chainID='I')
    '''

    sm = simulationModeller("large_8pcs.pdb")
    sm.trunctSimulationBox([95.522, 95.522, 180], output='final_box.pdb')
