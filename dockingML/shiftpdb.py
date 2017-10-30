import os
from .pdbIO import RewritePDB

class ShiftPDB(RewritePDB) :
    """
    Shift or rotate pdb files
    """

    def xyzChanger(self, pdbline, XYZ ):
        """
        shift x y z by a vector
        :param pdbline: str, a line in original pdb file
        :param XYZ: list of floats, added parameters for x y z
        :return: str, a new pdb line
        """
        # in PDB
        # X 31 - 38
        # Y 39 - 46
        # Z 47 - 54
        if len(XYZ) != 3 :
            print("Number of values for XYZ vector not equals to 3. Exit Now!")
            return ""

        head = pdbline[:30]
        tail = pdbline[54:]

        x = XYZ[0] + float(pdbline[30:38].strip())
        y = XYZ[1] + float(pdbline[38:46].strip())
        z = XYZ[2] + float(pdbline[46:54].strip())

        return head + '{:8.3f}'.format(x) + '{:8.3f}'.format(y) + '{:8.3f}'.format(z) + tail

    def xyzReverser(self, pdbline, XYZ):
        """
        rotate x y z by by a vector
        :param pdbline: str, a line in original pdb file
        :param XYZ: list of floats, multiplied parameters for x y z
        :return: str, a new pdb line
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

        return head + '{:8.3f}'.format(x) + '{:8.3f}'.format(y) + '{:8.3f}'.format(z) + tail

class SimulationModeller :

    def __init__(self, PDBin):
        self.pdb = PDBin

    def shiftXYZ(self, vectors, output="new.pdb",
                 append=False, atomGroup=[1],
                 chainID='A', residueList=[]) :
        """

        :param vectors:
        :param output:
        :param append:
        :param atomGroup:
        :param chainID:
        :param residueList:
        :return:
        """
        # duplicate the system
        if append :
            tofile = open(output,'w')
            with open(self.pdb) as lines :
                for s in lines :
                    if "END" not in s :
                        tofile.write(s)
            #tofile.close()
        else :
            tofile = open(output,'w')

        rwPBD = ShiftPDB(self.pdb)

        with open(self.pdb) as lines :
            for s in lines :
                if len(s.split()) > 5 and s[:6] in ["ATOM  ", "HETATM"] :
                    if int(s.split()[1]) in atomGroup :
                        if len(residueList) == 0 :
                            line = rwPBD.xyzChanger(s, vectors)
                            line = rwPBD.chainIDChanger(line, chainID)

                            tofile.write(line)

                        else :
                            if s[17:20] in residueList :
                                line = rwPBD.xyzChanger(s, vectors)
                                line = rwPBD.chainIDChanger(line, chainID)

                                tofile.write(line)
                    else :
                        pass
                        #tofile.write(s)
                else:
                    tofile.write(s)
        tofile.close()

    def reversXYZ(self, vectors=[1, 1, 1], output="new.pdb",
                  append=False, atomGroup=[1, ], chainID='A'):

        # duplicate the system
        if append:
            tofile = open(output, 'w')
            with open(self.pdb) as lines:
                for s in lines:
                    if "END" not in s:
                        tofile.write(s)
            tofile.close()
        else:
            tofile = open(output, 'w')

        rwPBD = ShiftPDB(self.pdb)

        with open(self.pdb) as lines:
            for s in lines:
                if len(s.split()) > 5 and \
                                s[:6] in ["ATOM  ", "HETATM"] and \
                                int(s.split()[1]) in atomGroup:
                    line = rwPBD.xyzReverser(s, vectors)
                    line = rwPBD.chainIDChanger(line, chainID)

                    tofile.write(line)
                else:
                    tofile.write(s)
        tofile.close()

    def trunctSimulationBox(self, pbcVector, output="box.pdb",
                            headAtom='P',
                            lipidRes=['DOE', 'DOP', 'LPS']):
        """
        remove the resides out of PBC box
        :param pbcVector:
        :param output:
        :param headAtom:
        :param lipidRes:
        :return:
        """
        tofile = open(output,'wb')
        removeResList = []
        removeAtomList = []
        with open(self.pdb) as lines :
            for s in lines :
                if len(s.split()) > 5 and \
                                s[:4] in ["ATOM", "HETA"] and \
                                s.split()[2] == headAtom and \
                                s.split()[3] in lipidRes :
                    if float(s[30:38].strip()) > pbcVector[0] or float(s[38:46].strip()) > pbcVector[1] :
                        removeResList.append(s.split()[3]+"_"+s.split()[4]+"_"+s.split()[5])
                        removeAtomList.append(s.split()[1])

                    if float(s[30:38].strip()) < 0 or float(s[38:46].strip()) < 0 :
                        removeResList.append(s.split()[3] + "_" + s.split()[4] + "_" + s.split()[5])
                        removeAtomList.append(s.split()[1])

        with open(self.pdb) as lines :
            for s in lines :
                if len(s.split()) > 5 and \
                                s[:4] in ["ATOM", "HETA"] and \
                                s.split()[3] in lipidRes :
                    if s.split()[3]+"_"+s.split()[4]+"_"+s.split()[5] not in removeResList and \
                                    s.split()[1] not in removeAtomList :
                        tofile.write(s)
        tofile.close()

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

    sm = SimulationModeller("large_8pcs.pdb")
    sm.trunctSimulationBox([95.522, 95.522, 180], output='final_box.pdb')
