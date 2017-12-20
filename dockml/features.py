# -*- coding: utf-8 -*-

from collections import defaultdict
import math
import argparse
from argparse import RawTextHelpFormatter
import sys, os
import numpy as np

class BindingFeature:

    def __init__(self):
        if os.path.exists("AtomType.dat") :
            self.atomtype = "AtomType.dat"
        else :
            PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
            DEFINITIONS_ROOT = os.path.join(PROJECT_ROOT, '../data/AtomType.dat')
            self.atomtype = DEFINITIONS_ROOT

    def getVdWParams(self, inputF) :
        '''
        Get vdw parameters, mainly sigma and epsilon term in energy function
        inputFile: str, a library file containing atomic vdw parameters
        :return:
        '''
        # obtain vdw radius for different atoms
        vdwParams = defaultdict(list)
        with open(inputF) as lines :
            for s in lines :
                #param = []
                if "#" not in s and ";" not in s :
                    #print s.split()[-2].split("e")[0]
                    ## in AtomType.dat nm
                    ## nm to angstrom by multiply 10.0
                    sigma   = float(s.split()[-2]) #.split("e")[0]) * ( 10 ** float(s.split()[-2].split("e")[1])) * 10.0
                    epsilon = float(s.split()[-1]) #.split("e")[0]) * ( 10 ** float(s.split()[-1].split("e")[1])) * 10.0
                    #param = [ sigma, epsilon ]
                    vdwParams[s.split()[0]] = [sigma, epsilon]

        # format of vdwParams : { "C":[C12, C6] }
        return vdwParams

    def getXYZCoord(self, pdbline) :
        '''
        from a pdb line to get xyz data
        :param pdbline:
        :return: a list of x y z data
        '''
        # obtain the X Y Z coordinates of a line from a pdb file, or a pdbqt file
        # both protein and ligand XYZ coordinates

        if "ATOM" not in pdbline and "HETATM" not in pdbline :
            print("Error! \nThis is not a PDB coordinate record! \nIgnore!")
            #sys.exit(1)
            return [0.0, 0.0, 0.0]
        else :
            # the x y z coordinates
            return [float(pdbline[30:38].strip()),
                    float(pdbline[38:46].strip()),
                    float(pdbline[46:54].strip())]

    def getAtomInfor(self, input, ligCode) :
        '''
        from a pdb file get ligand and receptor information
        input a pdbfile of the receptor-ligand complex
        return their xyz information and the detail information of each atom
        :param input: str, input receptor-ligand complex
        :param ligCode: str, res-name of the ligand
        :return: recInfor, ligInfor, recatomLine, ligatomLine
        '''
        recInfor = defaultdict(list)
        recatomLine = {}
        ligatomLine = {}
        ligInfor = defaultdict(list)

        with open(input) as lines :
            lines = [ s for s in lines if ("ATOM" == s.split()[0]
                                           or "HETATM" == s.split()[0]
                                           and len(s.split()) > 5)
                      ]
            for s in lines :
                #if "ATOM" == s.split()[0] or "HETATM" == s.split()[0] and len(s.split()) > 5 :
                ## atomID: AtomName_ResidueName_ResidueIndex_Chain
                atomID = s[6:11].strip() + "_" +s[17:20].strip() + "_" + s[22:26].strip()+"_"+s[21]

                ## receptor information
                if ligCode !=  s[17:20].strip() :
                    recInfor[ atomID ] = self.getXYZCoord(s)
                    recatomLine[atomID] = s

                ## ligand information
                elif ligCode ==  s[17:20].strip() :
                    ligInfor[ atomID ] = self.getXYZCoord(s)
                    ligatomLine[atomID] = s

        return recInfor, ligInfor, recatomLine, ligatomLine

    def atomDistance(self, xyz1, xyz2) :
        '''
        calculate distance between two atoms
        :param xyz1: list, x y z float coordiates
        :param xyz2: list, x y z float coordiates
        :return: float
        '''
        # input xyz for two atoms, return their distance

        if len(xyz1) == len(xyz2) :
            distance = sum(map(lambda x, y: (x-y)**2, xyz1, xyz2))
        else :
            print("Number of coordinates is not correct!")
            distance = 9999.0
            #sys.exit(1)
        return math.sqrt(distance)

    def atomDistMatrix(self, receptorXYZ, ligandXYZ) :
        '''
        calculate the distances of all the atom pairs
        receptorXYZ and ligandXYZ are dictionaries containing atom coordinates
        :param receptorXYZ: list 2D, a list of coordinates of atoms
        :param ligandXYZ: list 2D, a list of coordinates of atoms
        :return:
        '''
        allDistances = {}
        rkeys = receptorXYZ.keys()
        lkeys = ligandXYZ.keys()

        for rec in rkeys :
            for lig in lkeys :
                allDistances[ rec + "+" + lig ] = self.atomDistance( receptorXYZ[rec], ligandXYZ[lig])

        return allDistances

    def switchFuction(self, x, d0, m=12, n=6):
        """
        for countting, implement a rational switch function to enable a smooth transition
        the function is lik  s= [1 - (x/d0)^6] / [1 - (x/d0)^12]
        d0 is a cutoff, should be twice the larget than the distance cutoff
        :param x: float
        :param d0: distance cutoff, should be 2 times of normal cutoff
        :param m: int
        :param n: int
        :return: float
        """
        return(1.0 - math.pow((x / d0), n)) / (1.0 - math.pow((x / d0), m))

    def residueCounts(self, alldistpairs, atomDetailInfor, distanceCutoff = 3.5, shortRangeCutoff=12.0) :
        """
        input all the distances pairs
        return the number of contacts
        alldistpairs: the porteinID-LigID-Distance inforamtion
        atomdetailInfor : the atomID-PDBline information
        distance unit is angstrom
        :param alldistpairs: dictionary, all distances
        :param atomDetailInfor: dictionary, detail information
        :param distanceCutoff: float, distance cutoff
        :return:
        """

        atoms = alldistpairs.keys()
        recRes =[]

        backboneAtoms = ["C", "N", "O", "CA"]

        for atom in atoms :
            res = atom.split("+")[0].split("_")[1] + "_" + \
                  atom.split("+")[0].split("_")[2] + "_" + \
                  atom.split("+")[0].split("_")[3]

            # get the full residue_seq list
            if res not in recRes :
                recRes.append(res)

        # counts between backbone (or side-chain) and ligand
        backCount = {}  # format is : {"ARE175": 0 }
        sideCount = {}  # format is : {"ARE175": 0 }
        for res in recRes :
            backCount[res] = 0.0
            sideCount[res] = 0.0

        # find all the
        for atom in atoms :
            # first set a distance cutoff check, only short range distances are considered
            distance = alldistpairs[atom]
            if distance <= shortRangeCutoff :
                line = atomDetailInfor[atom.split("+")[0]]

                ### add more chain ID, in case you have several different chains
                ### better, you may use 1 2 3 4 as your chain id, here choose original chain id
                ### residueID  residueName_ResidueIndex_Chain
                residueID = atom.split("+")[0].split("_")[1] + "_" + \
                            atom.split("+")[0].split("_")[2] + "_" + \
                            atom.split("+")[0].split("_")[3]
                #residueID = line[17:20].strip() + "_" + line[22:26].strip()
                ## count contacts in backbones
                if line.split()[2] in backboneAtoms :
                    backCount[residueID] += self.switchFuction(distance, distanceCutoff*2.0)
                ## count contacts in sidechains
                elif line.split()[2] not in backboneAtoms :
                    sideCount[residueID] += self.switchFuction(distance, distanceCutoff*2.0)

        return recRes, backCount, sideCount

    def atomicVdWEnergy(self, atomtype1, atomtype2, distance, vdwParam) :
        """
        Calculate Van der Waals interaction energy
        VdW potential calculation, using combination rule 2
        vdwParam is a dictionary list, eg. { "O": [ 0.22617E-02, 0.74158E-0.6 ]; }
        :param atomtype1:
        :param atomtype2: dictionary, parameters of atom type and sigma epsilon
        :param distance: float
        :param vdwParam: dictionary of lists
        :return: float, short range vdw energy
        """

        Sigma_i   = vdwParam[atomtype1][0]
        Epsilon_i = vdwParam[atomtype1][1]
        Sigma_j   = vdwParam[atomtype2][0]
        Epsilon_j = vdwParam[atomtype2][1]

        # combination rule 2 of vdw MM force field
        sigma_ij   = 0.5 * ( Sigma_i + Sigma_j )
        epsilon_ij = math.sqrt( Epsilon_i * Epsilon_j )

        # Lennard-Jones Potential energy
        # the sigma is the distance at the lowest energy point, in unit of angstrom
        # (sigma / distance)^6
        C6  = 4.0 * epsilon_ij * (sigma_ij**6)
        C12 = 4.0 * epsilon_ij * (sigma_ij**12)
        d6  = distance ** 6
        d12 = distance ** 12

        energy = C12 / d12 - C6 / d6

        return energy

    def resVdWContribution(self, alldistpairs, recatomInfor, ligatomInfor,
                           vdwParams, maxCutoff, pdbqt=True) :
        """
        accumulated vdw for side-chain or backbone
        alldistpairs: the porteinID-LigID-Distance inforamtion
        atomdetailInfor : the atomID-PDBline information
        :param alldistpairs:
        :param recatomInfor:
        :param ligatomInfor:
        :param vdwParams:
        :param maxCutoff:
        :param pdbqt:
        :return:
        """

        atoms = alldistpairs.keys()
        recRes =[]

        backboneAtoms = ["C", "N", "O", "CA"]

        for atom in atoms :
            res = atom.split("+")[0].split("_")[1] + "_" + \
                  atom.split("+")[0].split("_")[2] + "_" + \
                  atom.split("+")[0].split("_")[3]
            if res not in recRes :
                recRes.append(res)  # get the full residue_seq_chain list

        backVdWEner = {}  # format is : {"ARE_175_A": -1.20 }
        sideVdWEner = {}  # format is : {"ARE_175_A": -1.20 }
        for res in recRes :
            backVdWEner[res] = 0.0
            sideVdWEner[res] = 0.0

        for atom in atoms :

            ## in recatomInfor, the whole original line of a pdb file is maintained, so as ligatomInfor
            recline = recatomInfor[atom.split("+")[0]]
            ligline = ligatomInfor[atom.split("+")[1]]

            residueID = atom.split("+")[0].split("_")[1] + "_" + atom.split("+")[0].split("_")[2] + "_" \
                        +  atom.split("+")[0].split("_")[3]

            # atom element, related to the vdw C6 and C12 terms
            #atomtype1 = proAtomTypes[ recline.split()[2] + "_" + recline[17:20].strip() ]
            # if this is a pdbqt file line
            if pdbqt :
                atomtype1 = recline.split()[-2]
                atomtype2 = ligline.split()[-2]
            else :
                # how to define the ligand atom type? Very difficult. May need to transfer file into Mol2 file
                # then decide the atom type from the mol2 file
                # atomtype2 = None  # atom element, related to the vdw radius
                # or, the case is, if there is no element information in the last column of the pdb line?
                # need to be completed here
                atomtype1 = recline.split()[-1]
                atomtype2 = ligline.split()[-1]

            if atomtype1 not in vdwParams.keys() :
                if atomtype1 == " N1+" :
                    atomtype1 = "N"
                atomtype1 = "DU"
            if atomtype2 not in vdwParams.keys() :
                if atomtype2 == " N1+":
                    atomtype2 = "N"
                atomtype2 = "DU"

            # calculate L J potential per residue
            if alldistpairs[atom] <= maxCutoff:
                # transfer angstrom to nanometer by multiplying 0.1
                energy = self.atomicVdWEnergy(atomtype1, atomtype2, alldistpairs[atom] * 0.1, vdwParams)
            else :
                energy = 0.0
                #energy = self.atomicVdWEnergy(atomtype1, atomtype2, maxCutoff, vdwParams)

            if recline.split()[2] in backboneAtoms :
                backVdWEner[residueID] += energy
            else :
                sideVdWEner[residueID] += energy

        return recRes, backVdWEner, sideVdWEner

    def contactsAtomtype(self, alldistpairs, recatomInfor, ligatomInfor,
                         vdwParams, distanceCutoff, pdbqt=True):
        # vdw counting based on different atom types
        # alldistpairs: the porteinID-LigID-Distance inforamtion
        # atomdetailInfor : the atomID-PDBline information
        # vdwParams is the C6 and C12 terms for different atom types
        atoms = alldistpairs.keys()
        atomTypeCounts = {}
        for atomT1 in vdwParams.keys() :
            for atomT2 in vdwParams.keys() :
                atomTypeCounts[ atomT1 + "_" + atomT2 ] = 0.0

        for atom in atoms :
            distance = alldistpairs[atom]
            if distance <= distanceCutoff * 4 :
                recline = recatomInfor[atom.split("+")[0]]
                ligline = ligatomInfor[atom.split("+")[1]]

                # unrecognized atomtype occur, "dummy" atomtype is used
                # create all the combinations of the atomtypes known
                if pdbqt :
                    ligAtom = ligline.split()[-2]
                else :
                    ligAtom = ligline.split()[-1]

                if ligAtom not in vdwParams.keys() :
                    if ligAtom == "N1+" :
                        ligAtom = "N"
                    else :
                        ligAtom = "DU"

                if pdbqt :
                    recAtom = recline.split()[-2]
                else :
                    recAtom = recline.split()[-1]

                if recAtom not in vdwParams.keys() :
                    if recAtom == "N1+" :
                        recAtom = "N"
                    else :
                        recAtom = "DU"

                atomTypeCombination = recAtom + "_" + ligAtom
                if atomTypeCombination not in atomTypeCounts.keys() :
                    # apply a switch function to smooth the transition
                    atomTypeCounts[atomTypeCombination] = self.switchFuction(distance, distanceCutoff*2.0)
                else :
                    atomTypeCounts[atomTypeCombination] += self.switchFuction(distance, distanceCutoff*2.0)

        return atomTypeCounts

    def coulombE(self, alldistpairs, recatomInfor, ligatomInfor,
                 maxCutoff, dielectric,
                 ):
        ## calculate eletrostatic interactions
        ## V = f (q1 * q2 ) / (epsilon * r12 )
        ## f = 1 / (4 * pi * epsilon) = 138.935 485 kJ * nm / (mol * e ^2)
        f = 138.935485

        ## all pair of atoms
        atoms = alldistpairs.keys()

        recRes = []

        backboneAtoms = ["C", "N", "O", "CA"]

        for atom in atoms:
            res = atom.split("+")[0].split("_")[1] + "_" + \
                  atom.split("+")[0].split("_")[2] + "_" + \
                  atom.split("+")[0].split("_")[3]
            if res not in recRes:
                recRes.append(res)
                # get the full residue_seq list

        backElectroEner, sideElectroEner = {}, {}  # format is : {"ARE175": -1.20 }
        for res in recRes:
            backElectroEner[res] = 0.0
            sideElectroEner[res] = 0.0

        for atom in atoms:
            if alldistpairs[atom] <= maxCutoff:

                recline = recatomInfor[atom.split("+")[0]]
                ligline = ligatomInfor[atom.split("+")[1]]
                #residueID = atom.split("+")[0].split("_")[1] + "_" + atom.split("+")[0].split("_")[2]
                res = atom.split("+")[0].split("_")[1] + "_" + \
                      atom.split("+")[0].split("_")[2] + "_" + \
                      atom.split("+")[0].split("_")[3]

                #q1 = RecCharges[ recline[16:20].strip() + "_" + recline[12:16].strip()]
                #q2 = LigCharges[ ligline[16:20].strip() + "_" + ligline[12:16].strip()]
                q1 = float(recline.split()[-2])
                q2 = float(ligline.split()[-2])

                if atom.split("+")[0].split("_")[0] in backboneAtoms :
                    # from angstrom to nanometer
                    backElectroEner[res] += f * q1 * q2 /(alldistpairs[atom] * 0.1 * dielectric)
                else :
                    # from angstrom to nanometer
                    sideElectroEner[res] += f * q1 * q2 /(alldistpairs[atom] * 0.1 * dielectric)

        return  recRes, backElectroEner, sideElectroEner

    def extractFeatures(self, inputfile="input.txt",
                        outputfile="output.dat",
                        vdwCountCutoff=6.0,
                        vdwEnerCutoff=12.0,
                        colEnerCutoff=12.0,
                        dielec=4.0,

                        ) :
        """
        extract necessary short range interaction features
        :param inputfile:
        :param outputfile:
        :param vdwCountCutoff:
        :param vdwEnerCutoff:
        :param colEnerCutoff:
        :param dielec:
        :return:
        """

        # the Van der Waals interaction parameters
        # the format is: { "C" : [0.339, 0.359824]; "DU" : [0.339, 0.359] }
        vdwParams = self.getVdWParams(self.atomtype)
        pdbfileLig = {}
        lines = open(inputfile)
        for s in lines :
            if s[0] != "#" :
                pdbfileLig[s.split()[0]] = s.split()[1]
        lines.close()

        allcases = []
        colnames = []
        for i in range(len(pdbfileLig.keys())) :
            print("Progress: no. %d file, out of %d"%(i, len(pdbfileLig.keys())))

            pdbfilename = sorted(pdbfileLig.keys())[i]

            receptorXYZ, ligandXYZ, recatomDetailInfor, ligatomDetailInfor = \
                self.getAtomInfor(pdbfilename, pdbfileLig[pdbfilename])
            alldistpairs = self.atomDistMatrix(receptorXYZ, ligandXYZ)
            # containing all the atom-atom distances data

            # counting contacts of atom pairs. Backbone and SideChain are seperated.
            # this function pass the test
            reslist2, backcount, sidecount = self.residueCounts(alldistpairs,
                                                                recatomDetailInfor,
                                                                distanceCutoff = vdwCountCutoff)
            # calculate Van der Waals contributions. Backbone and SideChain are seperated
            # this function works fine
            reslist2, backvan, sidevan = self.resVdWContribution(alldistpairs,
                                                                 recatomDetailInfor,
                                                                 ligatomDetailInfor,
                                                                 vdwParams ,
                                                                 maxCutoff = vdwEnerCutoff)

            # for all the atomtypes combinations, what are the contacts counts given a cutoff as 6.0 angstrom?
            # this function passes the test and works fine
            atomTypeCounts = self.contactsAtomtype(alldistpairs,
                                                   recatomDetailInfor,
                                                   ligatomDetailInfor,
                                                   vdwParams,
                                                   distanceCutoff = vdwCountCutoff )

            reslist2, backcol, sidecol = self.coulombE(alldistpairs,
                                                       recatomDetailInfor,
                                                       ligatomDetailInfor,
                                                       maxCutoff=colEnerCutoff,
                                                       dielectric=dielec,
                                                       )

            case = []

            # now, print the information into files
            if i == 0:
                reslist = sorted(reslist2)
                atomCombine = sorted(atomTypeCounts.keys())

                for cn in ["countback", "countside", "vdwback", "vdwside", "columbback", "columbside"] :
                    colnames += [cn+"_"+ x for x in reslist]

                colnames += [ "atomTCount" + "_" + x for x in atomCombine ]

            # countback and countside
            case += [ backcount[x] for x in reslist ]
            case += [ sidecount[x] for x in reslist ]

            # vdwback and vdwside
            case += [ backvan[x] for x in reslist ]
            case += [ sidevan[x] for x in reslist ]

            # columbback and cloumbback
            case += [ backcol[x] for x in reslist ]
            case += [ sidecol[x] for x in reslist ]

            # atomTcount
            case += [ atomTypeCounts[x] for x in atomCombine ]

            allcases.append(case)

        allcases = np.asarray(allcases)

        np.savetxt(outputfile, allcases, fmt="%.3f", delimiter=",",
                   header=",".join(colnames)
                   )

        return 1

if __name__ == "__main__" :
    d = '''
    Extract the binding features between receptor and ligand

    format of the input file:
    ------------------ input --------------
    pdbcomplex_1.pdb  LIG
    pdbcomplex_2.pdb  WAT
    pdbcomplex_3.pdb  STR
    ...
    ------------------ input --------------
    '''
    parser = argparse.ArgumentParser(description=d,formatter_class=RawTextHelpFormatter)
    parser.add_argument("-inp", default="input.dat",type=str,
                        help="The input file, containing names of proteins and their ligands as well.")
    parser.add_argument("-out",default="output.dat",type=str,
                        help="The output binding features file name.")
    parser.add_argument("-ac",type=int,default=1,
                        help="Atom counts between receptor and lignads.")
    parser.add_argument("-tc",type=int,default=1,
                        help="Counting of contacts by atom types. ")
    parser.add_argument("-vdwE", default=1, type=int,
                        help="Whether van del waals calculate energy per residue")
    parser.add_argument("-parm",default="AtomType.dat",type=str,
                        help="The parameters for VDW of each type of atoms.")

    args = parser.parse_args()

    binding = BindingFeature()

    binding.extractFeatures(args.inp, args.out)