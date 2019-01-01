# -*- coding: utf-8 -*-

from collections import defaultdict
import math
import argparse
from argparse import RawTextHelpFormatter
import sys, os
import numpy as np
import dockml.index as index
import dockml.pdbIO as pio
import mdanaly, dockml

class GridBasedFeature :

    def __init__(self, ligandPDB, receptorPDB, gridsize=1.0):
        '''

        :param gridsize: float, grid size, unit angstrom
        '''

        self.gridsize = gridsize
        self.ligpdb   = ligandPDB
        self.recpdb   = receptorPDB

        bf = BindingFeature()
        self.vdwparm = bf.getVdWParams()
        self.eleparm = bf.getElementParams()

        self.atominfor = pio.parsePDB().atomInformation(self.ligpdb)

        try :
            with open(self.ligpdb) as lines :
                self.ligndx = [ x.split()[1] for x in lines if "ATOM" in x ]
        except :
            self.ligndx = []

    def pocketFromRes(self, receptorPDB, reslist, chain='A'):
        '''
        get a list of residues' index, return their coordinates (angstrom)
        :param receptorPDB:
        :param reslist: list of int
        :param chain:
        :return:
        '''

        pndx = []
        for res in reslist:
            ndx = index.PdbIndex(reference=self.recpdb, chain=[chain, ],
                                 atomtype="all", resSeq=[res,])
            ndx.prepare_selection().res_index()
            pndx += ndx.atomndx

        pocketResCoord = pio.coordinatesPDB().getAtomCrdByNdx(singleFramePDB=self.recpdb,
                                                              atomNdx=pndx)
        return pocketResCoord

    def createCubicPocket(self, coordinates, extension=10.0):
        '''
        create the 6 element vector to describe a cubic pocket
        :param coordinates: ndarray, M*3 dim
        :param extension: float, extend the box and enlarge it
        :return: ([x, y, z], [x, y, z] )
        '''

        coordinates = np.asarray(coordinates)

        maxs = coordinates.max(axis=0) + extension
        mins = coordinates.min(axis=0) - extension

        return (maxs, mins)

    def generateGrids(self, boundaries, gridfile="GRID"):
        '''
        the grid ndarray, generated from the boundaries in x y z dimension
        :param boundaries: list of list, [ [x, y, z], [x, y, z]]
        :param gridfile: str, output grid file format
        :param gridsize: float, grid bin size
        :return: ndarray of floats, N * 3 dimmension
        '''

        upbound, lowbound = boundaries[0], boundaries[1]

        x_range = np.arange(lowbound[0], upbound[0], self.gridsize)
        y_range = np.arange(lowbound[1], upbound[1], self.gridsize)
        z_range = np.arange(lowbound[2], upbound[2], self.gridsize)

        xsize, ysize, zsize = x_range.shape[0], y_range.shape[0], z_range.shape[0]

        grid_x = np.repeat(x_range, ysize*zsize)

        grid_y = np.tile(np.repeat(y_range, zsize), xsize)

        grid_z = np.tile(np.tile(z_range, ysize), xsize)

        grid_center = np.concatenate((np.array([grid_x]),
                                      np.array([grid_y]),
                                      np.array([grid_z]),
                                      ), axis=0
                                     ).T

        np.savetxt(gridfile, grid_center, delimiter=" ", fmt="%8.3f")

        return grid_center

    def ligCoords(self, chain="B"):
        '''
        from ligandpdb file get ligand atoms coordinates
        :param chain: str
        :param resndx: list of int
        :return: ndarray, N*3 dimension
        '''

        coords = defaultdict(list)

        for ndx in self.ligndx :
            coords[ndx] = pio.coordinatesPDB().getAtomCrdByNdx(singleFramePDB=self.ligpdb,
                                                               atomNdx=[ndx,])[0]

        return coords

    def atomGridBoxId(self, grid_center, atomcrd, count_cutoff=0.01):
        '''
        calculate whether a point in grid bind centered around grid_center
        return the grid_id and counts
        :param grid_center:
        :param atomcrd:
        :param count_cutoff:
        :return:
        '''

        try :
            num_bin = grid_center.shape[0]
        except TypeError :
            num_bin = np.asarray(grid_center).shape[0]

        for num in range(num_bin) :
            grid = grid_center[num]
            count = self.pointInGrid(atomcrd, grid)
            if count > count_cutoff :
                print(num, count)
                return num, count
            else :
                pass

        return -1, 0.0

    def distToGrid(self, grid_point, physic_point):

        distance = mdanaly.ContactMap(self.ligpdb).atomDistance(list(grid_point),
                                                                list(physic_point),
                                                                sqrt=True
                                                                )

        return distance

    def pointInGrid(self, point, grid, cutoff=4.0):
        '''
        determine a point in a grid box or not
        :param point: float, a list of 3 items
        :param grid: float, a list of 3 items, the grid center
        :param gridsize: float, the grid size
        :return: float, scaled contact number
        '''

        '''if point[0] > grid[0] and point[0] <= grid[1] \
                        and point[1] > grid[2] and point[1] <= grid[3] \
                        and point[2] > grid[4] and point[2] <= grid[5] :
                    return True
                else :
                    pass

                return False'''

        dist = self.distToGrid(point, grid)
        dist_scaled = dockml.BasicAlgorithm().switchFuction(dist, d0=cutoff*2)

        return dist_scaled

    def gridAtomMap(self, grid_table, ligcoords):
        '''
        determine in each grid bin, whether ligand atom exists there
        :param grid_table: ndarray, N*3 dim
        :param ligcoords: defaultdict(list),{ key, dimension N*3},
                          the coordinates of ligand atoms
        :return: { grid_id: (atomndx, contactnumber) }
        '''

        gridMapper = defaultdict(list)

        for k in range(grid_table.shape[0]) :
            grid = grid_table[k]
            gridMapper[k] = [ (x, self.pointInGrid(ligcoords[x], grid)) for x in self.ligndx ]

            '''for ndx in self.ligndx :
                
                crd = ligcoords[ndx]
                #print(grid, crd)
                count = self.pointInGrid(crd, grid)
                #count = 0.0

                if count > count_cutoff :
                    #print("found", k, i, count)
                    try :
                        gridMapper[k].append((ndx, count))
                    except KeyError :
                        pass
            '''
        return gridMapper

    def gridBinProperty(self, gridMapper, properties, count_cutoff=0.2):
        '''
        get atomic based properties from the gridmapper
        :param gridMapper: defaultdict(list)
        :return: defaultdict(list), { grid_id, [dim=5]
        '''

        grid_features = defaultdict(list)
        grid_ids = sorted(gridMapper.keys())

        #properties = self.atomProperties(self.ligndx)

        for id in grid_ids :
            # mass vdw_epsilon, vdw_sigma, charge, negativity
            #p = np.array([ 0.0, 0.0, 0.0, 0.0, 0])
            grid_features[id] = [[ 0.0, 0.0, 0.0, 0.0, 0], ]

            grid_features[id] = [ list(np.array(properties[x[0]]) * x[1])
                                  for x in gridMapper[id]
                                  if (gridMapper[id] and len(x) == 2 and x[1] >= count_cutoff)
                                  ]

            grid_features[id] = np.sum(np.array(grid_features[id]), axis=0)
            #grid_features[id] = (list(p))

        return grid_features

    def getLigPartialCharges(self):
        '''
        get atom partial charges
        :return:
        '''

        if self.ligpdb[-5:] == "pdbqt" :
            charges = {}
            with open(self.ligpdb) as lines :
                for s in [ x for x in lines if "ATOM" in x ] :
                    charges[s.split()[1]] = float(s.split()[-2])
            return charges
        else :
            return {}

    def atomProperties(self, ligatomndx):
        '''
        get the properties of a atom
        mass, vdw_e, vdw_s, charge, negativity
        :param ligatomndx:
        :return:
        '''

        properties = defaultdict(list)

        for atomndx in ligatomndx :
            # mass vdw_epsilon, vdw_sigma, charge, negativity
            p = [0.0, 0.0, 0.0, 0.0, 0.0]

            atom_parms = self.atominfor[atomndx]

            ele = atom_parms[7]
            if ele == "A" :
                ele = "C"

            try :
                [e, s] = self.vdwparm[ele]
            except :
                e, s = 0.0, 0.0

            try :
                [eleid, negat] = self.eleparm[ele]
            except :
                eleid, negat = 0.0, 0.0

            try :
                charges = self.getLigPartialCharges()
                if atomndx in charges.keys() :
                    p[3] = charges[atomndx]
            except FileNotFoundError :
                pass

            p[1], p[2] = e, s
            p[0], p[4] = eleid, negat

            properties[atomndx] = p

        return properties

class BindingFeature :

    def __init__(self, vdwCutoff=6., vdwEnerCutoff=12.,
                 colEnerCutoff=12., contCutoff=12., disCutoff=3.5,
                 pdbqt=True,
                 ):
        if os.path.exists("AtomType.dat") and os.path.exists("elementNegativity.dat"):
            self.atomtype = "AtomType.dat"
            self.elenegat = "elementNegativity.dat"
        else :
            PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
            DEFINITIONS_ROOT = os.path.join(PROJECT_ROOT, '../data/AtomType.dat')
            self.atomtype = DEFINITIONS_ROOT
            self.elenegat = os.path.join(PROJECT_ROOT, '../data/elementNegativity.dat')

        # parameters
        self.vdwCountCutoff = vdwCutoff
        self.vdwEnerCutoff = vdwEnerCutoff
        self.colEnerCutoff = colEnerCutoff
        self.contactCutoff = contCutoff
        self.distCutoff    = disCutoff

        self.backboneAtoms = ["C", "N", "O", "CA"]
        self.pdbqt = pdbqt

        # attributes
        self.vdwParameters_ = self.getVdWParams()

    def getVdWParams(self) :
        '''
        Get vdw parameters, mainly sigma and epsilon term in energy function
        :inputFile: str, a library file containing atomic vdw parameters
        :return:
        '''
        # obtain vdw radius for different atoms
        vdwParams = defaultdict(list)
        with open(self.atomtype) as lines :
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

    def getElementParams(self):
        '''
        get element negativity and atomic mass
        :return: dict, { element: [atomic mass, negativity, ] }
        '''

        elemParams = defaultdict(list)
        with open(self.elenegat) as lines :
            for s in [ x for x in lines if (x[0] not in [";", "#"] and len(x.split()) > 3 )] :
                elemParams[ s.split()[1] ] = (int(s.split()[0]), float(s.split()[3]))
        return elemParams

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
        :return: recXYZ, ligXYZ, recatomLine, ligatomLine
        '''
        recXYZ = defaultdict(list)
        recatomLine = {}
        ligatomLine = {}
        ligXYZ = defaultdict(list)

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
                    recXYZ[ atomID ] = self.getXYZCoord(s)
                    recatomLine[atomID] = s

                ## ligand information
                elif ligCode ==  s[17:20].strip() :
                    ligXYZ[ atomID ] = self.getXYZCoord(s)
                    ligatomLine[atomID] = s

        return recXYZ, ligXYZ, recatomLine, ligatomLine

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
        :param receptorXYZ: list 2D, dimension N*3, a list of coordinates of atoms
        :param ligandXYZ: list 2D, dimension N*3, a list of coordinates of atoms
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

    def residueCounts(self, alldistpairs, atomDetailInfor) :
        """
        input all the distances pairs
        return the number of contacts
        alldistpairs: the porteinID-LigID-Distance inforamtion
        atomdetailInfor : the atomID-PDBline information
        distance unit is angstrom
        :param alldistpairs: dictionary, all distances
        :param atomDetailInfor: dictionary, detail information
        :return:
        """

        atoms = alldistpairs.keys()
        recRes =[]


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
            if distance <= self.contactCutoff :
                line = atomDetailInfor[atom.split("+")[0]]

                ### add more chain ID, in case you have several different chains
                ### better, you may use 1 2 3 4 as your chain id, here choose original chain id
                ### residueID  residueName_ResidueIndex_Chain
                residueID = atom.split("+")[0].split("_")[1] + "_" + \
                            atom.split("+")[0].split("_")[2] + "_" + \
                            atom.split("+")[0].split("_")[3]
                #residueID = line[17:20].strip() + "_" + line[22:26].strip()
                ## count contacts in backbones
                if line.split()[2] in self.backboneAtoms :
                    backCount[residueID] += self.switchFuction(distance, self.distCutoff*2.0)

                ## count contacts in sidechains
                elif line.split()[2] not in self.backboneAtoms :
                    sideCount[residueID] += self.switchFuction(distance, self.distCutoff*2.0)

        return recRes, backCount, sideCount

    def atomicVdWEnergy(self, atomtype1, atomtype2, distance) :
        """
        Calculate atomic level Van der Waals interaction energy
        VdW potential calculation, using combination rule 2
        vdwParam is a dictionary list, eg. { "O": [ 0.22617E-02, 0.74158E-0.6 ]; }
        unit of the parameters: nanometer
        :param atomtype1:
        :param atomtype2: dictionary, parameters of atom type and sigma epsilon
        :param distance: float
        :return: float, short range vdw energy
        """

        Sigma_i   = self.vdwParameters_[atomtype1][0]
        Epsilon_i = self.vdwParameters_[atomtype1][1]
        Sigma_j   = self.vdwParameters_[atomtype2][0]
        Epsilon_j = self.vdwParameters_[atomtype2][1]

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

        return C12 / d12 - C6 / d6

    def resVdWContribution(self, alldistpairs, recatomInfor, ligatomInfor,
                           vdwParams, repulsionMax=10.0 ) :
        """
        accumulated vdw for side-chain or backbone
        alldistpairs: the porteinID-LigID-Distance inforamtion
        atomdetailInfor : the atomID-PDBline information
        :param alldistpairs: dict, all distances
        :param recatomInfor:
        :param ligatomInfor: dict, the ligand-ligand_line information
        :param vdwParams:  dict, all vdw parameters
        :param pdbqt: bool, whether it is a pdbqt file
        :return:
        """

        atoms = alldistpairs.keys()
        recRes =[]

        #backboneAtoms = ["C", "N", "O", "CA"]

        for atom in atoms:
            res = atom.split("+")[0].split("_")[1] + "_" + \
                  atom.split("+")[0].split("_")[2] + "_" + \
                  atom.split("+")[0].split("_")[3]
            if res not in recRes:
                recRes.append(res)  # get the full residue_seq_chain list

        backVdWEner = {}  # format is : {"ARE_175_A": -1.20 }
        sideVdWEner = {}  # format is : {"ARE_175_A": -1.20 }
        for res in recRes:
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
            if self.pdbqt :
                atomtype1 = recline.split()[-1]
                atomtype2 = ligline.split()[-1]
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
                elif atomtype1 in ["L", "CL"]:
                    atomtype1 = "Cl"
                elif atomtype1 in ["B", "BR"]:
                    atomtype1 = "Br"
                else :
                    atomtype1 = "DU"
            if atomtype2 not in vdwParams.keys() :
                if atomtype2 == " N1+":
                    atomtype2 = "N"
                elif atomtype2 in ["L", "CL", "Cl"]:
                    atomtype1 = "Cl"
                elif atomtype2 in ["B", "BR", "R", "Br"]:
                    atomtype1 = "Br"
                else :
                    atomtype2 = "DU"

            # calculate L J potential per residue
            if alldistpairs[atom] <= self.vdwEnerCutoff :
                # transfer angstrom to nanometer by multiplying 0.1
                energy = self.atomicVdWEnergy(atomtype1, atomtype2, alldistpairs[atom] * 0.1, vdwParams)
            else:
                energy = 0.0
                #energy = self.atomicVdWEnergy(atomtype1, atomtype2, maxCutoff, vdwParams)

            if recline.split()[2] in self.backboneAtoms:
                backVdWEner[residueID] += energy
            else:
                sideVdWEner[residueID] += energy

        for key in backVdWEner.keys():
            if backVdWEner[key] > repulsionMax:
                backVdWEner[key] = repulsionMax
        for key in sideVdWEner.keys():
            if sideVdWEner[key] > repulsionMax:
                sideVdWEner[key] = repulsionMax

        return recRes, backVdWEner, sideVdWEner

    def contactsAtomtype(self, alldistpairs, recatomInfor, ligatomInfor,):

        """
        # vdw counting based on different atom types
        # alldistpairs: the porteinID-LigID-Distance inforamtion
        # atomdetailInfor : the atomID-PDBline information
        # vdwParams is the C6 and C12 terms for different atom types
        :param alldistpairs:
        :param recatomInfor:
        :param ligatomInfor:
        :param pdbqt: bool, whether it is a pdbqt file
        :return:
        """
        atoms = alldistpairs.keys()
        atomTypeCounts = {}
        for atomT1 in self.vdwParameters_.keys() :
            for atomT2 in self.vdwParameters_.keys() :
                atomTypeCounts[ atomT1 + "_" + atomT2 ] = 0.0

        for atom in atoms :
            distance = alldistpairs[atom]
            if distance <= self.distCutoff * 2 :
                recline = recatomInfor[atom.split("+")[0]]
                ligline = ligatomInfor[atom.split("+")[1]]

                # unrecognized atomtype occur, "dummy" atomtype is used
                # create all the combinations of the atomtypes known
                if self.pdbqt:
                    ligAtom = ligline.split()[-1]
                else :
                    #ligAtom = ligline.split()[-1]
                    ligAtom = ligline[13]

                if ligAtom not in self.vdwParameters_.keys():
                    if ligAtom == "N1+":
                        ligAtom = "N"

                    elif ligAtom in ["L", "CL", "Cl"] :
                        ligAtom = "Cl"
                    elif ligAtom in ["B", "BR", "Br", "R"] :
                        ligAtom = "Br"
                    else:
                        ligAtom = "DU"
                        print("DU {}".format(ligAtom), ligline)


                if self.pdbqt :
                    recAtom = recline.split()[-1]
                else :
                    #recAtom = recline.split()[-1]
                    recAtom = recline[13]

                if recAtom not in self.vdwParameters_.keys() :
                    if recAtom == "N1+" :
                        recAtom = "N"
                    elif recAtom == "CL" :
                        recAtom = "Cl"
                    else :
                        recAtom = "DU"

                atomTypeCombination = recAtom + "_" + ligAtom
                if atomTypeCombination not in atomTypeCounts.keys():
                    # apply a switch function to smooth the transition
                    atomTypeCounts[atomTypeCombination] = self.switchFuction(distance, self.distCutoff*2.0)
                else:
                    atomTypeCounts[atomTypeCombination] += self.switchFuction(distance, self.distCutoff*2.0)

        return atomTypeCounts

    def coulombE(self, alldistpairs, recatomInfor, ligatomInfor,
                 dielectric=1.0,
                 ):
        """
        calculate residue-ligand interaction coulomb energies
        :param alldistpairs: atom-distance information
        :param recatomInfor: dict, atom_index: atom_line
        :param ligatomInfor: dict, atom_index: atom_line
        :param maxCutoff: float, cutoff of the columbic eneriges
        :param dielectric: float, dielectric constant
        :return: tuple, (residues, backbone energies, sidechain energies)
        """
        ## calculate eletrostatic interactions
        ## V = f (q1 * q2 ) / (epsilon * r12 )
        ## f = 1 / (4 * pi * epsilon) = 138.935 485 kJ * nm / (mol * e ^2)
        f = 138.935485

        ## all pair of atoms
        atoms = alldistpairs.keys()
        recRes = []

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
            if alldistpairs[atom] <= self.vdwEnerCutoff:  #maxCutoff:

                recline = recatomInfor[atom.split("+")[0]]
                ligline = ligatomInfor[atom.split("+")[1]]
                #residueID = atom.split("+")[0].split("_")[1] + "_" + atom.split("+")[0].split("_")[2]
                res = atom.split("+")[0].split("_")[1] + "_" + \
                      atom.split("+")[0].split("_")[2] + "_" + \
                      atom.split("+")[0].split("_")[3]

                #q1 = RecCharges[ recline[16:20].strip() + "_" + recline[12:16].strip()]
                #q2 = LigCharges[ ligline[16:20].strip() + "_" + ligline[12:16].strip()]
                if self.pdbqt :
                    q1 = float(recline.split()[-2])
                    q2 = float(ligline.split()[-2])
                else :
                    q1, q2 = 0., 0.

                if recline.split()[2] in self.backboneAtoms :
                    #if atom.split("+")[0].split("_")[0] in backboneAtoms :
                    # from angstrom to nanometer
                    backElectroEner[res] += f * q1 * q2 /(alldistpairs[atom] * 0.1 * dielectric)
                else :
                    # from angstrom to nanometer
                    sideElectroEner[res] += f * q1 * q2 /(alldistpairs[atom] * 0.1 * dielectric)

        return recRes, backElectroEner, sideElectroEner

    def extractFeatures(self, inputfile="input.txt",
                        outputfile="output.dat",
                        dielec=4.0,
                        feature_terms = [ True, True, True, True],
                        ) :
        """
        extract necessary short range interaction features
        :param inputfile: str, file containing two cols, filename and ligand code
        :param outputfile: str, file containing all energy terms
        :param vdwCountCutoff: float
        :param vdwEnerCutoff: float
        :param colEnerCutoff: float
        :param dielec:
        :return:
        """

        error_log = open("./error.log", 'w')
        # the Van der Waals interaction parameters,
        # the format is: { "C" : [0.339, 0.359824]; "DU" : [0.339, 0.359] }
        vdwParams = self.getVdWParams()

        with open(inputfile) as lines :
            fnames_ligs = [ (x.split()[0], x.split()[1]) for x in lines if "#" not in x ]
            pdbfileLig = dict(fnames_ligs)

        allcases = []
        colnames = []
        num_cols = 0
        case_names = []

        for i in range(len(pdbfileLig.keys())) :
            print("Progress: no. %d file, out of %d"%(i, len(pdbfileLig.keys())))

            # obtain pro-lig complex file name
            pdbfilename = fnames_ligs[i][0]

            # get atom information
            receptorXYZ, ligandXYZ, recatomDetailInfor, ligatomDetailInfor = \
                self.getAtomInfor(input=pdbfilename, ligCode=fnames_ligs[i][1])

            alldistpairs = self.atomDistMatrix(receptorXYZ, ligandXYZ)
            # containing all the atom-atom distances data

            # counting contacts of atom pairs. Backbone and SideChain are seperated.
            # this function pass the test
            reslist1, backcount, sidecount = self.residueCounts(alldistpairs,
                                                                recatomDetailInfor,
                                                                )
            # calculate Van der Waals contributions. Backbone and SideChain are seperated
            # this function works fine
            reslist2, backvan, sidevan = self.resVdWContribution(alldistpairs,
                                                                 recatomDetailInfor,
                                                                 ligatomDetailInfor,
                                                                 vdwParams ,
                                                                 pdbqt=False,
                                                                 )

            # for all the atomtypes combinations, what are the contacts counts given a cutoff as 6.0 angstrom?
            # this function passes the test and works fine
            atomTypeCounts = self.contactsAtomtype(alldistpairs,
                                                   recatomDetailInfor,
                                                   ligatomDetailInfor,
                                                   vdwParams,
                                                    )

            reslist3, backcol, sidecol = self.coulombE(alldistpairs,
                                                       recatomDetailInfor,
                                                       ligatomDetailInfor,
                                                       dielectric=dielec,
                                                       )

            case = []

            # now, print the information into files
            if i == 0:
                reslist = sorted(list(set.intersection(set(reslist1), set(reslist2), set(reslist3))))
                atomCombine = sorted(atomTypeCounts.keys())

                for cn in ["countback", "countside", "vdwback", "vdwside", "columbback", "columbside"]:
                    colnames += [cn+"_" + x for x in reslist]

                colnames += ["atomTCount" + "_" + x for x in atomCombine]

                num_cols = len(colnames)
                # case_names.append(pdbfilename)

            try :

                # countback and countside
                case += [backcount[x] for x in reslist]
                case += [sidecount[x] for x in reslist]

                # vdwback and vdwside
                case += [backvan[x] for x in reslist]
                case += [sidevan[x] for x in reslist]

                # columbback and cloumbback
                case += [backcol[x] for x in reslist]
                case += [sidecol[x] for x in reslist]

                # atomTcount
                case += [atomTypeCounts[x] for x in atomCombine]
                case_names.append(pdbfilename)
                allcases.append(case)

            except :
                if num_cols == 0 :
                    num_cols = len(allcases[-1])
                #case = [ 0.0 ] * num_cols
                error_log.write("%s errors here \n"%pdbfilename)
                print("Errors: %s "%pdbfilename)

        tofile = open("ligand_features_name.dat", "w")
        tofile.write(",".join(case_names)+"\n")

        error_log.write(",".join(case_names)+"\n")
        allcases = np.asarray(allcases)

        # output the features
        np.savetxt(outputfile, allcases, fmt="%.3f", delimiter=",",
                   header=",".join(colnames)
                   )

        error_log.close()

        return 1

class LigandFingerPrints(BindingFeature) :

    def elementCount(self, ligand):
        from dockml import pdbIO

        elements = self.getVdWParams().keys()
        atominfor = pdbIO.parsePDB(ligand).atomInformation(ligand)

        elem_count = dict(zip(elements, np.zeros(len(elements))))

        for atom in atominfor.keys():

            if atominfor[atom][7] not in elements:
                elem_count["DU"] += 1
            else :
                elem_count[atominfor[atom][7]] += 1

        return elem_count

def main() :
    d = '''
    Extract the binding features between receptor and ligand
    
    Usage: 
    features.py -inp input.dat -out features.dat 

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

    '''
    parser.add_argument("-ac",type=bool,default=True,
                        help="Atom counts between receptor and lignads.")
    parser.add_argument("-tc",type=bool,default=True,
                        help="Counting of contacts by atom types. ")
    parser.add_argument("-vdwE", default=True, type=bool,
                        help="Whether van del waals calculate energy per residue")
    parser.add_argument("-columb", type=bool, default=True,
                        help="Whether include columbic interactions. Default is True")
    parser.add_argument("-parm",default="AtomType.dat",type=str,
                        help="The parameters for VDW of each type of atoms.") '''

    args = parser.parse_args()

    binding = BindingFeature()

    binding.extractFeatures(args.inp, args.out, )