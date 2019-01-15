#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import subprocess as sp
from glob import glob
from collections import *
import linecache
import glob
from .gentop import GenerateTop

# import modeller for loop refinement
try:
    from modeller import *
    from modeller.automodel import *
    MODELLER_EXIST = True
except ImportError :
    print("Warning: Modeller is not well imported, some features may not work. ")
    MODELLER_EXIST = False


class MolDocking(object):

    def __init__(self):
        pass

    def atomicChargesITP(self, ff_file="protein.itp"):
        '''
        Read the atomic charge from the itp file
        :param ff_file:
        :return:
        '''

        condition = False
        atomic_charges = {}

        with open(ff_file) as lines :
            for s in lines :
                if "[ atoms ]" in s :
                    condition = True
                if "[ bonds ]" in s :
                    condition = False

                if len(s.split()) and s[0] != ";" and condition :
                    # resname + "_" + atomname
                    id = s.split()[3] + "_" + s.split()[4]
                    if id not in atomic_charges :
                        atomic_charges[id] = float(s.split()[6])

        return atomic_charges

    def atomicChargesLig(self, obabelexe='obabel', ff_file="ligand.itp", netcharges=None, ):
        """

        :param obabelexe:
        :param ff_file:
        :param netcharges: if none, the total net charges will be deduced from the mol2 or pdbqt file
                            otherwise, default value 0 will be given.
        :return:
        """
        atomic_charges = {}
        if ".itp" == ff_file[-4:] :
            return self.atomicChargesITP(ff_file)

        elif ".mol2" == ff_file[-5:] :
            # generate a shell file for antechamber
            gentop = GenerateTop()
            gentop.runAntechamber(ff_file, netcharges)

            # run antechamber
            job = sp.Popen("sh antechamber.sh %s %s" % (ff_file, ff_file.split(".")[0]), shell=True)
            job.communicate()

            # run gentop
            if not os.path.exists(ff_file.split(".")[0]+".pdb") :
                gentop.runObabel(obabelexe, ff_file, ff_file.split(".")[0]+".pdb")

            gentop.gmxTopBuilder(ff_file.split(".")[0]+".pdb", ff_file.split(".")[0],
                                 frcmodFile="frcmod."+ff_file.split(".")[0],
                                 prepfile="prep." + ff_file.split(".")[0], boxEdge=0,
                                 FField=["gaff",],
                                 )
            # you will have a *.top and *.itp file
            try :
                return self.atomicChargesITP(ff_file.split(".")[0]+".top")
            except :
                return self.atomicChargesITP(ff_file.split(".")[0] + ".itp")

        if ".pdb" == ff_file[-4:]  :
            gentop = GenerateTop()
            if not os.path.exists(ff_file.split(".")[0]+".mol2") :
                gentop.runObabel(obabelexe, ff_file, ff_file.split(".")[0]+".mol2")

            gentop.runAntechamber(ff_file.split(".")[0]+".mol2", netcharges)

            # run antechamber
            job = sp.Popen("sh antechamber.sh %s %s" % (ff_file.split(".")[0]+".pdb", ff_file.split(".")[0]),
                           shell=True)
            job.communicate()

            # gen top here
            gentop.gmxTopBuilder(ff_file, ff_file.split(".")[0],
                                 frcmodFile="frcmod." + ff_file.split(".")[0],
                                 prepfile="amberff.prep." + +ff_file.split(".")[0], boxEdge=0,
                                 FField=["gaff", ],
                                 )

            # return charges
            return self.atomicChargesITP(ff_file.split(".")[0] + ".itp")

    def getProperty(self, csvFile, ligandid,head=True, property="IC50", delimator=",", position=-1):
        '''
        From a mol2 file read the property of a molecule
        :param csvFile:
        :param ligandid:
        :param head:
        :param property:
        :param delimator:
        :param position:
        :return:
        '''
        fields  =[]
        prop = ''

        if head :
            fields = linecache.getline(csvFile, 1)

        with open(csvFile) as lines :
            for s in lines :
                if len(s.split()) and s.split()[0] == ligandid :
                    if len(fields) and not position :
                        prop = s.split(delimator)[fields.index(property)]
                    else :
                        prop = s.split(delimator)[position]
                else :
                    pass
        return prop

def runGenTop(AMBERHOME) :
    '''
    instant a gmxTopBuilder class, and input some parameters,
        to generate amber and gromacs topology and coordinates

    if a small ligand provide, if not prep and frcmod exist,
        then am1/bcc charge based on amber antechamber sqm method
        will be calculated. in this process, tleap will be used to
        get bonded and non-bonded parameters
    :parameter: None
    :return: None
    '''
    top = GenerateTop()
    args = top.arguments()

    try :
        os.environ["AMBERHOME"] = AMBERHOME
    except :
        os.system("export AMBERHOME=%s"%AMBERHOME)

    # define the input coordinates information, a pdb, or mol2, or a sequence
    # of amino acids
    if args.aaseq:
        structure = args.aaseq
    elif args.inp:
        structure = args.inp
    else:
        structure = ''
        print('Error: No input structure defined.\n Exit Now!')
        sys.exit(0)

    if not args.resp:
        # not calculate atomic charges
        top.gmxTopBuilder(args.frcmod, args.prep,
                          structure, args.out,
                          args.ion, args.nion,
                          args.bt, args.size,
                          args.ff)
    else:
        # prepare resp charges for small ligand with antechamber
        # give a netcharge: args.nc
        sh = top.runAntechamber(args.inp, args.nc)

        try:
            # set env for AMBERHOME
            os.environ["AMBERHOME"] = AMBERHOME
            job = sp.Popen("sh %s %s %s" % (sh, args.inp, args.out), shell=True)
            job.communicate()
        except :
            print("Antechamber generateing RESP and atomic parameters error! \n "
                  "Double check your input structure.")
            sys.exit(1)

        top.gmxTopBuilder("frcmod." + args.out,
                          "prep." + args.out,
                          structure, args.out,
                          args.ion, args.nion,
                          args.bt, args.size,
                          args.ff)

    if args.ion[0] == "X+" or args.bt == "NA":
        # parse a gromacs .top to a .itp file
        # print "ITP file created!"
        pass
        #top.removeMolInfor(args.out)

