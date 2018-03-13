#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess as sp
from glob import glob
from collections import *
import linecache
import glob
from .gentop import GenerateTop

from mdanaly import extract

# import modeller for loop refinement
try:
    from modeller import *
    from modeller.automodel import *
    MODELLER_EXIST = True
except ImportError :
    print("Warning: Modeller is not well imported, some features may not work. ")
    MODELLER_EXIST = False

class MolDocking :
    def __init__(self):
        pass

    def runVina(self, vinaexe, receptor, ligand,
                output='result.pdbqt', logfile='result.log',
                ncpu=1, exhaustiveness=32,
                center=[0,0,0], sizes=[40,40,40],
                no_modes = 20, en_range = 5, seed=-1,
                ):
        '''
        perform molecular docking using autodock vina

        :param vinaexe: executable vina binary file
        :param receptor: receptor name in pdbqt format
        :param ligand: ligand file name in pdbqt format
        :param output: ligand binding pose, in pdbqt format
        :param logfile: docking results log
        :param ncpu: number of cups
        :param exhaustiveness: how accurate
        :param center: binding pocket center
        :param sizes: size of the binding pocket
        :param no_modes: output number of modes
        :param en_range: energy range
        :param seed: random seed
        :return:
        '''
        try :
            job = sp.Popen('%s --receptor %s --ligand %s '
                           '--center_x %f --center_y %f --center_z %f '
                           '--size_x %f --size_y %f --size_z %f '
                           '--log %s --out %s --cpu %d '
                           '--exhaustiveness %d --num_modes %d '
                           '--energy_range %d --seed %d ' %
                           (
                               vinaexe, receptor, ligand,
                               center[0], center[1], center[2],
                               sizes[0], sizes[1], sizes[2],
                               logfile, output, ncpu,
                               exhaustiveness, no_modes,
                               en_range, seed),
                           shell= True
                           )
            job.communicate()
            job.terminate()
        except IOError :
            print("Docking molecule %s to %s using %s failed. \n"
                  "Check your input and logfile.")

            job = sp.check_output('%s --help'%vinaexe)
            print(job)

        return(1)

    def rankVinaResults(self, logfileList):
        '''
        obtain the binding energy score from vina log files

        :param logfileList:
        :return: a list of tuples, each key matches with a result list (top 3 results)
        '''
        vinaResults = defaultdict(list)

        for resultfile in logfileList :
            condition = -1
            vinaResults[resultfile] = []
            with open(resultfile) as lines :
                for s in lines :

                    if "Refining results ... done" in s :
                        condition += 1
                    elif "affinity" in s :
                        condition += 1
                    else :
                        pass

                    if condition :
                        if len(s.split()) and s.split()[0] in ['1','2','3'] :
                            vinaResults[resultfile].append(float(s.split()[1]))

        return(sorted(vinaResults.items(), key=lambda x: x[1]))

    def findOriginalLig(self, filename, ligid, type="mol2"):
        if type == filename.split(".")[-1] and os.path.exists(filename) :
            with open(filename) as lines :

                ligcontent = []
                condition = 0
                for s in lines:
                    '''in a ligand content, condition = True '''
                    if "<TRIPOS>MOELECULES" in s:
                        condition += 1

                    elif condition == 1 :
                        if ligid == s.split()[0]:
                            condition += 1
                            ligcontent.append("<TRIPOS>MOELECULES \n")
                            ligcontent.append(ligid + "  \n")
                        else :
                            condition = 0

                    elif condition == 2:
                        ligcontent.append(s)
                        condition = 0

                    if condition >= 3 :
                        break

            return ligcontent
        else :
            print("Error! Mol2 file or mol2 type error!")
            return []

    def findDockedLig(self, filepattern, ligid, filenum=1):

        filelist = glob.glob(filepattern)
        thefiles = []

        ''' going to find the gold results file based on ligand id '''
        for filename in filelist:
            with open(filename) as lines:

                for s in lines:
                    if len(s.split()) > 0 and s.split()[0] in ligid.strip('./'):
                        thefiles.append(filename)
                        break

            if len(thefiles) > filenum:
                break

        return thefiles

    def getLigandID(self, ligidinfor, pathname):
        """
        Determine whether the input is a file or a list
        :param ligidinfor:
        :param pathname:
        :return:
        """
        if len(ligidinfor) > 1:
            '''multiple ligand id directly'''
            ligandlist = ligidinfor
        else:
            ''' supposing the result ligands ids in file, each record per line '''
            ligandlist = []
            if ligidinfor in glob.glob(pathname):
                lines = open(ligidinfor[0])
                # ligandlist = []
                for s in lines:
                    if len(s.split()) > 0 and "#" not in s:
                        ligandlist.append(s.split()[0])
                lines.close()
            else:
                ligandlist.append(ligidinfor[0])

        return ligandlist

    def getLigand(self, ligidfile, topNum, order):

        # ligidfile is in relative path
        if ligidfile in os.listdir("./"):
            # this file exist and the format is similar to bestranking lst file
            # docking score in descending order
            linecount = 0
            ligand = OrderedDict()

            if order == "d" or order[0] == "d":
                linecount = 1
                while linecount <= topNum:
                    s = linecache.getline(ligidfile, linecount)[:-1]  # remove the "\n"
                    ligand[s.split()[-1].strip("\'")] = s.split()[-2].strip("\'")
                    # increase line number to goto next line
                    linecount += 1
            elif order == "a" or order[0] == "a":
                lines = open(ligidfile)
                nl = len(lines.readlines())
                linecount = nl
                if nl < topNum:
                    bot = 0
                else:
                    bot = nl - topNum
                while linecount > bot:
                    s = linecache.getline(ligidfile, linecount)[:-1]  # remove the "\n"
                    ligand[s.split()[-1].strip("\'")] = s.split()[-2].strip("\'")
                    linecount -= 1

            else:
                print("Selecting score order in file %s not successful!" % ligidfile)
                sys.exit(1)
        else:
            print( "No such file  %s  in current folder! \n" % ligidfile)
            sys.exit(1)
        # ligand locations are in absolute paths
        return ligand

    def sepf2Cplx(self, receptorf, ligandf, outf, obabelexe):

        gentop = GenerateTop()

        if ".pdb" not in receptorf :
            gentop.runObabel(obabelexe=obabelexe, input=receptorf, output=receptorf.split(".")[0]+".pdb", verbose=True)

        if ".pdb" not in ligandf :
            gentop.runObabel(obabelexe=obabelexe, input=ligandf, output=ligandf.split(".")[0]+".pdb", verbose=True)

        with open(outf, "wb") as tofile :
            with open(receptorf.split(".")[0]+".pdb",) as linesr :
                for s in linesr :
                    if "ATOM" in s or "HETATM" in s and "**" not in s :
                        tofile.write(s)
            tofile.write("TER \n")
            with open(ligandf.split(".")[0]+".pdb",) as linesl :
                for s in linesl :
                    if "ATOM" in s or "HETATM" in s and "**" not in s :
                        tofile.write(s)

        return 1

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

