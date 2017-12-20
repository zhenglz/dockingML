#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import urllib
import time
import argparse
import subprocess as sp
from argparse import RawTextHelpFormatter
import math
from glob import glob
from collections import *
import linecache
import numpy as np
import glob
import time

from dockingML import extract

# import modeller for loop refinement
try:
    from modeller import *
    from modeller.automodel import *
    MODELLER_EXIST = True
except ImportError :
    print("Warning: Modeller is not well imported, some features may not work. ")
    MODELLER_EXIST = False

class RewritePDB :
    def __init__(self, filename):
        self.pdb = filename

    def pdbRewrite(self, output, chain, atomStartNdx, resStartNdx):
        # change atom id, residue id and chain id
        resseq = resStartNdx
        atomseq = int(atomStartNdx)
        chainname = chain
        lines = open(self.pdb)
        newfile = open(output,'w')
        resseq_list = []

        for s in lines :
            if "ATOM" in s and len(s.split()) > 6 :
                atomseq += 1
                newline = s
                newline = self.atomSeqChanger(newline, atomseq)
                newline = self.chainIDChanger(newline, chainname)
                if len(resseq_list) == 0 :
                    newline = self.resSeqChanger(newline, resseq)
                    resseq_list.append(int(s[22:26].strip()))
                else :
                    if resseq_list[-1] == int(s[22:26].strip()) :
                        newline = self.resSeqChanger(newline, resseq)
                    else :
                        resseq += 1
                        newline = self.resSeqChanger(newline, resseq)
                    resseq_list.append(int(s[22:26].strip()))
                newfile.write(newline)
            else :
                newfile.write(s)
        print "Completed!"

    # Change the ID of residue, or index of residue
    def resSeqChanger(self, inline, resseq):
        resseqstring = " "*(4- len(str(resseq)))+str(resseq)
        newline = inline[:22] + resseqstring + inline[26:]
        return newline

    def atomSeqChanger(self, inline, atomseq) :
        atomseqstring = " " * (5 - len(str(atomseq))) + str(atomseq)
        newline = inline[:6] + atomseqstring + inline[11:]
        return newline

    def resNameChanger(self, inline, resname) :
        resnamestr = " " * ( 4 - len(str(resname) ) ) + str(resname)
        newline = inline[:17] + resnamestr + inline[20:]
        return newline

    def chainIDChanger(self, inline, chainid) :
        newline = inline[:21] + str(chainid) + inline[22:]
        return newline

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



class AutoRunMD :
    def __init__(self, topFile, taskName, grompp, mdrun, verbose=True, qsub=False):
        self.top = topFile
        self.task= taskName
        self.verbose = verbose
        self.grompp = grompp
        self.mdrun = mdrun

    def modifyMDP(self, inputMDPFile, outputMDPFile, parameters):
        tofile = open(outputMDPFile, 'wb')
        with open(inputMDPFile) as lines :
            for s in lines :
                if len(s.split()) > 0 and s[0] != ";" :
                    if s.split()[0] in parameters.keys() \
                            and len(parameters[s.split()[0]]) > 0 :
                        tofile.write("%s    = ")
                        for item in parameters[s.split()[0]] :
                            tofile.write("%s "%str(item))
                        tofile.write(" \n")
                    else:
                        tofile.write(s)
                else :
                    tofile.write(s)
        return(outputMDPFile)

    def energyMinimization(self, emMDPFile, groFile, outTprFile, np=4, qsub=False):
        if self.verbose :
            print("Start energy minimization for task %s."%self.task)

        # generate GMX tpr file
        job = sp.check_output("%s -f %s -c %s -o %s -p %s "%
                              (self.grompp, emMDPFile, groFile, outTprFile, self.top),
                              shell=True )
        if self.verbose :
            print(job)
            print("Generating TPR file for EM completed.")

        # run MD with qsub or on the terminal
        cmd = "mpirun -np %d %s -deffnm %s " % (np, self.mdrun, outTprFile)
        if qsub :
            pscript = PrepScripts('sample_qsub.sh')
            pscript.prepareQsub('em_qsub.sh', cmd)
            job = sp.Popen("qsub em_qsub.sh", shell=True)
            job.communicate()

        else :
            job = sp.Popen(cmd, shell=True)
            job.communicate()

        return(outTprFile+".gro")

    def runMD(self, MDPFile, groFile, outTprFile, mdpOptions,
              sampleScript='sample_qsub.sh',
              qsub=False, np=4,
              index="index.ndx",
              ):

        # prepare mdp file for MD
        current_mdp = MDPFile
        if len(mdpOptions.keys()) :
            self.modifyMDP(inputMDPFile=MDPFile,
                           outputMDPFile="modified_"+MDPFile,
                           parameters=mdpOptions)
            current_mdp = "modified_"+MDPFile

        # prepare tpr file for md
        job = sp.check_output("%s -f %s -c %s -o %s -p %s -n %s"
                              %(self.grompp, current_mdp, groFile, outTprFile, self.top, index),
                              shell=True)
        if self.verbose :
            print(job)

        # run md
        cmd = "mpirun -np %d %s -deffnm %s " % (np, self.mdrun, outTprFile[-4:])
        if qsub :
            pscript = PrepScripts(sampleScript)
            oscript = str(len(glob("./*_"+sampleScript)))+"_"+sampleScript
            pscript.prepareQsub(oscript, cmd)
            job = sp.Popen("qsub %s"%oscript, shell=True)
            job.communicate()
        else :
            job = sp.Popen(cmd+" &", shell=True)
            job.communicate()

        return(outTprFile[-4:]+".gro")

def runGenTop() :
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
        os.environ["AMBERHOME"] = "/home/liangzhen/software/amber14/"
    except :
        os.system("export AMBERHOME=/home/liangzhen/software/amber14/")

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
        sh = top.runAntechamber(args.nc)

        try:
            # set env for AMBERHOME
            os.environ["AMBERHOME"] = "/home/liangzhen/software/amber14/"
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
        top.removeMolInfor(args.out)

if __name__ == "__main__" :
    # cd to PWD
    os.chdir(os.getcwd())

    if len(sys.argv) <= 2 :
        help_information = '''
        The tool is used for automating Gromacs based MD simulation.
        Several independent tools are provided.

        Separated tools including index, gentop, autoMD etc.

        Usage:

        python autoMD.py -h
        python autoMD.py gentop -h
        python autoMD.py index -h
        python autoMD.py extract -h
        '''
        print(help_information)
        sys.exit(1)
    else :
        if str(sys.argv[1]) == "index" :
            #print("Show help message: ")
            index = PdbIndex()
            index.genGMXIndex()
        #if 1 :
        elif str(sys.argv[1]) == "gentop" :
            runGenTop()

        elif str(sys.argv[1]) == "extract" :
            runExtract()
