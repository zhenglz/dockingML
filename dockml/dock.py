#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import subprocess as sp
from glob import glob
from collections import *
import linecache
import glob
from dockml import convert


class VinaDocking(object):

    def __init__(self, vinaexe="vina"):
        self.vina = vinaexe

    def runVina(self, vinaexe, receptor, ligand,
                output='result.pdbqt', logfile='result.log',
                ncpu=1, exhaustiveness=32,
                center=[0, 0, 0], sizes=[40, 40, 40],
                no_modes=20, en_range=5, seed=-1,
                ):
        """
        Perform molecular docking using autodock vina.

        Parameters
        ----------
        vinaexe: str,
            executable vina binary file
        receptor: str,
            receptor file name in pdbqt format
        ligand: str,
            ligand file name in pdbqt format
        output: str,
            ligand binding pose, in pdbqt format
        logfile: str,
            docking results log
        ncpu: int,
            number cores for cpu calculation
        exhaustiveness: int,
            how accurate the vina docking
        center: list,
            the x y z center coordinates
        sizes: list,
            size of the binding pocket
        no_modes: int,
            number of binding poses in output
        en_range: int,
            the vina output energy range
        seed: int,
            random state seed number

        Returns
        -------
        self

        """

        try:
            job = sp.Popen('%s --receptor %s --ligand %s '
                           '--center_x %f --center_y %f --center_z %f '
                           '--size_x %f --size_y %f --size_z %f '
                           '--log %s --out %s --cpu %d '
                           '--exhaustiveness %d --num_modes %d '
                           '--energy_range %d --seed %d ' %
                           (self.vina, receptor, ligand,
                            center[0], center[1], center[2],
                            sizes[0], sizes[1], sizes[2],
                            logfile, output, ncpu,
                            exhaustiveness, no_modes,
                            en_range, seed),
                           shell=True
                           )
            job.communicate()
            job.terminate()
        except IOError:
            print("Docking molecule %s to %s using %s failed. \n"
                  "Check your input and logfile.")

            job = sp.check_output('%s --help ' % self.vina)
            print(job)

        return None

    def rankVinaResults(self, logfileList):
        """
        Obtain the binding energy score from vina log files

        Parameters
        ----------
        logfileList: list,
            a list of output log file names

        Returns
        -------
        a list of tuples, each key matches with a result list (top 3 results)
        """

        vinaResults = defaultdict(list)

        for resultfile in logfileList:
            condition = -1
            vinaResults[resultfile] = []
            with open(resultfile) as lines:
                for s in lines:

                    if "Refining results ... done" in s:
                        condition += 1
                    elif "affinity" in s:
                        condition += 1
                    else:
                        pass

                    if condition:
                        if len(s.split()) and s.split()[0] in ['1', '2', '3']:
                            vinaResults[resultfile].append(float(s.split()[1]))

        return sorted(vinaResults.items(), key=lambda x: x[1])


class GoldDocking(object):

    def __init__(self):
        pass

    def findOriginalLig(self, filename, ligid, type="mol2"):

        if type == filename.split(".")[-1] and os.path.exists(filename):
            with open(filename) as lines :

                ligcontent = []
                condition = 0
                for s in lines:
                    '''in a ligand content, condition = True '''
                    if "<TRIPOS>MOELECULES" in s:
                        condition += 1

                    elif condition == 1:
                        if ligid == s.split()[0]:
                            condition += 1
                            ligcontent.append("<TRIPOS>MOELECULES \n")
                            ligcontent.append(ligid + "  \n")
                        else :
                            condition = 0

                    elif condition == 2:
                        ligcontent.append(s)
                        condition = 0

                    if condition >= 3:
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
            print("No such file  %s  in current folder! \n" % ligidfile)
            sys.exit(1)
        # ligand locations are in absolute paths
        return ligand

    def sepf2Cplx(self, receptorf, ligandf, outf, obabelexe="obabel"):

        """
        Concatenate ligand-pose file and receptor file to a complex file.

        Parameters
        ----------
        receptorf: str,
            receptor file name
        ligandf: str,
            ligand pose file name
        outf: str,
            output complex file name
        obabelexe: str
            openbabel obabel command

        Returns
        -------
        self
        """

        converter = convert.Convert(obabel=obabelexe)

        if ".pdb" not in receptorf:
            converter.convert(input=receptorf, output=receptorf.split(".")[0] + ".pdb", verbose=True)

        if ".pdb" not in ligandf:
            converter.convert(input=ligandf, output=ligandf.split(".")[0] + ".pdb", verbose=True)

        with open(outf, "wb") as tofile:
            with open(receptorf.split(".")[0] + ".pdb") as linesr:
                for s in linesr:
                    if "ATOM" in s or "HETATM" in s and "**" not in s:
                        tofile.write(s)
            tofile.write("TER \n")
            with open(ligandf.split(".")[0] + ".pdb",) as linesl :
                for s in linesl :
                    if "ATOM" in s or "HETATM" in s and "**" not in s :
                        tofile.write(s)

        return self

