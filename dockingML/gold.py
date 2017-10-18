#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Post processing of the GOLD docking results
"""

import argparse
import glob
import os
import sys
import linecache
import subprocess as sp

from .mol2IO import Mol2IO

from argparse import RawTextHelpFormatter
from collections import  OrderedDict


class GoldResults :
    def __init__(self):
        pass

    def findOriginalLig(self, input, ligid, type="mol2"):
        """
        get ligand data using ligand id from a multi-frame mol2 file
        :param input: str, a multi-frame mol2 file
        :param ligid: str, the id of the ligand
        :param type: the extension of the input
        :return: a list of lines
        """
        if type == input[-4:] and os.path.exists(input) :
            with open(input) as lines :
                ligcontent = []
                condition = 0
                for s in lines :
                    '''in a ligand content, condition = True '''
                    if "<TRIPOS>MOELECULES" in s :
                        condition += 1
                    if condition == 1 and len(s.split()) > 0 and ligid == s.split()[0] :
                        condition += 1
                        ligcontent.append("<TRIPOS>MOELECULES \n")
                        ligcontent.append(ligid+"  \n")
                    if "SUB" in s and "TEMP" in s and "ROOT" in s :
                        if condition == 2 :
                            ligcontent.append(s)
                        condition = 0
                    if condition == 2 :
                        ligcontent.append(s)
                lines.close()
        else :
            print("The file type is %s, while filename is %s"%(type, input))
            ligcontent = []

        return ligcontent

    def findDockedLig(self, filepath, ligid, filenum=1) :
        """
        from the GOLD docking results, select only some of the important files.
        based on their ligand id
        :param filepath: str, path of the gold docking results
        :param ligid: str, the id of the ligand
        :param filenum:
        :return: files together with their absoult path
        """

        filelist = glob.glob(filepath+"/gold_soln*.mol2")
        thefiles = []

        ''' going to find the gold results file based on ligand id '''
        for filename in filelist :
            with open(filename) as lines :

                for s in lines :
                    if len(s.split()) > 0 and s.split()[0] in ligid.strip('./') :
                       thefiles.append(filename)

            if len(thefiles) > filenum :
                break

        return thefiles

    def rankingLst(self, outfile, shell=True, lstpath='output*/*.lst') :
        """
        sort the bestranking.lst files
        :param outfile: str, output file name
        :param shell: bool, using shell to rank
        :param lstpath: str, the input of the lst files
        :return:
        """
        if shell :
            # using shell for ranking all the lst files exisited
            sp.Popen('awk \'$1 !~ /#/ && $1 !~ /-/ {print $0} %s\' | uniq | sort > %s '
                     %(lstpath, outfile),shell=True)
            return 1
        else :
            # pythonic style ranking
            filedirect = os.listdir("./")
            for fd in filedirect :
                if lstpath.split("/")[0][:-1] in fd and os.path.isdir(fd) :
                    files = os.listdir(fd)
                    for f in files :
                        if ".lst" in f and not os.path.isdir(fd + "/" + f) :
                            pass
            return 1

    def getLigandID(self, ligidinfor, pathname) :
        '''
        Determine whether the input is a file or a list
        :param ligidinfor:
        :param pathname:
        :return:
        '''

        if len(ligidinfor ) > 1 :
            '''multiple ligand id directly'''
            ligandlist = ligidinfor
        else :
            ''' supposing the result ligands ids in file, each record per line '''
            ligandlist = []
            if ligidinfor in glob.glob(pathname) :
                lines = open(ligidinfor[0])
                #ligandlist = []
                for s in lines :
                    if len(s.split()) > 0 and "#" not in s :
                        ligandlist.append(s.split()[0])
                lines.close()
            else :
                ligandlist.append(ligidinfor[0])

        return ligandlist

    def getLigand(self, ligidfile, topNum, order):
        '''
        get top N ligands from the lst file
        :param ligidfile:
        :param topNum:
        :param order:
        :return:
        '''

        # ligidfile is in relative path
        if ligidfile in os.listdir("./") :
            # this file exist and the format is similar to bestranking lst file
            # docking score in descending order
            linecount = 0
            ligand = OrderedDict()

            if order == "d" or order[0] == "d" :
                linecount = 1
                while linecount <= topNum :
                    s = linecache.getline(ligidfile, linecount)[:-1]  # remove the "\n"
                    ligand[s.split()[-1].strip("\'")] = s.split()[-2].strip("\'")
                    # increase line number to goto next line
                    linecount += 1
            elif order == "a" or order[0] == "a" :
                lines = open(ligidfile)
                nl = len(lines.readlines())
                linecount = nl
                if nl < topNum :
                    bot = 0
                else :
                    bot = nl - topNum
                while linecount > bot :
                    s = linecache.getline(ligidfile, linecount)[:-1]  # remove the "\n"
                    ligand[s.split()[-1].strip("\'")] = s.split()[-2].strip("\'")
                    linecount -= 1

            else :
                print( "Selecting score order in file %s not successful!" % ligidfile)
                sys.exit(1)
        else :
            print( "No such file  %s  in current folder! \n" % ligidfile)
            sys.exit(1)

        # ligand locations are in absolute paths
        return ligand

    def catGOLDResult(self, input_receptor, input_ligand ):
        """
        cat two mol2 files, or both pdb files are accepted
        :param input_receptor: str, file name
        :param input_ligand: str, filename
        :return:
        """

        mio = Mol2IO()

        if all([os.path.exists(input_receptor), os.path.exists(input_ligand)]) :
            mio.catenateMol(input_receptor, input_ligand, False)
            return 1
        else :
            print("input file (or files) missing. exit now")
            sys.exit(0)

    def removeLonePair(self, input, output, shell=True):
        '''
        input a pdb file, remove the lone pair electron lines in the file
        :param input: str, input file name
        :param output: str, output filename
        :return:
        '''
        if shell :
            if os.path.exists(input) :
                job = sp.Popen("awk \'$1 ~ /HETATM/ && $3 !~ /XX/ {print $0}\' %s > %s " % (input, output), shell=shell)
                job.communicate()
            else :
                print("input file %s not exist. exit now"%input)
                sys.exit(0)
        else :
            tofile = open(output, "wb")
            with open(input) as f :
                for s in f :
                    if "XX" not in s and "**" not in s :
                        tofile.write(s)
                    else :
                        pass
            tofile.close()

        return 1


if __name__ == "__main__" :

    description = ''' The script is used to extract specific structures after docking. '''

    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-id', type=str, default="ligands_id.dat",
                        help="Input file of ligands id for extracting.\n"
                             "Default is ligands_id.dat\n")
    parser.add_argument('-ord', type=str, default="d", help="The sort order in the ligid file.\n"
                                                            "Options: a(ascending) or d(decending).\n"
                                                            "Default is d (descending)\n")
    parser.add_argument('-lib', type=str, default="ligandlibrary.mol2",
                        help="The original ligand library file for screening. \n"
                             "If NA or NULL provided, do nothing. \n"
                             "Default is ligandlibrary.mol2 \n")
    parser.add_argument('-rnk', type=str, default='output*/*.lst',help="Ranking the lst files in output folders.\n"
                                                                       "You could enter output*/*.lst \n"
                                                                       "Default set is NA.")
    parser.add_argument('-dpath', type=str, default='./output/', help="The path of the docking results. ")
    parser.add_argument('-opath', type=str, default='./moloutput/', help="The path of the selected results. ")
    parser.add_argument('-topn',type=int, default=100, help="How many top ranking ligands for extracting. Default is 10. ")

    args = parser.parse_args()
    gr   = GoldResults()

    if len(sys.argv) < 3 :
        parser.print_help()
        parser.print_usage()

    outputpath = args.opath
    if "/" != outputpath[-1]:
        outputpath += "/"
    dockedpath = args.dpath

    ''' Prepare ligands for extracting '''

    if len(args.id) == 0 :
        print( "Error! Ligand file or ID missing!")
        liglistDict = {}

        sys.exit(1)
    elif args.id not in os.listdir("./") :
        if args.rnk != "NA" :
            # not find the ligand id file and try to rank all the lst files
            print( "Not find file %s, ranking the lst files firstly.\n" % args.id)
            ligId = gr.rankingLst("rankingall.lst",args.rnk)
            liglistDict = gr.getLigand(ligId, args.topn, "a")
        else :
            print( "Getting ranking ligands\' information failed!\nExit Now!\n")
            sys.exit(1)
    elif args.id in os.listdir("./"):
        liglistDict = gr.getLigand(args.id, args.topn, args.ord)
    else :
        print( "Getting ranking ligands\' information failed!\nExit Now!\n")
        sys.exit(1)

    ''' Extract ligand files from original ligand library file '''
    if True :

        solutions = []
        if len(liglistDict.keys()) > 0 :
            for ligid in liglistDict.keys() :
                print( "\nNow working on ligand %s \n" % ligid)

                if args.lib != "NULL" and args.lib != "NA" :
                    tofile = open(outputpath + ligid + "_orig.mol2", 'w')
                    ligcontent = gr.findOriginalLig(args.lib, ligid, type="mol2")
                    for s in ligcontent :
                        tofile.write(s)

                    tofile.close()
                else :
                    ## not copy original docking ligand file
                    pass

                ## Get Docked ligands and give the name of the docked file to solnfile and cp it to a new file.
                solnfile = liglistDict[ligid]

                sp.Popen("cp " + solnfile + " "+ outputpath+ligid+"_docked.mol2", shell=True)
                solutions.append(outputpath+ligid+"_docked.mol2")

            if "combined_soultions.mol2" in os.listdir(outputpath) :
                os.system("rm "+outputpath+"combined_soultions.mol2 -f")
            tofile = open(outputpath+"combined_soultions.mol2",'w')
            for file in solutions :
                with open(file) as lines :
                    for s in lines :
                        tofile.write(s)
            tofile.close()
            print("\n" + outputpath + "combined_soultions.mol2 file created! \nExit Now!\n")
    else :
        print( "Error! Ligand Library file should be a mol2 file!")