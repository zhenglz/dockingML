#!/usr/bin/env pypy
# Post processing gold docking results
# Copyright @ ZHENG Liangzhen, astrozheng@gmail.com, Nov.8, 2016

import argparse
import glob
import os
import sys
import shutil
import subprocess as sp
from argparse import RawTextHelpFormatter
from collections import  OrderedDict
import linecache

def findOriginalLig(filename, ligid, type="mol2"):
    if type == filename[-4:] :
        lines = open(filename)

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

        return ligcontent

def findDockedLig(filepath, ligid, filenum=1) :

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

def rankingLst(outfile, shell=True, lstpath='output*/*.lst') :
    if shell :
        # using shell for ranking all the lst files exisited
        sp.Popen('awk \'$1 !~ /#/ && $1 !~ /-/ {print $0} %s\' | uniq | sort > %s '
                 %(lstpath, outfile),shell=True)
        return outfile
    else :
        # pythonic style ranking
        filedirect = os.listdir("./")
        for fd in filedirect :
            if lstpath.split("/")[0][:-1] in fd and os.path.isdir(fd) :
                files = os.listdir(fd)
                for f in files :
                    if ".lst" in f and not os.path.isdir(fd + "/" + f) :
                        pass
                        

def getLigandID(ligidinfor, pathname) :
    '''
    Determine whether the input is a file or a list
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

def getLigand(ligidfile, topNum, order):

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
            print "Selecting score order in file %s not successful!" % ligidfile
            sys.exit(1)
    else :
        print "No such file  %s  in current folder! \n" % ligidfile
        sys.exit(1)
    # ligand locations are in absolute paths
    return ligand

if __name__ == "__main__" :
    pwd = sp.check_output("pwd", shell=True)
    os.chdir(pwd[:-1])

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

    outputpath = args.opath
    if "/" != outputpath[-1]:
        outputpath += "/"

    dockedpath = args.dpath
    #if "/" != dockedpath[-1]:
    #    dockedpath += "/"

    ''' Prepare ligands for extracting '''

    if len(args.id) == 0 :
        print "Error! Ligand file or ID missing!"
        liglistDict = {}

        sys.exit(1)
    elif args.id not in os.listdir("./") :
        if args.rnk != "NA" :
            # not find the ligand id file and try to rank all the lst files
            print "Not find file %s, ranking the lst files firstly.\n" % args.id
            ligId = rankingLst("rankingall.lst",args.rnk)
            liglistDict = getLigand(ligId, args.topn, "a")
        else :
            print "Getting ranking ligands\' information failed!\nExit Now!\n"
            sys.exit(1)
    elif args.id in os.listdir("./"):
        liglistDict = getLigand(args.id, args.topn, args.ord)
    else :
        print "Getting ranking ligands\' information failed!\nExit Now!\n"
        sys.exit(1)

    ''' Extract ligand files from original ligand library file '''
    if True :

        solutions = []
        if len(liglistDict.keys()) > 0 :
            for ligid in liglistDict.keys() :
                print "\nNow working on ligand %s \n" % ligid

                if args.lib != "NULL" and args.lib != "NA" :
                    tofile = open(outputpath + ligid + "_orig.mol2", 'w')
                    ligcontent = findOriginalLig(args.lib, ligid, type="mol2")
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
            print "\n" + outputpath + "combined_soultions.mol2 file created! \nExit Now!\n"
    else :
        print "Error! Ligand Library file should be a mol2 file!"