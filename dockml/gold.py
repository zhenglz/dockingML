#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Post processing of the GOLD docking results
"""

import argparse
import glob
import os
import sys
import numpy as np
import linecache
import subprocess as sp
from .mol2IO import Mol2IO
from argparse import RawTextHelpFormatter
from collections import OrderedDict
from collections import defaultdict

class GoldResults :
    def __init__(self, bestranking):
        self.bestrank = bestranking
        self.results = self.lstResults(self.bestrank)

    def lstResults(self, lst):
        '''
        read the bestranking.lst file
        duplicated results will be overided by last result record
        results: dict, key: ligand values: scores and file paths
        :return: defaultdict(list), the results
        '''

        results = defaultdict(list)

        with open(lst) as lines :
            for s in [ x for x in lines if ("#" not in x and len(x.split())) ] :
                results[ s.split()[-1] ] = s.split()[:-1]

        # { ligid: lstline }
        return results

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

    def findDockedLig(self, ligid) :
        """
        from the GOLD docking results, select the result pose mol2 file
        based on their ligand id
        :param ligid: str, the id of the ligand
        :return: filename
        """

        return self.results[ligid][-1]

    def getGoldScore(self, ligid):
        '''

        :param ligid: str, ligand name
        :return: float
        '''

        return float(self.results[ligid][0])

    def getElementWeight(self, inputFile="/home/liangzhen/bin/AtomType.dat"):
        mwElement = defaultdict(list)

        with open(inputFile) as lines :
            for s in lines :
                #mwparam = 0.0
                if "#" not in s and ";" not in s and len(s.split()) > 3 :
                    mwElement[s.split()[0]] = float(s.split()[2])
        return mwElement

    def molecularWeight(self, ligDockedFile):
        ## calculate molecular weight of a mol2 file
        mwElement = self.getElementWeight()
        molecularW = 0.0
        heavyAtomCount = 0

        atomsField = False
        #print "LIG DOCKED FILE HERE"
        if os.path.exists(ligDockedFile) :

            with open(ligDockedFile) as lines :
                for s in lines :
                    if len(s.split()) > 0 and "#" not in s :
                        #print "read lines here"
                        #print s

                        if "@<TRIPOS>ATOM" in s :
                            atomsField = True
                        elif "@<TRIPOS>BOND" in s :
                            atomsField = False

                        elif atomsField and "****" not in s :
                            ## heavy atom number count
                            heavyAtomCount += 1

                            if "H" not in s.split()[5] :
                                ## heavy atom number count
                                heavyAtomCount += 1
                            else :
                                pass

                            if s.split()[5].split(".")[0] not in mwElement.keys():
                                molecularW += mwElement["DU"]
                            else:
                                molecularW += mwElement[s.split()[5].split(".")[0]]

            return heavyAtomCount, molecularW
        else :
            print(ligDockedFile, "  Not exist")
            return -10000.0, 10000.0

    def rankingLst(self, outfile, shell=True, lstpath='output*/*.lst', decending=False) :
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
            files = glob.glob(lstpath)

            contents = []
            for fn in files :
                if os.path.exists(fn) :
                    with open(fn) as lines :
                        contents += [x.split() for x in lines.readlines() if "#" not in x ]

            contents.sort(key=lambda x: float(x[0]), reverse=decending )

            lines = [ " ".join(x) for x in contents ]

            with open(outfile) as tofile :
                map(lambda x: tofile.write(x), lines)

            return 1

    def sortResult(self, reverse=False):

        contents = []
        if os.path.exists(self.bestrank):
            with open(self.bestrank) as lines:
                for s in [ x for x in lines if ("#" not in x and len(x.split()))] :
                    contents.append(s.split())
        #ligfiles = [ x[-2].strip("\'") for x in contents ]

        #mws = map(self.molecularWeight, ligfiles)
        #print(list(mws))

        #scores = [ float(x[0]) for x in contents ]
        #print(scores)

        #new_scores = map(lambda x, y: x/y[-1], scores, mws)

        return sorted(contents, key=lambda x: float(x[0]), reverse=reverse)

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

                ligandlist.append( Mol2IO().properties(ligidinfor, "MOLECULE")[0] )

                '''
                lines = open(ligidinfor[0])
                #ligandlist = []
                for s in lines :
                    if len(s.split()) > 0 and "#" not in s :
                        ligandlist.append(s.split()[0])
                lines.close()
            else :
                ligandlist.append(ligidinfor[0])'''

        return ligandlist

    def getLigand(self, lst, topNum, reverse=False):
        '''
        get top N ligands' file pathnames from the lst file
        as well as their names
        :param ligidfile:
        :param topNum:
        :param reverse:
        :return:
        '''

        self.rankingLst("sorted_list", lst, decending=reverse)

        sorted = self.lstResults("sorted_list")

        ligands = sorted.keys()[-topNum:]
        filenames = [ x[-1] for x in sorted.values()[-topNum:]]

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
        '''

        # ligand locations are in absolute paths
        return filenames, ligands

    def copyLigandPoseFile(self, pose, out):
        try :

            job = sp.Popen("cp {} {} ".format(pose, out), shell=True)
            job.communicate()
        except FileNotFoundError :
            print(pose, " Not Found!")

        return 1

    def catGOLDResult(self, input_receptor, input_ligand):
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

class RankResults :

    def __init__(self, lst):
        self.bestrank = lst

    def sortLst(self, reverse=False):
        '''
        sort a bestranking.lst file
        :param lst: str, file name
        :param reverse:
        :return: a 2d array, [ [score, ..., ligid], ]
        '''
        return GoldResults(self.bestrank).sortResult(reverse=reverse)

    def getTopLigandsID(self, topNum=100):
        '''
        return the ligand ids in the top ranking list
        :param topNum:
        :return:
        '''

        return [ x[-1] for x in self.sortLst(reverse=True)[: topNum]]

    def commonTopLigandsID(self, lst_list, topNum=1000):
        '''
        input a list of bestranking files, return their common ranking ligands
        :param lst_list:
        :param topNum:
        :return:
        '''

        com_top = set([])

        for i in range(len(lst_list)) :
            topligs = RankResults(lst_list[i]).getTopLigandsID(topNum)

            if i == 0 :
                com_top = set(topligs)
            else :
                com_top = com_top & set(topligs)

        return com_top

'''
if __name__ == "__main__" :

    description = # The script is used to extract specific structures after docking. #

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
    gr   = GoldResults(args.rnk)

    if len(sys.argv) < 3 :
        parser.print_help()
        parser.print_usage()

    outputpath = args.opath
    if "/" != outputpath[-1]:
        outputpath += "/"
    dockedpath = args.dpath

    # Prepare ligands for extracting #

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

    # Extract ligand files from original ligand library file #
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

'''