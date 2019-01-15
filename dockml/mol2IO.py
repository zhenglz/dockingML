#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from glob import glob
import subprocess as sp
from random import random

class Mol2IO :

    """
    Manipulating mol2 files

    Notes:
        property is the keyword after the @<TRIPOS> flag
    """

    def __init__(self):

        self.startProp = "MOLECULE"

    def properties(self, input, property="", lineindex=1):
        """
        extract the properties flag of a mol2 file
        :return:
        """
        props = []
        with open(input) as lines :
            lines = [x for x in lines if len(x.split())]
            condition = False
            count     = -1
            for x in lines :
                if "@<TRIPOS>" + property == x.split()[0] :
                    condition = True
                    count = 0

                if condition and "@<" not in x :
                    count += 1
                    if count == lineindex :
                        props.append(x.split()[0])

        return props



    def triposInfor(self, input, property):
        """
        get the property information of the input file
        :param input:
        :param property:
        :return:
        """

        with open(input) as lines :
            condition = False
            inf       = []
            lines = [x for x in lines if len(x.split())]
            for x in lines :
                if "@<TRIPOS>" in x and property in x :
                    condition = True

                if  "@<TRIPOS>" in x and property not in x :
                    condition = False

                if condition and "@<" not in x:
                    inf += x.split()
        return inf

    def splitMol2(self, input, output_base):
        """
        split the large mol2 file into single frame mol2 files
        :param input:
        :param output_base:
        :return:
        """
        count = 0
        if os.path.exists(input) :
            tofile = open(output_base + str(count) + ".mol2", "wb")
            with open(input) as lines :
                for s in [ x for x in lines if x[0] != "#"] :
                    if "@<TRIPOS>" + self.startProp == s.strip() :
                        if count :
                            tofile.close()

                        tofile = open(output_base+str(count)+".mol2", "wb")
                        tofile.write(s)

                        count += 1
                    else :
                        tofile.write(s)
            tofile.close()

        return count

    def selectFrame(self, inputs, molID):
        """
        select the frame with the same molid
        :param inputs: a list of files
        :param molID: the molecule name
        :return: the filename of the selected frame
        """

        for filen in inputs :
            if os.path.exists(filen) :
                inf = self.triposInfor(filen, self.startProp)
                if isinstance(inf[0], list) :
                    if inf[0][0] == molID :
                        return filen

                else :
                    if inf[0] == molID :
                        return filen

        return ""

    def catenateMol(self, input1, input2, output, shell=False):
        """
        catenate two mol2 files, or pdb files
        :param input1: str, filename
        :param input2: str, filename
        :param output: str, filename
        :return:
        """

        #tofile = open(output, "wb")


        if os.path.exists(input1) and os.path.exists(input2) :
            if shell :
                if output in [input2, input1] :
                    temp = "temp"+ str(random()*1000)

                    if os.path.exists(output) :
                        os.rename(output, temp)
                    else :
                        job = sp.Popen("touch %s"%temp, shell=True)
                        job.communicate()

                    job = sp.Popen("cat %s %s > %s " %(input1, input2, temp), shell=shell)
                    job.communicate()
                    os.rename(temp, output)

                else :
                    job = sp.Popen("cat %s %s > %s " %(input1, input2, output), shell=shell)
                    job.communicate()

            else :
                content = []
                with open(input1) as lines :
                    content += lines
                with open(input2) as lines :
                    content += lines

                tofile = open(output, "wb")
                for s in content :
                    tofile.write(s)
                tofile.close()

        else :
            print("Files %s and %s not exist !"%(input1, input2))

        return 1

    def removeDuplicated(self, input, output, property, remove_splited=True):

        """
        remove duplicated frames in the input muti-frame mol2 file
        :param input: str, mutiple-frame mol2 file
        :param output: str, mutiple-frame mol2 file
        :return:
        """

        # get a list of molecule id in the muti-frame mol2
        mol_ids = set(self.properties(input, property))
        print(mol_ids)

        # create seperated frame files
        if not os.path.exists("./splited/") :
            os.makedirs("./splited/")
        self.splitMol2(input, output_base="./splited/S_")

        # get necessary files
        uni_files = []
        mol2 = glob("./splited/S_*.*")
        for filen in mol2 :
            id = self.triposInfor(filen, property)[0]
            if id in mol_ids :
                uni_files.append(filen)
                mol_ids.remove(id)

        # cat files all together
        os.system("touch ./temp")
        print(uni_files)
        for filen in uni_files :
            if os.path.exists(filen) :
                self.catenateMol("./temp", filen, "./temp", shell=True)

        # rename file
        os.rename("./temp", output)
        if remove_splited :
            sp.Popen("rm -rf ./splited/* ", shell=True)

        return 1

    def changingProperty(self, input, output, property, index, newInfor):

        tofile = open(output, "wb")

        with open(input) as lines :
            condition = False
            count     = -1
            for s in lines :
                if "@<TRIPOS>" + property == s.split()[0] :
                    condition = True
                    count    += 1
                    tofile.write(s)

                elif "@<TRIPOS>" in s.split()[0] and property not in s.split()[0] :
                    condition = False
                    tofile.write(s)

                elif condition :
                    count += 1
                    if count == index :
                        tofile.write(newInfor+" \n")
                    else:
                        tofile.write(s)
                else :
                    tofile.write(s)
        tofile.close()

        return 1

