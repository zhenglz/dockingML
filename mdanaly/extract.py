#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Extract PDB files
"""

import sys
import argparse
from argparse import RawTextHelpFormatter


class ExtractPDB(object):
    def __init__(self, filename=""):
        self.fn = filename

    def extract_pdb(self, filename, structname, first_n_frame):
        """From a big pdb file to extract single PDB structure file

        Parameters
        ----------
        filename: str,
            the input filename
        structname: str,
            the output separated filename header
        first_n_frame: int,
            only output the first n frames

        Returns
        -------

        """

        lines = open(filename)
        file_no = 1
        pdbline = open(structname+'_'+str(file_no)+'.pdb', 'w')

        for s in lines :

            if "MODEL" in s:
                if file_no != 1:
                    pdbline = open(structname+'_'+str(file_no)+'.pdb', 'w')
                pdbline.write(s)
                print("Start Model " + str(file_no))
                file_no += 1

            elif "AT" in s.split()[0] or "TER" == s.split()[0]:
                pdbline.write(s)

            elif "ENDMDL" in s:
                pdbline.write(s)
                pdbline.close()

            else:
                pass

            if first_n_frame != 0 and file_no > first_n_frame + 1 :
                print("Extract only the first "+str(first_n_frame))
                break

        print("Finished!\n\n")

    def extract_frame(self, filename, structname, no_frames=[]):
        """Extract only selected frames.

        Parameters
        ----------
        filename: str,
            the input pdb filename
        structname: str,
            the output separated filename header
        no_frames: list,
            the list of frames ids to be extracted

        Returns
        -------

        """

        lines = open(filename)
        print("The models of the pdb file is : ")
        for s in lines:
            if "MODEL" in s:
                print("    "+s[:-1])
        lines.close()

        if not len(no_frames):
            try:
                print("Which frames would you want to extract ? ")
                frames = input("Input the frame number(s) here (multi numbers are accepted):  ")
                frame_list = [int(x) for x in frames.split()]
            except IOError:
                print("You haven't select correct frames.")
                frame_list = []
        else:
            frame_list = no_frames

        for frame_no in frame_list:

            lines = open(filename)
            condition = False
            for s in lines:
                if "MODEL" in s and int(frame_no) == int(s.split()[1]):
                    newline = open(structname+"_"+str(frame_no)+".pdb", "w")
                    newline.write(s)
                    condition = True
                elif ("ATOM" in s or "HETATM" in s) and condition:
                    newline.write(s)
                elif condition and "ENDMDL" in s:
                    condition = False
                elif "MODEL" in s and int(frame_no)+1 == int(s.split()[1]):
                    condition = False
                    break
                else:
                    condition = False
            newline.close()
            lines.close()

        print("Finished writing frames to separated files!\n\n")
        return 1

    def printinfor(self):
        """
        Print instructions

        Returns
        -------

        """

        d = """
        What would you like to do now ?
        1. Extract all the frames from the input pdb file;
        2. Extract selected frames from the input file;
        3. Extract the first N frames from the input file;
        4. Do nothing and exit now.
        """

        print(d)

    def indexCoord(self, filename):
        """Provide a large coordination file, eg, a pdb file or a mol2 file,
        return the indexing of the frames

        Parameters
        ----------
        filename: str,
            the file name of a large multiple-frames coordination file
            either a pdb, a pdbqt or a mol2 file

        Returns
        -------
        indexing: list,
            the indexing of the first line of a multiple-frame file

        """

        indexing = []
        lineNumber = -1

        extention = filename.split(".")[-1]

        with open(filename) as lines:
            if extention in ['pdb', 'pdbqt']:
                for s in lines :
                    lineNumber += 1
                    if len(s.split()) > 1 and "MODEL" == s.split()[0]:
                        indexing.append(lineNumber)
            elif extention in ['mol2']:
                for s in lines :
                    lineNumber += 1
                    if "@<TRIPOS>MOLECULE" in s:
                        indexing.append(lineNumber)
            else:
                print("Only a pdb file, a pdbqt file or a mol2 file supported.")
                sys.exit(0)

        return indexing

    def extract_all(self, filename, structname):

        filenum = 1
        tofile = open(structname + "_" + str(filenum) + '.pdb', 'w')

        with open(filename) as lines :
            for s in lines :
                if "MODEL" in s :
                    tofile = open(structname + "_" + str(filenum) + '.pdb', 'w')
                    print("START FRAME: %d " % (filenum))

                elif "ENDMDL" in s:
                    tofile.write(s)
                    tofile.close()
                    filenum += 1
                    print("START FRAME: %d " % (filenum))
                elif "ATOM" in s :
                    tofile.write(s)

        return None

    def arguments(self):
        """
        argument prarser
        :return:
        """
        # PWD, change directory to PWD
        #os.chdir(os.getcwd())

        d = '''
        This script try to extract PDB frames from a long trajectory file from GMX trjconv or g_cluster.
        Any questions, contact Liangzhen Zheng, astrozheng@gmail.com

        Examples :
        Get help information
        python autoMD.py extract -h
        Extract frames in a multiple-frames PDB file, which contains "MODEL 1", "MODEL 2"

        autoMD.py extract -i md.pdb -o splitted

        after which interactive mode is entered and you are promoted to choose one of the 4 options:
            1. Extract all the frames from the input pdb file;
            2. Extract selected frames from the input file;
            3. Extract the first N frames from the input file;
            4. Do nothing and exit now.
        '''
        parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
        parser.add_argument('-i', '--input', type=str, default="yourpdb.pdb",
                            help="Input PDB file name. Default is yourpdb.pdb.")
        parser.add_argument('-o', '--output', type=str, default="out",
                            help="Output file format. Default is out_*.pdb.")
        parser.add_argument('-m', '--allMol2', type=bool, default=False,
                            help='Extract all the frames in a mol2 file. \n'
                                 'Options: True, False \n'
                                 'Default is False.')
        options, args = parser.parse_known_args()

        if len(sys.argv) < 2:
            parser.print_help()
            sys.exit(1)
        else:
            parser.print_help()

        return(options)

    def runExtract(self):
        """run ExtractPDB class
        extract frames from a multiple-frames pdb file or mol2 file

        Returns
        -------

        """

        args = self.arguments()

        pdbfile = args.input
        structname = args.output

        if args.allMol2:
            # if we choose to extract all frames in a mol2 file
            self.extract_all(args.input, args.output)
        else:
            self.printinfor()
            command = input("Your choice:  ")
            while command in ['0', '1', '2', '3', '4']:
                if command == "1":
                    self.extract_all(pdbfile, structname)
                    command = '0'
                elif command == "2":
                    self.extract_frame(pdbfile, structname)
                    command = '0'
                elif command == "3":
                    fnframe = input("  The first N frames output : N = ")
                    self.extract_pdb(pdbfile, structname, int(fnframe))
                    command = '0'
                elif command == "4" or command not in "01234":
                    print("\nExit now. \nHave a nice day!")
                    command = '5'
                    sys.exit(1)
                    # return None
                elif command == "0":
                    self.printinfor()
                    command = input("Your choice:  ")


def main():

    ext = ExtractPDB("")
    ext.runExtract()

