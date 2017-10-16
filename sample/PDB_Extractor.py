#!/usr/bin/env python

import os
import sys
import argparse
import subprocess as sp
from argparse import RawTextHelpFormatter

class ExtractPDB :
    def __init__(self):
        pass
    
    def extract_all(self,filename,structname, first_n_frame):
        """ from a big pdb file to extract single PDB structure file """
    
        lines = open(filename)
        file_no = 1
        pdbline = open(structname+'_'+str(file_no)+'.pdb','w')
    
        for s in lines :
        
            if  "MODEL" in s :
                if file_no != 1 :
                    pdbline = open(structname+'_'+str(file_no)+'.pdb','w')
                pdbline.write(s)
                print "Start Model " + str(file_no)
                file_no += 1
                
            elif "ATOM" == s.split()[0] or "TER" == s.split()[0] :
                pdbline.write(s)
        
            elif "ENDMDL" in s :
                pdbline.write(s)
                pdbline.close()
            
            else :
                pass
            
            if first_n_frame != 0 and file_no > first_n_frame + 1 :
                print "Extract only the first "+str(first_n_frame)
                break
    
        print "Finished!\n\n"
    
    def extract_frame(self,filename,structname):
        lines = open(filename)
        print "The models of the pdb file is : "
        for s in lines :
            if "MODEL" in s :
                print "    "+s[:-1]
        lines.close()
        
        print "Which frames would you want to extract ? "
        frames = raw_input("Input the frame number(s) here (multi numbers are accepted):  ")
        
        frame_list = frames.split()
        for frame_no in frame_list :
            
            lines  = open(filename)
            condition = False
            for s in lines :
                if "MODEL" in s and int(frame_no) == int(s.split()[1])  :
                    newline = open(structname+"_"+str(frame_no)+".pdb","w")
                    newline.write(s)
                    condition = True
                elif "ATOM" in s and condition :
                    newline.write(s)
                elif condition and "ENDMDL" in s :
                    condition = False
                elif "MODEL" in s and int(frame_no)+1 == int(s.split()[1]) :
                    condition = False
                    break 
                else :
                    condition = False
            newline.close()
            lines.close()
        
        print "Finished writing frames to separated files!\n\n"

    def printinfor(self):
        print "What would you like to do now ?"
        print "  1. Extract all the frames from the input pdb file;"
        print "  2. Extract selected frames from the input file;"
        print "  3. Extract the first N frames from the input file;"
        print "  4. Do nothing and exit now. "

if __name__ ==  "__main__" :
    # PWD, change directory to PWD
    pwd = sp.check_output("pwd", shell=True)
    os.chdir(pwd[:-1])

    d = '''
    This script try to extract PDB frames from a long trajectory file from GMX trjconv or g_cluster.
    Any questions, contact Liangzhen Zheng, astrozheng@gmail.com

    Examples :
    Get help information
    PDB_Extractor.py -h
    Extract frames in md.pdb file, which contains "MODEL 1", "MODEL 2"

    PDB_Extractor.py -i md.pdb -o splitted
    after which interactive mode is entered and you are promoted to choose one of the 4 options.
        1. Extract all the frames from the input pdb file;
        2. Extract selected frames from the input file;
        3. Extract the first N frames from the input file;
        4. Do nothing and exit now.
    '''
    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-i', '--input',type=str, default="yourpdb.pdb",help="Input PDB file name. Default is yourpdb.pdb.")
    parser.add_argument('-o', '--output',type=str, default="out", help="Output file format. Default is out_*.pdb.")
    args = parser.parse_args()

    if len(sys.argv) < 2 :
        parser.print_help()
        sys.exit(1)
    else :
        parser.print_help()

    pdbfile = args.input
    if ".pdb" not in pdbfile :
        pdbfile += ".pdb"
    if pdbfile not in os.listdir("./") :
        print "\nNo such file %s in current directory! \nTry again!\nExit now!" % pdbfile
        sys.exit(0)

    structname = args.output

    pdb = ExtractPDB()
    pdb.printinfor()
    command = raw_input("Your choice:  ")
    while command in ['0','1','2','3','4'] :
        if command == "1" :
            pdb.extract_all(pdbfile, structname, 0)
            command = '0'
        elif command == "2" :
            pdb.extract_frame(pdbfile,structname)
            command = '0'
        elif command == "3" :
            fnframe=raw_input("  The first N frames output : N = ")
            pdb.extract_all(pdbfile, structname, int(fnframe))
            command = '0'
        elif command == "4" or command not in "01234" :
            print "\nExit now. \nHave a nice day!"
            command = '5'
            sys.exit(0)
            #return None
        elif command == "0" :

            pdb.printinfor()
            command = raw_input("Your choice:  ")

