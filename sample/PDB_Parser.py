import os
import sys
import time
import numpy as np
import pandas as pd
from collections import defaultdict
import urllib

class PDB_Parser :
    def __init__(self, inPDB='1a28.pdb'):
        self.pdbin = inPDB

    def pdbListInfor(self, pdbList):
        #pdbInfor = defaultdict(list)
        pdbCode, resolution, ligandName, Ktype, Kvalue = [], [], [], [], []
        with open(pdbList) as lines :
            for s in lines :
                if len(s.split()) > 3 and "#" not in s :
                    pdbCode.append(s.split()[0])
                    resolution.append(s.split()[1])
                    ligandName.append(s.split()[-2])
                    Ktype.append(s.split()[-1][:2])
                    Kvalue.append(s.split()[-1][3:])
        pdbInfor = pd.DataFrame({
            "PDBCode" : pdbCode,
            "Resolution" : np.asarray(resolution ),
            "LigandName" : ligandName,
            "AffinityType" : Ktype,
            "AffinityValue": Kvalue,
        })

        return pdbInfor

    def parsePDB(self, pdbin, proteinChains, pdbOut, ligandInfor=["SUB","A"]):
        pdbout = open(pdbOut, 'w')

        with open(pdbin) as lines :
            for s in lines :
                if "ATOM" in s :
                    if s[21] in proteinChains :
                        pdbout.write(s)
                elif "HETATM" in s :
                    if s[21] == ligandInfor[1] and s[17:20].strip() == ligandInfor[0] :
                        pdbout.write(s)
                else :
                    pass

        pdbout.close()
        return 1

if __name__ == "__main__" :
    ## PWD
    pwd = os.getcwd()
    os.chdir(pwd)

    pp = PDB_Parser()

    #infor = pp.pdbListInfor('sorted_infor_HIV1.txt')
    infor = pd.read_csv('./Clean_PDB_List.csv', sep=',', header=0)

    ML_input = open("ML_input.dat",'w')
    for index in range(len(infor["PDBID"])) :
        pdbid = infor["PDBID"][index]
        if not os.path.exists(pdbid + ".pdb" ) :
            tofile = open(pdbid + ".pdb", 'w')
            ## PDB file not exist, download from RCSB Protein Data Bank
            url = "http://www.rcsb.org/pdb/files/%s.pdb" % pdbid
            tofile.write(urllib.urlopen(url).read())
            tofile.close()

        proteinchains = ["A","B"]
        pdbout = pdbid + "_clean" + ".pdb"
        ligandName = infor["Chain_Lig"][index].split("|")[1].split("_")[0]
        ligandChain= infor["Chain_Lig"][index].split("|")[1].split("_")[1]

        pp.parsePDB(pdbid + ".pdb", proteinChains=proteinchains,pdbOut=pdbout, ligandInfor=[ligandName, ligandChain])

        ML_input.write("%s   %s  \n" % (pdbout, ligandName))

    ML_input.close()
