#!/usr/bin/env python

from .autoMD import GenerateTop, MolDocking, CleanPDB, SummaryPDB
import os, sys
import subprocess as sp

"""
Combine the mol2 receptor and mol2 ligand files,
and generate a complex PDB file
And generate topology files for MD simulation
"""
def prepITP(oldmol2) :
    # collect the mol2 files, and prepare their pdb file
    # calculate the bcc charges of the molecules
    mold = MolDocking()
    gentop = GenerateTop()
    newpdb = oldmol2.split(".")[0]+".pdb"
    print("Generate a new pdb file")
    gentop.runObabel('obabel', oldmol2, newpdb)

    cpdb = CleanPDB(newpdb, "obabel")
    newpdb2 = "cleaned_" + oldmol2.split(".")[0]+".pdb"
    cpdb.removeLonePair(newpdb, newpdb2)
    newmol2 = newpdb2.split(".")[0]+".mol2"
    gentop.runObabel("obabel", newpdb2, newmol2)

    spdb = SummaryPDB(newpdb2, aminoLibFile="amino-acid.lib")
    nc = spdb.netCharges(oldmol2)

    # now generate new itp file
    mold.atomicChargesLig('obabel', newmol2, netcharges=nc)
    # this will generating a itp file, name "cleaned_" + oldmol2.split(".")[0] + ".itp"

if __name__ == "__main__" :
    # pwd
    os.chdir(os.getcwd())

    mol2d = MolDocking()
    count = 0
    input = open("input.dat", "wb")
    with open("bestranking.lst") as lines :
        for s in lines :
            if "#" not in s :
                count += 1
                filen = s.split()[-2].strip("\'").split("/")[-1]
                newn  = filen.split("_")[-2]

                job = sp.Popen("cp ../fda_3000_docking/fda_lig/%s ./%s.mol2 " % (filen, newn), shell=True)
                job.communicate()

                prepITP(newn)
                mol2d.sepf2Cplx("protein.pdb", "cleaned_%s.pdb"%newn, 'cplx_%s.pdb'%newn, 'obabel')

                input.write("%s   LIG   %s  \n" % ('cplx_%s.pdb'%newn, newn))

                if count % 10 == 0 :
                    print("Now Processing %d ligand!  \n" * 5)

                if count > 3 :
                    input.close()
                    sys.exit(0)
    input.close()
