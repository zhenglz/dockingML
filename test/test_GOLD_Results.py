#!/usr/bin/env python

from dockml import gold
import sys, os
from dockml import convert
from dockml import pdbIO
from automd import cleanpdb

if __name__ == "__main__" :

    gout = gold.GoldResults(sys.argv[1])

    results = gout.results

    convert = convert.Convert()
    pio     = pdbIO.rewritePDB("")

    all_ligs = list(results.keys())

    for i in range(len(all_ligs)) :
        lig = all_ligs[i]

        print("LIG {}, {} out of {}".format(lig, i, len(all_ligs)))
        pose = results[lig][-1].strip("\'").replace( "data", "scratch")

        gout.copyLigandPoseFile(pose, lig.strip("\'")+".mol2")

        convert.convert(pose, lig.strip("\'")+".pdb")

        try :
            cleanpdb.CleanPDB(pose).removeLonePair(lig.strip("\'")+".pdb", "noL"+lig.strip("\'")+".pdb")
        except :
            pass

        try :
            os.system("cat {} {} > {}".format(sys.argv[2], "noL"+lig.strip("\'")+".pdb",
                                              "complex_%s.pdb"%lig.strip('\'')))

            #pio.combinePDBFromLines("complex_%s.pdb"%lig.strip('\''), pro_lines+lig_lines)
        except FileNotFoundError :
            pass

        os.system("rm -f {} {} {}".format(lig.strip("\'")+".pdb",
                                          "noL"+lig.strip("\'")+".pdb",
                                          lig.strip("\'") + ".mol2"))



