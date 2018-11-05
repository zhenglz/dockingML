#!/usr/bin/env python

from automd import gentop, fixpdb
from mdanaly import gmxcli
import sys, os

if __name__ == "__main__":

    d = """
    Generate gromacs style toplogy files from pdb files. 
    
    """

    parser = gmxcli.GromacsCommanLine(d=d)

    parser.arguments()

    parser.parser.add_arguments("-p", type="str", default="topol.top",
                                help="Output, optional. Output format: .top \n"
                                     "The output file name of the generated gromacs \n"
                                     "style topology files. ")
    parser.parser.add_arguments("-ff", type="str", default=["gaff", ], nargs="+",
                                help="Output, optional. Output format: .top \n"
                                     "The output file name of the generated gromacs \n"
                                     "style topology files. ")
    parser.parser.add_arguments("-amberhome", type="str", default="/usr/local/amber/",
                                help="Input, optional. \n"
                                     "The AMBERHOME environment for topology generation. \n")

    args = parser.parse_arguments()

    cmd = "export PATH=%s:$s:$PATH" % (args.amberhome, args.amberhome + "/bin")
    os.system(cmd)

    fpdb = fixpdb.FixPDB()
    fpdb.addhydrogenReduce(args.f, "H_"+args.f, args.amberhome+"/bin/reduce", flipH=False, verbose=True)

    pdbin = "H_"+args.f

    #Chem.rdmol
    gtop = gentop.GenerateTop()
    spdb = fixpdb.SummaryPDB("H_"+args.f)
    net_charge = spdb.netCharges(pdbin, "UNK")
    ftype = "AMBER"

    # ligand prep and frcmod files
    gtop.runAntechamber(pdbin, net_charge, ftype)

    frcmod = ftype + ".frcmod"
    prep   = ftype + ".prep"

    gtop.gmxTopBuilder(
        PDBFile=args.f, outputName=args.p, frcmodFile=frcmod,
        prepFile=prep, FField=args.ff, ionName=[], ionNum=[],
        solveBox=None, amberhome=args.amberhome, verbose=True,
    )




