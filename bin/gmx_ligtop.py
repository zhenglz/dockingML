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

    parser.parser.add_argument("-p", type=str, default="topol.top",
                               help="Output, optional. Output format: .top \n"
                                    "The output file name of the generated gromacs \n"
                                    "style topology files. ")
    parser.parser.add_argument("-ff", type=str, default=["gaff", ], nargs="+",
                               help="Output, optional. Output format: .top \n"
                                    "The output file name of the generated gromacs \n"
                                    "style topology files. ")
    parser.parser.add_argument("-amberhome", type=str, default="/usr/local/amber/",
                               help="Input, optional. \n"
                                    "The AMBERHOME environment for topology generation. \n")
    parser.parser.add_argument("-addH", type=lambda x: (str(x).lower() == "true"), default=False,
                               help="Input, optional. Default is False. \n"
                                    "Whether add hydrogens before generating gmx topologies. \n")
    parser.parser.add_argument("-cplx", type=lambda x: (str(x).lower() == "true"), default=False,
                               help="Input, optional. Default is False. \n"
                                    "Whether the input is a ligand-receptor complex. ")
    parser.parser.add_argument("-ion", type=str, nargs="+", default=[],
                               help="The ions to be added to the system. Mutiple Ions supported. \n"
                                    "Options are Na+, Cl-\n"
                                    "Default is X+. \n")
    parser.parser.add_argument("-nion", type=int, nargs="+", default=[],
                               help="Number of ions to be added. Default is empty. \n"
                                    "Number of args must be the same as -ion. \n")
    parser.parser.add_argument("-neutral", type=lambda x: (str(x).lower() == "true"), default=False,
                               help="Input, optional. Default is False. "
                                    "Whether neutralize the system with Na+ and Cl- . ")
    parser.parser.add_argument("-bt", type=str, default="",
                               help="Box type, if you wish to solvate the mol in water box, you should\n"
                                    "provide an option: NA, TIP3PBOX, TIP4PEW, or TIP5PBOX. \n"
                                    "Default choice is empty. \n")
    parser.parser.add_argument("-size", type=float, default=0,
                               help="The size of your solvation box. Default is 0 angstrom. \n")

    parser.parse_arguments()
    args = parser.args

    print("AMBERHOME: ")
    print(args.amberhome)

    cmd = "export PATH=%s:%s:$PATH" % (args.amberhome, args.amberhome + "/bin")
    os.system(cmd)

    if args.addH:
        #fpdb = fixpdb.FixPDB()
        #fpdb.addhydrogenReduce(args.f, "H_"+args.f, args.amberhome+"/bin/reduce", flipH=False, verbose=True)
        pdbin = "H_"+args.f
        gtop = gentop.GenerateTop()
        try:
            gtop.addHydrogen(args.f, output=pdbin, method='rdkit')
        except RuntimeError:
            fpdb = fixpdb.FixPDB()
            fpdb.addhydrogenReduce(args.f, "H_" + args.f, args.amberhome + "/bin/reduce",
                                   flipH=False, verbose=True)
    else:
        pdbin = args.f

    # ligand prep and frcmod files
    ftype = "gaff" # amber gaff ff
    frcmod = "frcmod." + ftype
    prep   = "prep." + ftype

    gtop = gentop.GenerateTop()

    if not (os.path.exists(prep) and os.path.exists(frcmod)):
        # formal charge or net charge of a molecule
        spdb = fixpdb.SummaryPDB(pdbin)
        net_charge = spdb.netCharges(pdbin, "UNL")
        gtop.runAntechamber(pdbin, net_charge, ftype, amberhome=args.amberhome)

    gtop.gmxTopBuilder(
        PDBFile=[args.f, ], outputName=args.p, frcmodFile=frcmod,
        prepFile=prep, FField=args.ff, ionName=args.ion, ionNum=args.nion,
        solveBox=args.bt, boxEdge=args.size, amberhome=args.amberhome, verbose=True,
        neutral=args.neutral,
    )

    if len(args.bt) and args.cplx:
        if os.path.exists("topol_GMX.top"):
            tmp = open("tmp_top", 'w')
            with open("topol_GMX.top") as lines:
                for s in lines:
                    if "    1       IP      1          NA+         NA+       1      1     22.9898" in s:
                        tmp.write("    1       Na+      1          NA+         NA+       1      1     22.9898  \n")
                    elif "    1       IM      1         CL-           CL-      1     -1     35.45300" in s:
                        tmp.write("    1       Cl-      1         CL-           CL-      1     -1     35.45300  \n")
                    else:
                        tmp.write(s)
            tmp.close()
            os.remove("topol_GMX.top")
            os.rename("tmp_top", "topol_GMX.top")
        else:
            pass

