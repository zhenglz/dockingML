#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
import subprocess as sp
import argparse
from argparse import RawTextHelpFormatter
from automd.utils.gentop import GenerateTop, AcpypeGenTop
import mdtraj as mt


def arguments():
    """
    Parse arguments in the command line mode.

    GMX-AMBER, topology builder.\n
    Provide a PDB file, or sometimes prep and frcmod files, itp and prmtop file will be given.
    Copyright @ Liangzhen Zheng, contact astrozheng@gmail.com for any technique support. \n

    Examples
    --------
    Show help information
      python autoMD.py gentop -h

    Generate Amberff14SB for a protein molecule
      python autoMD.py gentop -inp your.pdb -out yourpdb -ff ff14

    Generate Amberff14SB for peptide sequences
      python autoMD.py gentop -aaseq ACE ALA CYS ALA HIS ASN NME -ff AMBER99SBildn -out peptide


      For small molecules,
      python autoMD.py gentop -inp ligand.mol2 -resp True -out LIG
      python autoMD.py gentop -inp lig.pdb -out lig -ff gaff -prep amber.prep.lig -frcmod frcmod.log
      python autoMD.py gentop -inp cplx.pdb -out cplx_box -ff gaff ildn
                  -bt TIP3PBOX -d 1.2 -prep amber.prep.lig -frcmod frcmod.log

    Returns
    -------
    args: ArgumentParser object
    """

    # go current directory, pwd
    os.chdir(os.getcwd())

    d = '''
    GMX-AMBER, topology builder.\n
    Provide a PDB file, or sometimes prep and frcmod files, itp and prmtop file will be given.
    Copyright @ Liangzhen Zheng, contact astrozheng@gmail.com for any technique support. \n
    Examples:
    Show help information
      python autoMD.py gentop -h

    Generate Amberff14SB for a protein molecule
      python autoMD.py gentop -inp your.pdb -out yourpdb -ff ff14

    Generate Amberff14SB for peptide sequences
      python autoMD.py gentop -aaseq ACE ALA CYS ALA HIS ASN NME -ff AMBER99SBildn -out peptide


      For small molecules,
      python autoMD.py gentop -inp ligand.mol2 -resp True -out LIG
      python autoMD.py gentop -inp lig.pdb -out lig -ff gaff -prep amber.prep.lig -frcmod frcmod.log
      python autoMD.py gentop -inp cplx.pdb -out cplx_box -ff gaff ildn
                  -bt TIP3PBOX -d 1.2 -prep amber.prep.lig -frcmod frcmod.log
    '''
    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-inp", type=str, default="protein_xxx.pdb",
                        help="The PDB file or a mol2 file for protein topology generating. "
                             "Default is None\n")
    parser.add_argument("-lig", type=str, default="xxx_ligand.mol2",
                        help="The ligand file (mol2, pdb, pdbqt, smi) for protein topology generating. "
                             "Default is None\n")
    parser.add_argument("-out", type=str, default="OUT", help="The output name. Default is OUT.\n")
    parser.add_argument("-ion", type=str, nargs="+", default=["X+", ],
                        help="The ions to be added to the system. Mutiple Ions supported. \n"
                             "Options are Na+, Cl-\n"
                             "Default is X+. \n")
    parser.add_argument("-nion", type=int, nargs="+", default=[-1, ],
                        help="Number of ions to be added. Default is -1. \n"
                             "Number of args must be the same as -ion. \n")
    parser.add_argument("-bt", type=str, default="TIP3PBOX",
                        help="Box type, ff you wish to solvate the mol in water box, you should\n"
                             "provide an option: NA, TIP3PBOX, TIP4PEW, or TIP5PBOX. \n"
                             "Default choice is TIP3PBOX. \n")
    parser.add_argument("-size", type=float, default=10,
                        help="The size of your solvation box. Default is 10 angstrom. \n")
    parser.add_argument("-ff", type=str, nargs='+', default=["AMBER99SB", ],
                        help="The force field for simulation. Multiple force filed files were supported.\n"
                             "Options including: 99SB, 99SBildn, ff14, gaff\n"
                             "or use leaprc.xxx files here. \n"
                             "For protein, default is AMBER99SB. \n"
                             "For ligand, default is gaff. \n")
    parser.add_argument("-aaseq", type=str, nargs='+', default=None,
                        help="Amino acid sequences in case no PDB file provide.\n"
                             "Default is None. It is an optional argument.\n")
    parser.add_argument("-c", default='bcc', type=str,
                        help="If it is a small molecule, whose atomic charges \n"
                             "not defined, neither the bonding angles, we could \n"
                             "use RESP with bcc AM1 to determine charges and \n"
                             "prepare parameters. \n"
                             "Defualt is bcc. \n")
    parser.add_argument("-nc", default=None, type=int,
                        help="This argument only works with \n"
                             "the -c argument. Default is 0. \n")
    parser.add_argument("-H", default="~/anaconda3", type=str,
                        help="AMBERHOME path. Default is ~/anaconda3.")

    # args include the options in the argparser
    args, unknown = parser.parse_known_args()
    if len(sys.argv) < 3:
        parser.print_help()
        print("\n\nNumber of arguments are not correct! Exit Now!\n\n")
        sys.exit(1)

    return args


def runGenTop(AMBERHOME="") :
    '''Create a gmxTopBuilder class, and input some parameters,
        to generate amber and gromacs topology and coordinates

    if a small ligand provide, if not prep and frcmod exist,
        then am1/bcc charge based on amber antechamber sqm method
        will be calculated. in this process, tleap will be used to
        get bonded and non-bonded parameters

    '''

    top = GenerateTop()
    args = arguments()

    if AMBERHOME == "":
        AMBERHOME = args.H

    if os.path.isdir(AMBERHOME) and os.path.exists(os.path.join(AMBERHOME, "bin/tleap")):
        try:
            os.environ["AMBERHOME"] = AMBERHOME
        except:
            os.system("export AMBERHOME=%s"%AMBERHOME)
    else:
        print("INFO: AMBERHOME %s not set or not found " % AMBERHOME)
        sys.exit(0)

    # define the input coordinates information, a pdb, or mol2, or a sequence
    # of amino acids
    if not os.path.exists(args.lig):
        if args.aaseq:
            structure = args.aaseq
        else:
            structure = [args.inp, ]

        # not calculate atomic charges
        top.gmxTopBuilder(structure, args.out, amberhome=AMBERHOME,
                          ionName=args.ion, ionNum=args.nion,
                          solveBox=args.bt, boxEdge=args.size,
                          FField=args.ff)
    else:
        print("INFO: processing ligand %s" % args.lig)
        # prepare resp charges for small ligand with acpype
        lig_top = AcpypeGenTop(args.lig)
        #lig_top.run_acpype("ligand")

        frcmod = os.path.join("ligand.acpype", "ligand_AC.frcmod")
        off = os.path.join("ligand.acpype", "ligand_AC.lib")
        lig_pdb = os.path.join("ligand.acpype", "ligand_NEW.pdb")

        if os.path.exists(args.inp):
            # keep only protein
            pdb = mt.load_pdb(args.inp)
            protein_ndx = pdb.topology.select("protein")
            pdb = pdb.atom_slice(protein_ndx)
            pdb.save_pdb(os.path.basename(args.inp)+"_temp.pdb")

            # remove hydrogen
            cmd = "reduce -Trim %s > %s" %(os.path.basename(args.inp)+"_temp.pdb", os.path.basename(args.inp)+"_noHy.pdb")
            job = sp.Popen(cmd, shell=True)
            job.communicate()

            # make a complex pdb
            cmd = "cat %s %s | grep \"ATOM\" > %s" % (os.path.basename(args.inp)+"_noHy.pdb", lig_pdb,
                                                      os.path.basename(args.inp)+"_complex.pdb")
            job = sp.Popen(cmd, shell=True)
            job.communicate()

            structure = [os.path.basename(args.inp)+"_complex.pdb", ]

        else:
            structure = [lig_pdb, ]

        top.gmxTopBuilder(structure, args.out,
                          frcmodFile=frcmod,
                          prepFile="",
                          offFile=off, amberhome=AMBERHOME,
                          ionName=args.ion, ionNum=args.nion,
                          solveBox=args.bt, boxEdge=args.size,
                          FField=args.ff)


if __name__ == "__main__":

    runGenTop()

