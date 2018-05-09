# -*- coding: utf-8 -*-

import sys, os
from datetime import time
import subprocess as sp
from automd import sumpdb
import argparse
from argparse import RawDescriptionHelpFormatter, RawTextHelpFormatter

class GenerateTop :
    def __init__(self):
        PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
        self.antechamber = PROJECT_ROOT + "/../data/sample_antechamber.sh"

    def gmxTopBuilder(self, PDBFile, outputName, frcmodFile=None,
                      prepFile=None, ionName=[''], ionNum=[0],
                      solveBox=None, boxEdge=12,
                      FField=["AMBER99SB", ], ACPYPE="acpype",
                      verbose=True):
        '''
        provide a coordination file, output a amber and gromacs topology file.

        required parameters

        :param frcmodFile: not required if the molecule is a protein, or DNA, RNA
        :param prepfile: not required if the molecule is a protein, or DNA, RNA
        :param PDBFile: input, coordination file
        :param outputName: output, output file name
        :param ionName: optional
        :param ionNum: optional
        :param solveBox: optional, type of solvation box, default is TIP3PBOX
        :param boxEdge: optional
        :param FField: a set of amber force fields, generally, we could choose
            99sb, 99sbildn, gaff, 12, 14sb
        :param verbose:
        :return: None
        '''

        # check AMBERHOME PATH
        AMBERHOME = sp.check_output("echo $AMBERHOME", shell=True)
        #AMBERHOME = AMBERHOME[:-1] + "/"
        if verbose :
            print("Here is amber home path:")
            print(AMBERHOME)

        AMBERHOME = "/app/amber16"

        # multiple FF supoorted here
        leapin = open("leapIn.in", 'wb')
        for ff in FField:
            # create tleap input file
            if "gaff" in ff:
                leapin.write("source leaprc.gaff \n".encode())
            elif "ildn" in ff or "ILDN" in ff:
                leapin.write("source oldff/leaprc.ff99SBildn  \n".encode())
            elif ff == "AMBER99SB" or "amber99sb" == ff or "99sb" in ff.lower() :
                leapin.write("source oldff/leaprc.ff99SB  \n".encode())
            elif ff == "AMBER14SB" or "14" in ff:
                leapin.write("source leaprc.ff14SB  \n".encode())
            elif "leaprc." in ff:
                leapin.write(("source %s  \n" % ff).encode())
            else:
                print( "Load Force Field File Error! \nExit Now!")
                sys.exit(1)

        # load amber frcmod and prep files
        if frcmodFile :
            leapin.write(("loadamberparams  " + frcmodFile + "  \n").encode())
        if prepFile :
            leapin.write(("loadamberprep " + prepFile + "  \n").encode())

        if frcmodFile and prepFile :
            leapin.write("pdb = LIG \n".encode())
        else :
            # prepare PDB file and load it
            if ".pdb" in PDBFile:
                leapin.write(("pdb = loadPDB  " + PDBFile + "  \n").encode())
            elif ".mol2" in PDBFile :
                # convert the mol2 file to pdb file using obabel
                job = sp.Popen("obabel %s -O %s"%(PDBFile, PDBFile[:-4]+"pdb"), shell=True)
                job.communicate()
                leapin.write(("pdb = loadpdb " + PDBFile[:-4]+"pdb" + "  \n").encode())
            elif len(PDBFile) >= 2 and ".pdb" not in PDBFile[0] :
                leapin.write("pdb = sequence{ ".encode())
                for item in PDBFile:
                    leapin.write((item + " ").encode())
                leapin.write(" } \n".encode())
            else:
                print( "Loading PDB file or Sequence file error!")
                sys.exit(1)

        # save a amber off file of the molecule
        leapin.write(("saveoff pdb %s.lib \n" % outputName).encode())

        # add counter ions and solvate solute into water box
        if ionName and ionNum and ionNum[0] > 0 and ionName[0] != "X+":
            if len(ionNum) == len(ionName):
                for i in range(len(ionNum)):
                    leapin.write(("addions2 pdb %s %d \n" % (ionName[i], ionNum[i])).encode())
            else:
                print( "\nAdd ions not successful!\n")
        else:
            "\nNot adding any ions!\n"
        if solveBox :
            if boxEdge :
                leapin.write(("solvatebox pdb %s %f \n" % (solveBox, boxEdge)).encode())
            else:
                print( "\nBOX size not correctly set.\nExit Now!\n")
                sys.exit(1)
        else:
            print( "\nNot setting simulation box!\n")

        # check object
        leapin.write("check pdb \n".encode())
        leapin.write(("savepdb pdb %s  \n" % (outputName + ".pdb")).encode())
        leapin.write(("saveoff pdb %s  \n"%(outputName+".lib")).encode())
        leapin.write(("saveamberparm pdb %s.prmtop %s.prmcrd \n" % (outputName, outputName)).encode())
        leapin.write("quit \n".encode())
        leapin.close()

        if verbose :
            print("Generating a leap input file for tleap topology processing.")

        # run tleap
        try :
            out = sp.check_output("%s/bin/tleap -f leapIn.in  \n" % AMBERHOME, shell=True)

            if verbose:
                print(out.decode())

        except :
            print("tleap loading failed!")

        # convert AMBER format to GMX format
        #time.sleep(2)
        if len(ACPYPE) :
            try :
                out = sp.check_output("%s -b %s -x %s.prmcrd -p %s.prmtop "
                                      %(ACPYPE, outputName, outputName, outputName)
                                      )
            except :
                "Converting AMBER files using ACPYPE to GMX failed! "

        else :
            try:
                out = sp.check_output("Amb2gmx.pl --prmtop %s.prmtop --crd %s.prmcrd --outname gmx_%s " \
                                      % (outputName, outputName, outputName), shell=True)
                if verbose:
                    print(out.decode() + "\n\nGMX and Amber topologies created!")

            except :
                print("Converting AMBER files using AMB2GMX to GMX failed!")

        return 1

    def pdb2gmx(self, pdb2gmx, pdbin,
                groout, topolout,
                protein=6, tip3p=1,
                verbose=True):
        """
        call system GMX pdb2gmx function
        :param pdb2gmx:
        :param pdbin:
        :param groout:
        :param topolout:
        :param verbose:
        :return:
        """
        cmd = "%s -f %s -o %s -p %s -ignh" % (pdb2gmx, pdbin, groout, topolout)
        job = sp.Popen(cmd, shell=True)
        job.communicate(input="%d \n %d \n"%(protein, tip3p))

        return 1

    def addWatIon(self, editconf, genbox,
                  genion, grompp, top,
                  groin='in.gro', distance=1.2,
                  conc=0.15, spc="spc903.gro",
                  mdpfile = "em_sol.mdp",
                  verbose=True
                  ):
        """
        add water and ions
        :param editconf:
        :param genbox:
        :param genion:
        :param grompp:
        :param top:
        :param groin:
        :param distance:
        :param conc:
        :param spc:
        :param mdpfile:
        :param verbose:
        :return:
        """

        # genbox
        cmd1 = "%s -f %s -o %s -d %f -c " % (editconf, groin, "box_" + groin, distance)
        job = sp.Popen(cmd1, shell=True)
        job.communicate()
        cmd2 = "%s -cs %s -cp %s -o %s -p %s " % (genbox, spc, "box_" + groin, "wat_" + groin, top)
        job = sp.Popen(cmd2, shell=True)
        job.communicate()

        # genion
        cmd3 = "%s -f %s -p %s -c %s -o addion -maxwarn 50 " %( grompp, mdpfile, top, "wat_"+groin)
        job = sp.Popen(cmd3,shell=True)
        job.communicate()

        cmd4 = "%s -s addion -p %s -o %s -neutral -conc %f " % (genion, top, "ion_"+groin, conc)
        job = sp.Popen(cmd4, shell=True)
        job.communicate()

        return 1

    def runObabel(self, obabelexe, input, output, verbose=True):
        """
        run openbabel
        :param obabelexe:
        :param input:
        :param output:
        :param verbose:
        :return:
        """
        if os.path.exists(input) :
            status = sp.check_output("%s %s -O %s "%(obabelexe, input, output), shell=True)
            if verbose :
                print(status)
        else :
            print( input + " is not exist!!!")
            sys.exit(0)

        return(1)

    def top2itp(self, outputName, topFileName=None, verbose=True):
        """
        prepare a itp file from top file
        :param outputName:
        :param topFileName:
        :param verbose:
        :return:
        """
        # generate itp file from gmx.top file
        itpfile = open(outputName + ".itp", 'w')
        if not topFileName :
            topFileName = "gmx_" + outputName + ".top"

        with open(topFileName) as lines:
            start = 0
            stop = 0
            for s in lines:
                if "[ moleculetype" in s:
                    start = 1
                if "[ system" in s:
                    stop = 1

                if start and not stop:
                    if "solute" in s:
                        itpfile.write("%-6s            3 \n" % outputName)
                    else:
                        itpfile.write(s)
        itpfile.close()
        if verbose :
            print("ITP file for %s generated! " % outputName)
        return(1)

    def runAntechamber(self, infile, netCharge=None):
        '''
        run antechamber to generate RESP am1/bcc atomic charges
        :param netcharge: input, number of netcharges of the molecule
        :param antechamber: input, optional, a sample antechamber shell script
        :return: the file name of the antechamber shell script
        '''

        antechamber = self.antechamber

        if not os.path.exists(antechamber) :
            with open(antechamber, 'w') as tofile :
                content = '''
if [ $# -ne 2 ]; then
echo "Usage: input RED III.1  output, charge & spin & residuNAME read from Mol2"
echo "Format: file_in(Mol2); atom_type(gaff or amber)"
exit
fi

antechamber -i $1 -fi mol2 -o prep.$2 -fo prepi -at $2 -pf y -s 2
#antechamber -i $1 -fi mol2 -o prep.$2 -fo prepi -at $2 -pf y -s 2 -c resp
parmchk -i prep.$2 -f prepi -o frcmod.$2 -a Y
grep need frcmod.$2

if test -f leap.in; then
   rm leap.in
fi

echo "============================================="
echo "dmfff = loadamberparams frcmod.$2" >> leap.in
echo "loadamberprep prep.$2" >> leap.in

echo "prepared run_parm.sh for tleap"
echo "tleap -s -f ${AMBERHOME}dat/leap/cmd/leaprc.ff99 -f leap.in"
echo "tleap -s -f ${AMBERHOME}dat/leap/cmd/leaprc.gaff -f leap.in"
echo "general AMBER : leaprc.gaff"
echo "AMBER : leaprc.ff90"
                '''
                tofile.write(content)
        # run antechamber here, may first generate a sh script and then wrap it up
        if netCharge is None :
            try:
                spdb = sumpdb.SummaryPDB(pdbfile=infile)
                netcharge = spdb.netCharges(inputMol=infile)
            except :
                print("Getting netcharge error!")
                netcharge = 0
        else :
            netcharge = netCharge

        tofile = open("antechamber.sh", 'w')
        with open(antechamber) as lines:
            for s in lines:
                if len(s.split()) > 0 and s.split()[0] == "antechamber":
                    tofile.write(
                        "antechamber -i $1 -fi mol2 -o prep.$2 -fo prepi -at $2 -pf y -s 2 -c bcc -nc %d \n" % netcharge)
                else:
                    tofile.write(s)
        return "antechamber.sh"

    def trimHydrogen(self, reduce, pdbin, pdbout, verbose=False) :
        """
        remove hydrogens in pdb file
        :param reduce:
        :param pdbin:
        :param pdbout:
        :param verbose:
        :return:
        """
        job = sp.check_output('%s -Trim %s > %s '%(reduce, pdbin, pdbout), shell=True)

        if verbose :
            print(job)

        return 1

    def arguments(self):
        '''
        Parse arguments in the command line mode.

        :return: args, a parser.argument object
        '''

        # go current directory, pwd
        os.chdir(os.getcwd())

        if 0 :
            frcmodFile = "fit_R2.frcmod"
            PDBFile = "trpcage.pdb"
            outputName = "TRP"
            FField = "AMBER99SB"
            self.gmxTopBuilder(frcmodFile, PDBFile, outputName, FField)

        else:
            d = '''
            GMX-AMBER, topology builder.\n
            Provide a PDB file, or sometimes prep and frcmod files, itp and prmtop file will be given.
            Copyright @ Liangzhen Zheng, contact astrozheng@gmail.com for any technique support. \n
            Examples:
              python autoMD.py gentop -h
              python autoMD.py gentop -inp your.pdb -out yourpdb -ff ff14
              python autoMD.py gentop -aaseq ACE ALA CYS ALA HIS ASN NME -ff AMBER99SBildn -out peptide
              python autoMD.py gentop -inp lig.pdb -out lig -ff gaff -prep amber.prep.lig -frcmod frcmod.log
              python autoMD.py gentop -inp cplx.pdb -out cplx_box -ff gaff ildn
                                      -bt TIP3PBOX -d 1.2 -prep amber.prep.lig -frcmod frcmod.log
              python autoMD.py gentop -inp ligand.mol2 -resp True -nt 0 -out LIG
            '''
            parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
            parser.add_argument("-inp", type=str, default=None,
                                help="The PDB file or a mol2 file for topology generating. Default is None\n")
            parser.add_argument("-prep", type=str, default="amberff.prep",
                                help="Prep file (amber format) stored coordinates and charges.\n"
                                     "Default is None.")
            parser.add_argument("-frcmod", type=str, default="frcmod",
                                help="The additional parameters, stored in frcmod file (Amber format).\n")
            parser.add_argument("-out", type=str, default="OUT", help="The output name. Default is OUT.\n")
            parser.add_argument("-ion", type=str, nargs="+", default=["X+", ],
                                help="The ions to be added to the system. Mutiple Ions supported. \n"
                                     "Options are Na+, Cl-\n"
                                     "Default is X+. \n")
            parser.add_argument("-nion", type=int, nargs="+", default=[-1, ],
                                help="Number of ions to be added. Default is -1. \n"
                                     "Number of args must be the same as -ion. \n")
            parser.add_argument("-bt", type=str, default=None,
                                help="Box type, ff you wish to solvate the mol in water box, you should\n"
                                     "provide an option: NA, TIP3PBOX, TIP4PEW, or TIP5PBOX. \n"
                                     "Default choice is TIP3PBOX. \n")
            parser.add_argument("-size", type=float, default=0,
                                help="The size of your solvation box. Default is 0 angstrom. \n")
            parser.add_argument("-ff", type=str, nargs='+', default=["AMBER99SB", ],
                                help="The force field for simulation. Multiple force filed files were supported.\n"
                                     "Options including: 99SB, 99SBildn, ff14, gaff\n"
                                     "or use leaprc.xxx files here. Default is AMBER99SB. \n")
            parser.add_argument("-aaseq", type=str, nargs='+', default=None,
                                help="Amino acid sequences in case no PDB file provide.\n"
                                     "Default is None. It is an optional argument.\n")
            parser.add_argument("-resp", default=False, type=bool,
                                help="If it is a small molecule, whose atomic charges \n"
                                     "not defined, neither the bonding angles, we could \n"
                                     "use RESP with bcc AM1 to determine charges and \n"
                                     "prepare parameters. Options: True, False. \n"
                                     "Defualt is False. \n")
            parser.add_argument("-nc", default=None, type=int,
                                help="This argument only works with \n"
                                     "the -resp argument. Default is 0. \n")

            # args include the options in the argparser
            args, unknown = parser.parse_known_args()
            if len(sys.argv) < 3 :
                parser.print_help()
                print("Number of arguments are not correct! Exit Now!")
                sys.exit(1)

            return(args)