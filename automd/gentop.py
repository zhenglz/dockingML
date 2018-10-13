# -*- coding: utf-8 -*-

import sys, os
import subprocess as sp
import argparse
from argparse import RawDescriptionHelpFormatter, RawTextHelpFormatter
from dockml import convert
import time
from automd import fixpdb

class GenerateTop:
    """

    Parameters
    ----------

    Attributes
    ----------
    antechamber
    acpype

    """

    def __init__(self):
        PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
        self.antechamber = PROJECT_ROOT + "/../data/sample_antechamber.sh"

        # acpype for topology type convert
        self.acpype = "acpype"
        self.leapin = "leapIn.in"

    def gmxTopBuilder(self, PDBFile, outputName, frcmodFile=None,
                      prepFile=None, ionName=['Na'], ionNum=[0],
                      solveBox="TIP3PBOX", boxEdge=15,
                      FField=["AMBER99SB", ], amberhome="/app/amber16",
                      verbose=True):
        """
        provide a coordination file, output a amber and gromacs topology file.
        this is the entry point for generating a gmx topology.

        Parameters
        ----------
        PDBFile
        outputName: str, output,
            output topology file
        frcmodFile: str, input, optional
            not required if the molecule is a protein, or DNA, RNA
        prepFile: str, input
            prepare file, generated from tleap and antechamber
        ionName: list,
            the names of the ions
        ionNum: list,
            the number of ions to be added
        solveBox: str,
            the solvent water box, default is TIP3PBOX
        boxEdge: float,
            the edge of water box
        FField: list,
            the force field for topology generation
            a set of amber force fields, generally, we could choose
            99sb, 99sbildn, gaff, 12, 14sb
        amberhome: str,
            the AMBERHOME path
        verbose: bool,
            print details

        Returns
        -------

        """

        # check AMBERHOME PATH
        try:
            AMBERHOME = sp.check_output("echo $AMBERHOME", shell=True)
        except EnvironmentError:
            print("AMBERHOME is not defined, using the default instead.")
            AMBERHOME = amberhome
            #sys.exit(0)

        if verbose:
            print("Here is amber home path: ", AMBERHOME)

        if os.path.exists(self.leapin):
            os.remove(self.leapin)
        leap_contents = []

        # multiple FF supoorted here. Prepare a leap.in file for leap.
        leap_contents += self.writeFF2File(self.leapin, FField)

        # load amber frcmod and prep files
        if frcmodFile and prepFile:
            self.writePrepFrcmod2File(prepFile, frcmodFile)
        else:
            # prepare PDB file and load it
            if ".pdb" in PDBFile:
                leap_contents.append(("pdb = loadPDB  " + PDBFile + "  \n"))

            elif ".mol2" in PDBFile:
                # convert the mol2 file to pdb file using obabel
                convert.Convert(obabel="obabel").convert(PDBFile, PDBFile[:-4]+"pdb", verbose=verbose)
                leap_contents.append(("pdb = loadpdb " + PDBFile[:-4]+"pdb" + "  \n"))

        # input residue sequence to generate topology
        if len(PDBFile) >= 2 and ".pdb" not in PDBFile[0]:
            leap_contents += self.writeSequence2File(PDBFile)

        # add counter ions and solvate solute into water box
        self.writeIonInfo2File(ionName, ionNum)

        if solveBox and boxEdge:
            leap_contents.append("solvatebox pdb %s %f \n" % (solveBox, boxEdge))
        else:
            print("\nNot setting simulation box!\n")

        # save a amber off file of the molecule
        leap_contents += self.writeSaveParam2File(outputName)

        # write information to a leapin file
        with open(self.leapin, 'w') as tofile:
            for s in leap_contents:
                tofile.writable(s)

        if verbose:
            print("Generating a leap input file for tleap topology processing.")

        # run tleap
        try:
            out = sp.check_output("%s/bin/tleap -f %s  \n" % (AMBERHOME, self.leapin), shell=True)
            print("Generating amber topology files! ")
            if verbose:
                print(out)
        except SystemExit:
            print("tleap loading failed! Exit now!")
            sys.exit(0)

        # convert AMBER format to GMX format
        time.sleep(2)
        if len(self.acpype):
            try:
                out = sp.check_output("%s -b %s -x %s.prmcrd -p %s.prmtop "
                                      %(self.acpype, outputName, outputName, outputName)
                                      )
                if verbose:
                    print(out)
            except SystemExit:
                print("Converting AMBER files using ACPYPE to GMX failed! ")
        else:
            try:
                out = sp.check_output("Amb2gmx.pl --prmtop %s.prmtop --crd %s.prmcrd --outname gmx_%s "
                                      % (outputName, outputName, outputName), shell=True)
                if verbose:
                    print(out.decode() + "\n\nGMX and Amber topologies created!")

            except SystemExit:
                print("Converting AMBER files using AMB2GMX to GMX failed!")

        return None

    def writeFF2File(self, filename, ff_list):
        leapin = []
        for ff in ff_list :
            # create tleap input file
            if "gaff" in ff:
                leapin.append("source leaprc.gaff \n")
            elif "ildn" in ff or "ILDN" in ff:
                leapin.append("source oldff/leaprc.ff99SBildn  \n")
            elif ff == "AMBER99SB" or "amber99sb" == ff or "99sb" in ff.lower() :
                leapin.append("source oldff/leaprc.ff99SB  \n")
            elif ff == "AMBER14SB" or "14" in ff:
                leapin.append("source leaprc.ff14SB  \n")
            elif "leaprc." in ff:
                leapin.append("source %s  \n" % ff)
            else:
                print("ff not in amber ff library. Pass")

        return leapin

    def writePrepFrcmod2File(self, prepFile, frcmodFile):
        leap_contents = list()
        leap_contents.append(("loadamberprep " + prepFile + "  \n"))
        leap_contents.append(("loadamberparams  " + frcmodFile + "  \n"))
        leap_contents.append("pdb = LIG \n")

        return leap_contents

    def writeSequence2File(self, seqs):
        leap_contents = list()
        leap_contents.append("pdb = sequence{ ")
        for item in seqs:
            leap_contents.append((item + " "))
            leap_contents.append(" } \n")

        return leap_contents

    def writeSaveParam2File(self, outputName):
        leapin = list()
        leapin.append("check pdb \n")
        leapin.append("savepdb pdb %s.pdb  \n" % outputName)
        leapin.append("saveoff pdb %s.lib  \n"%outputName)
        leapin.append("saveamberparm pdb %s.prmtop %s.prmcrd \n" % (outputName, outputName))
        leapin.append("quit \n")

        return leapin

    def writeIonInfo2File(self, ionnames, ionnum):
        leapin = list()
        if len(ionnames) == len(ionnum) and len(ionnames):
            for i in range(len(ionnum)):
                leapin.append("addions2 pdb %s %d \n" % (ionnames[i], ionnum[i]))
        else:
            print("\nNot writing any ion information\n")

        return leapin


    def pdb2gmx(self, pdb2gmx, pdbin,
                groout, topolout,
                protein=6, tip3p=1,
                verbose=True):
        """

        Parameters
        ----------
        pdb2gmx: str,
            system gmx pdb2gmx commond
        pdbin: str,
            input pdb file name
        groout: str,
            output gmx gro file name
        topolout: str,
            output gmx top file name
        protein: bool,
            this parameter is not used
        tip3p: str,
            this parameter is not used
        verbose: bool,
            print details

        Returns
        -------

        """
        cmd = "%s -f %s -o %s -p %s -ignh" % (pdb2gmx, pdbin, groout, topolout)
        out = sp.check_output(cmd, shell=True)
        if verbose :
            print(out)
        #job.communicate(input="%d \n %d \n"%(protein, tip3p))

        return None

    def addWatIon(self, editconf, genbox,
                  genion, grompp, top,
                  groin='in.gro', distance=1.2,
                  conc=0.15, spc="spc903.gro",
                  mdpfile = "em_sol.mdp",
                  verbose=True
                  ):
        """
        add water and ions

        Parameters
        ----------
        editconf
        genbox
        genion
        grompp
        top
        groin
        distance
        conc
        spc
        mdpfile
        verbose

        Returns
        -------

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

        return None

    def top2itp(self, outputName, topFileName="", verbose=True):
        """
        generate a .itp file from the .top file

        Parameters
        ----------
        outputName
        topFileName
        verbose

        Returns
        -------

        """
        # generate itp file from gmx.top file
        itpfile = open(outputName + ".itp", 'w')
        if not topFileName:
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
        if verbose:
            print("ITP file for %s generated! " % outputName)
        return None

    def runAntechamber(self, infile, netCharge=None):
        """
        run antechamber to generate RESP am1/bcc atomic charges

        Parameters
        ----------
        infile
        netCharge

        Returns
        -------

        """

        antechamber = self.antechamber

        if not os.path.exists(antechamber):
            print("Sample antachamber file not found. Exit now!")
            sys.exit(1)

        # run antechamber here, may first generate a sh script and then wrap it up
        if not bool(netCharge):
            try:
                spdb = fixpdb.SummaryPDB(infile, "")
                netcharge = spdb.netCharges(inputMol=infile)
            except :
                print("Getting netcharge error! Using default number 0!")
                netcharge = 0

        tofile = open("antechamber.sh", 'w')
        with open(antechamber) as lines:
            for s in lines:
                if len(s.split()) > 0 and s.split()[0] == "antechamber":
                    tofile.write("antechamber -i $1 -fi mol2 -o prep.$2 -fo prepi -at $2 -pf y -s 2 -c bcc -nc %d \n"
                                 % netCharge)
                else:
                    tofile.write(s)
        tofile.close()

        # run the antechamber script
        try:
            job = sp.Popen("sh ./antechamber.sh", shell=True)
            job.communicate()
        except SystemExit :
            print("Run antachamber failed.")

        return None

    def trimHydrogen(self, reduce, pdbin, pdbout, verbose=False) :
        """
        remove hydrogens in pdb file. call amber reduce function.

        Parameters
        ----------
        reduce
        pdbin
        pdbout
        verbose

        Returns
        -------

        """

        job = sp.check_output('%s -Trim %s > %s '%(reduce, pdbin, pdbout), shell=True)

        if verbose:
            print(job)

        return None

    def arguments(self):
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
        if len(sys.argv) < 3:
            parser.print_help()
            print("Number of arguments are not correct! Exit Now!")
            sys.exit(1)

        return args
