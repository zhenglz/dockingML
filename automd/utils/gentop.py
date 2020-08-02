# -*- coding: utf-8 -*-

import sys, os
import subprocess as sp
from dockml import convert
import time
try:
    from rdkit import Chem
except:
    print("INFO: rdkit not found ...")
from automd.utils import fixpdb


class PrepareLigand(object):

    def __init__(self, ligand="ligand.mol2", output="output.mol2"):

        self.ligand = ligand
        self.output = output

    def _addHydrogen(self, method="rdkit"):
        """
        Add hydrogens to a molecule

        Parameters
        ----------
        pdbfile: str, input
            the pdb file for hydrogen adding
        output: str, output
            the output file name
        method: str, default is rdkit
            the hydrogen adding method

        Returns
        -------

        """

        if method == "rdkit":
            mol = Chem.MolFromPDBFile(self.ligand)
            mol2 = Chem.AddHs(mol, addCoords=True)

            writer = Chem.PDBWriter(self.output)
            writer.write(mol2)

            writer.close()

        else:
            cmd = "reduce %s > %s " % (self.ligand, self.output)

            job = sp.Popen(cmd, shell=True)
            job.communicate()

        return 0


class GenerateTop(object):
    """

    Parameters
    ----------

    Attributes
    ----------
    PROJECT_ROOT

    antechamber: str,
        the sample antechamber shell script
    acpype: str,
        the executable acpype commond, for topology format coverting
    leapin: str, default is leapIn.in
        the output file, feeding it for tleap to generating topologies

    """

    def __init__(self):
        self.PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
        self.antechamber = self.PROJECT_ROOT + "/../data/sample_antechamber.sh"

        # acpype for topology type convert
        self.acpype = "acpype"
        self.leapin = "leapIn.in"

    def gmxTopBuilder(self, PDBFile, outputName, frcmodFile=None,
                      prepFile=None, offFile=None, ionName=['Na', ], ionNum=[0, ],
                      solveBox="TIP3PBOX", boxEdge=15,
                      FField=["AMBER99SB", ], amberhome="/app/amber16",
                      neutral=True,
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
        AMBERHOME = amberhome

        if verbose:
            print("Here is amber home path: ", AMBERHOME)

        if os.path.exists(self.leapin):
            os.remove(self.leapin)
        leap_contents = []

        # multiple FF supoorted here. Prepare a leap.in file for leap.
        leap_contents += self.writeFF2File(self.leapin, FField)

        # load amber frcmod and prep files
        if os.path.exists(frcmodFile) or os.path.exists(prepFile) \
                or os.path.exists(offFile):
            leap_contents += self.writePrepFrcmod2File(prepFile, frcmodFile, offFile)

        # prepare PDB file and load it
        if ".pdb" in PDBFile[0]:
            leap_contents.append(("pdb = loadPDB  " + PDBFile[0] + "  \n"))

        elif ".mol2" in PDBFile[0]:
            # convert the mol2 file to pdb file using obabel
            convert.Convert(obabel="obabel").convert(PDBFile[0], PDBFile[0][:-4]+"pdb", verbose=verbose)
            leap_contents.append(("pdb = loadpdb " + PDBFile[0][:-4]+"pdb" + "  \n"))

        # input residue sequence to generate topology
        if ".pdb" not in PDBFile[0] and len(PDBFile) > 2 :
            leap_contents += self.writeSequence2File(PDBFile)

        # add counter ions and solvate solute into water box
        leap_contents += self.writeIonInfo2File(ionName, ionNum)

        if neutral:
            leap_contents.append("addions pdb Na+ 0 \n")
            leap_contents.append("addions pdb Cl- 0 \n")

        if len(solveBox) and boxEdge > 0 :
            leap_contents.append("solvatebox pdb %s %f \n" % (solveBox, boxEdge))
        else:
            print("\nNot setting simulation box!\n")

        # save a amber off file of the molecule
        leap_contents += self.writeSaveParam2File(outputName)

        # write information to a leapin file
        with open(self.leapin, 'w') as tofile:
            for s in leap_contents:
                tofile.write(s)

        if verbose:
            print("Generating a leap input file for tleap topology processing.")

        # run tleap
        try:
            job = sp.Popen("%s/bin/tleap -f %s  \n" % (AMBERHOME, self.leapin), shell=True)
            print("Generating amber topology files! ")
            job.communicate()
        except SystemError:
            print("tleap loading failed! Exit now!")
            sys.exit(0)

        # convert AMBER format to GMX format
        time.sleep(2)
        if len(self.acpype):
            try:
                out = sp.Popen("%s -b %s -x %s.inpcrd -p %s.prmtop "%
                               (self.acpype, outputName, outputName, outputName),
                               shell=True)
                out.communicate()
            except SystemError:
                print("Converting AMBER files using ACPYPE to GMX failed! ")
        else:
            try:
                out = sp.check_output("Amb2gmx.pl --prmtop %s.prmtop --crd %s.prmcrd --outname gmx_%s "
                                      % (outputName, outputName, outputName), shell=True)
                if verbose:
                    print(out.decode() + "\n\nGMX and Amber topologies created!")

            except SystemError:
                print("Converting AMBER files using AMB2GMX to GMX failed!")

        return None

    def writeFF2File(self, filename, ff_list):
        leapin = []
        for ff in ff_list :
            # create tleap input file
            if "gaff" in ff:
                leapin.append("source leaprc.gaff2 \n")
            elif "ildn" in ff or "ILDN" in ff:
                leapin.append("source oldff/leaprc.ff99SBildn  \n")
            elif ff == "AMBER99SB" or "amber99sb" == ff or "99sb" in ff.lower() :
                leapin.append("source oldff/leaprc.ff99SB  \n")
            elif ff == "AMBER14SB" or "14" in ff:
                leapin.append("source leaprc.protein.ff14SB  \n")
            elif "leaprc." in ff:
                leapin.append("source %s  \n" % ff)
            else:
                print("ff not in amber ff library. Pass")

            leapin.append("source leaprc.water.tip3p\n")
        return leapin

    def writePrepFrcmod2File(self, prepFile, frcmodFile, offFile):
        leap_contents = list()
        if os.path.exists(prepFile):
            leap_contents.append(("loadamberprep " + prepFile + "  \n"))
        if os.path.exists(frcmodFile):
            leap_contents.append(("loadamberparams  " + frcmodFile + "  \n"))
        if os.path.exists(offFile):
            leap_contents.append("loadoff %s \n" % offFile)

        #leap_contents.append("pdb = %s \n" % ligand)

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
        leapin.append("saveamberparm pdb %s.prmtop %s.inpcrd \n" % (outputName, outputName))
        leapin.append("quit \n")

        return leapin

    def writeIonInfo2File(self, ionnames, ionnum):
        leapin = list()
        leapin.append("source leaprc.water.tip3p \n")
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
            system gmx pdb2gmx command
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

    def top2itp(self, outputName, topFileName=None, verbose=True):
        """
        generate a .itp file from the .top file, the [ system ] part is
        removed

        Parameters
        ----------
        outputName: str,
            output itp file
        topFileName: str,
            input topology file
        verbose: bool, default is True
            whether print detail information

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

    def runAntechamber(self, infile, netCharge=None, filetype="gaff", amberhome="/app/amber16"):
        """
        run antechamber to generate RESP am1/bcc atomic charges, meanwhile the gaff
        parms for bonding and angles are also generated in prep and frcmod files
        details could be found here: ambermd.org/doc12/Amber16.pdf

        Parameters
        ----------
        infile: str,
            input pdb file or mol2 file, which is used for sqm calculation for am1/bcc
            charges
        netCharge: int,
            the netcharge of the molecule, or the total formal charges

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
            except:
                print("Getting netcharge error! Using default number 0!")
                netcharge = 0
        else:
            netcharge = netCharge

        tofile = open("antechamber.sh", 'w')
        with open(antechamber) as lines:
            for s in lines:
                if len(s.split()) > 0 and "antechamber" in s.split()[0]:
                    tofile.write("export AMBERBIN=%s/bin \n" % (amberhome))
                    # antechamber -i in.pdb -fi pdb -o prep.AMBER -fo prepi -at AMBER -pf y -s 2 -c bcc -nc 1
                    tofile.write("$AMBERBIN/antechamber -i $1 -fi %s -o prep.$2 -fo prepi -at $2 -pf y -s 2 -c bcc -nc %d \n"
                                 % (infile.split(".")[-1], netcharge))
                else:
                    tofile.write(s)
        tofile.close()

        # run the antechamber script
        try:
            job = sp.Popen("sh ./antechamber.sh %s %s" % (infile, filetype), shell=True)
            job.communicate()
            print("Run antachamber completed. ")
        except SystemError:
            print("Run antachamber failed.")

        return None


    def trimHydrogen(self, reduce, pdbin, pdbout, verbose=False):
        """
        remove hydrogens in pdb file. call amber reduce function.

        Parameters
        ----------
        reduce: str,
            the amber reduce tool, the executable command
        pdbin: str,
            the input pdb file name
        pdbout: str,
            the output pdb file
        verbose: bool, default is False
            whether print detail information

        Returns
        -------

        """

        job = sp.check_output('%s -Trim %s > %s ' % (reduce, pdbin, pdbout), shell=True)

        if verbose:
            print(job)

        return None


class AcpypeGenTop(object):

    def __init__(self, ligand="xxx.lig"):
        self.ligand = ligand

    def run_acpype(self, out_base='ligand'):

        cmd = "acpype -i %s -b %s -d -a gaff2" % (self.ligand, out_base)

        job = sp.Popen(cmd, shell=True)
        job.communicate()

        return None
