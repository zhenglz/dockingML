#!/usr/bin/env python

import os, sys
import subprocess as sp
import argparse
from argparse import RawDescriptionHelpFormatter


class AutoRunMD :

    def __init__(self):
        self.top = ""

        self.PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))

    def run_suprocess(self, cmd):
        job = sp.Popen(cmd, shell=True)
        job.communicate()

        return self

    def generate_top(self, inpdb, ff="amber99sb-ildn", water='tip3p',
                     outgro="output", top="topol", ignh=True):

        cmd = "gmx pdb2gmx -f %s -o %s -p %s -ff %s -water %s " %\
              (inpdb, outgro, top, ff, water)
        if ignh: cmd += "-ignh "

        self.run_suprocess(cmd)

        return self

    def modify_mdp(self, inmdp, outmdp, parameters):

        tofile = open(outmdp, 'wb')
        with open(inmdp) as lines :
            for s in lines :
                if len(s.split()) > 0 and s[0] != ";" :
                    if s.split()[0] in parameters.keys() \
                            and len(parameters[s.split()[0]]) > 0 :
                        tofile.write("%s    = ")
                        for item in parameters[s.split()[0]] :
                            tofile.write("%s "%str(item))
                        tofile.write(" \n")
                    else:
                        tofile.write(s)
                else :
                    tofile.write(s)

        return self

    def add_box(self, ingro, outgro, ):
        cmd = "gmx editconf -f %s -o %s -c -d 1.2 -bt cubic" % (ingro, outgro)

        self.run_suprocess(cmd)

        return self

    def add_solvent(self, ingro, outgro, intop="topol", spc="spc903.gro"):

        if not os.path.exists(spc):
            spc = os.path.join(self.PROJECT_ROOT, "data/spc903.gro")

        cmd = "gmx solvate -cp %s -cs %s -o %s -p %s " % (ingro, spc, outgro, intop)
        self.run_suprocess(cmd)

        return self

    def add_ions(self, ingro, outgro, emmdp="em_sol.mdp", intop="topol", ion_conc=0.15, sol_index=13):
        if not os.path.exists(emmdp):
            emmdp = os.path.join(self.PROJECT_ROOT, "data/em_sol.mdp")

        cmd1 = "gmx grompp -f %s -c %s -p %s -o %s -maxwarn 100" % (emmdp, ingro, intop, "addion")
        self.run_suprocess(cmd1)

        cmd2 = "echo \"%d 0 0 \" | gmx genion -s %s -p %s -o %s -neutral -conc %.2f " %\
               (sol_index, "addion", intop, outgro, ion_conc)
        self.run_suprocess(cmd2)

        self.top = intop

        return self

    def minimize(self, ingro, outgro, emmdp="em_sol.mdp", intop="topol", nt=4):
        if not os.path.exists(emmdp):
            emmdp = os.path.join(self.PROJECT_ROOT, "data/em_sol.mdp")

        cmd1 = "gmx grompp -f %s -c %s -p %s -o %s -maxwarn 100" % (emmdp, ingro, intop, outgro)
        self.run_suprocess(cmd1)

        cmd2 = "gmx mdrun -deffnm %s -nt %d -v -gpu_id 0" % (outgro, nt)
        self.run_suprocess(cmd2)

        return self

    def md(self, ingro, outgro, nptmdp="npt.mdp", intop="topol",
           nt=4, gpu_ids="1", restraints=False):
        """
        Run MD simulation with this function.

        Examples
        --------
        >>> app = AutoRunMD()
        >>> app.md("em-peptide.gro", "npt_1", "npt_production_run.mdp",
        "topol.top", nt=4, gpu_ids="0,1,2,3")
        >>> print("MD completed.")

        Parameters
        ----------
        ingro : str,
            Input gro file name.
        outgro : str,
            Output file name.
        nptmdp : str, default is npt.mdp (in project root data folder)
            The gromacs configurational file, mdp format.
        intop : str, default is topol.top
            The topology file for MD simulation
        nt : int, default is 4.
            Number of CPU cores for simulations.
        gpu_ids : str, default is "1"
            The GPU ids, comma seperated without space
        restraints : bool, default is False
            Apply restraints to simulation not.

        Returns
        -------

        """

        if not restraints:
            mdp = nptmdp
            if not os.path.exists(nptmdp):
                mdp = os.path.join(self.PROJECT_ROOT+"/data/", nptmdp)
        else:
            if not os.path.exists(nptmdp):
                mdp = "npt_rest.mdp"
                mdp = os.path.join(self.PROJECT_ROOT+"/data/", mdp)

        cmd1 = "gmx grompp -f %s -c %s -p %s -o %s -maxwarn 100" % (mdp, ingro, intop, outgro)
        self.run_suprocess(cmd1)

        cmd2 = "gmx mdrun -deffnm %s -nt %d -v -gpu_id %s" % (outgro, nt, gpu_ids)
        self.run_suprocess(cmd2)

        print("MD Simulation completed. ")

        return self

    def run_app(self, inpdb, outname, mode="solvated",
                ff="amber99sb-ildn", production_run=True,
                preparation=False,
                gpu_ids="1", nt=4, sol_index=13, restraints=False):
        """
        Run the simulation, the final step after data preparation.

        Parameters
        ----------
        inpdb : str,
            The input pdb file name.
        outname : str,
            The naming pattern of output files. For example, if you add water
            to the simulation box, you will have wat_outname.gro file generated.
        mode : str, default = 'solvated'
            The simulation mode.
            gbsa: run implict MD simulation without consider water effect, thus much faster but less accurate.
            solvated: run explict water MD simulation with full-solvated water environmentm, it is slow but accurate.
        production_run : bool, default = True
            Run the product simulation.

        Returns
        -------
            self : an instance of itself

        """

        if restraints and production_run:
            npt_mdp = "npt_rest.mdp"
        else:
            npt_mdp = "npt.mdp"

        if mode == "solvated":
            if preparation:
                self.generate_top(inpdb, outgro=outname, top="topol", ff=ff)
                self.add_box(ingro=outname, outgro="box_"+outname)
                self.add_solvent(ingro="box_"+outname, outgro="wat_"+outname, )
                self.add_ions(ingro="wat_"+outname, outgro="ion_"+outname, sol_index=sol_index)
                self.minimize(ingro="ion_"+outname, outgro="em_"+outname)

            if production_run:
                self.md("em_"+outname, outgro="npt_"+outname, nptmdp=npt_mdp,
                        gpu_ids=gpu_ids, nt=nt, restraints=restraints)
 
        elif mode == "gbsa":
            self.generate_top(inpdb, outgro=outname, top="topol", ff=ff)
            self.add_box(ingro=outname, outgro="box_"+outname)
            self.minimize(ingro="box_"+outname, outgro="em_"+outname, emmdp="em_sol.mdp")

            if production_run:
                mdp = os.path.join(self.PROJECT_ROOT, "data/gbsa.mdp")
                self.md("em_"+outname, outgro="npt_"+outname, nptmdp=mdp, gpu_ids=gpu_ids, nt=nt)

        return self

if __name__ == "__main__":

    d = """
    Run gromacs simulation in one-liner. 
    
    Note: the performance of gromacs simulations is less than openmm simulations.
    
    Examples:
    
        autoRunMD_gmx.py -h
        
        # run simulations with energy minimization followed by production run
        autoRunMD_gmx.py -f input.pdb -o output -gbsa False -product True -ff amber14sb
        
    """
    parser = argparse.ArgumentParser(description=d, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-f", type=str, default="input.pdb",
                        help="Input, str. The input pdb file. Default is input.pdb")
    parser.add_argument("-o", type=str, default="protein",
                        help="Input, str, optional. The output file name pattern. \n"
                             "Default is protein ")
    parser.add_argument("-gbsa", type=lambda x: (str(x).lower() == "true"), default=False,
                        help="Input, bool, optional. Run GBSA implict water simulation. \n"
                             "Default is False. ")
    parser.add_argument("-product_only", type=lambda x: (str(x).lower() == "true"), default=True,
                        help="Input, bool, optional. Run production MD simulations. \n"
                             "Default is True. ")
    parser.add_argument("-prepare_only", type=lambda x: (str(x).lower() == "true"), default=True,
                        help="Input, bool, optional. Run preparation simulations only. \n"
                             "Default is True. ")
    parser.add_argument("-ff", type=str, default="amber99sb-ildn",
                        help="Input, str, optional. Force field type, default is amber99sb-ildn. \n"
                             "Options: amber99sb, amber99sb-ildn, amber12sb, amber14sb and so on. ")
    parser.add_argument("-nt", default=4, type=int, help="Input, int, optional. Number of CPU cores to use.")
    parser.add_argument("-gpuids", default="1", type=str,
                        help="Input, str, optional. The gpu ids to use in the simulation. \n"
                             "Example: 1,0 or 0,1,2,3 or 4")
    parser.add_argument("-SOL_index", default=13, type=int,
                        help="Input, int, optional. The SOL group index, which is used in adding ion to \n"
                             "solvation box step. Default is 13, but could be 15 sometimes.")

    args = parser.parse_args()

    app = AutoRunMD()

    if args.gbsa:
        mode = "gbsa"
    else:
        mode = "solvated"

    app.run_app(inpdb=args.f, outname=args.o, mode=mode, preparation=args.prepare_only,
                production_run=args.product_only,
                ff=args.ff, nt=args.nt, gpu_ids=args.gpuids, sol_index=args.SOL_index)

