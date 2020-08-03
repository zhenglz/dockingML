#!/usr/bin/env python

import os, sys, shutil
import subprocess as sp
import argparse
from argparse import RawDescriptionHelpFormatter


class AutoRunMD:
    def __init__(self, top, gro, grompp_exe, mdrun_exe):
        self.top = top
        self.gro = gro
        self.mdp = ""
        self.grompp = grompp_exe
        self.mdrun = mdrun_exe
        self.PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
        self.PROJECT_DATA = os.path.join(self.PROJECT_ROOT, "data")

    def _run_suprocess(self, cmd):
        print("INFO: Now runing the following command: \n%s" % cmd)
        job = sp.Popen(cmd, shell=True)
        job.communicate()

        return self

    def _modify_mdp(self, inmdp, outmdp, parameters):
        tofile = open(outmdp, 'w')
        with open(inmdp) as lines:
            for s in lines:
                if len(s.split()) > 0 and s[0] != ";":
                    if s.split()[0] in parameters.keys() \
                            and len(parameters[s.split()[0]]) > 0:
                        tofile.write(s.split()[0]+"          =  ")
                        for item in parameters[s.split()[0]]:
                            tofile.write("%s " % str(item))
                        tofile.write("; modified by automd \n")
                    else:
                        tofile.write(s)
                else:
                    tofile.write(s)
        return self

    def minimize(self, ingro, out, mdp="em_sol.mdp", nt=4, nsteps=100):
        if not os.path.exists(mdp):
            mdp = os.path.join(self.PROJECT_DATA, "em_sol.mdp")

        # modify mdp file here
        temp_mdp = "xxx_tmp.mdp" #TemporaryFile(prefix="/tmp/")
        params = {
            "nsteps": [nsteps, ],
        }

        new_mdp = "em.mdp"
        shutil.copy(mdp, temp_mdp)
        self._modify_mdp(temp_mdp, new_mdp, params)

        cmd1 = "%s -f %s -c %s -p %s -o %s -maxwarn 100" % (self.grompp, new_mdp, ingro, self.top, out)
        self._run_suprocess(cmd1)
        cmd2 = "%s -deffnm %s -nt %d -v" % (self.mdrun, out, nt)
        self._run_suprocess(cmd2)

        os.remove(temp_mdp)
        return self

    def md(self, ingro, out, mdp="npt.mdp", intop="topol",
           nt=4, gpu_ids="1", restraints=False,
           continue_run=False, nsteps=1000, temperature=300,
           groups=['Water', 'non-Water'], mode="md"):
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
        out : str,
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
            if not os.path.exists(mdp):
                if mode == "md":
                    mdp = os.path.join(self.PROJECT_DATA, "npt.mdp")
                else:
                    mdp = os.path.join(self.PROJECT_DATA, "nvt_eq.mdp")
            else:
                mdp = mdp
        else:
            if not os.path.exists(mdp):
                if mode == "md":
                    mdp = "npt_rest.mdp"
                else:
                    mdp = "nvt_rest_eq.mdp"

                mdp = os.path.join(self.PROJECT_DATA, mdp)
            else:
                mdp = mdp

        # modify mdp file here
        genval = "no" if continue_run else "yes"
        params = {
            "nsteps": [nsteps,],
            "tc-grps": groups,
            "ref_t": [temperature, ] * len(groups),
            "gen_val": [genval, ],
            "energygrps": groups,
        }

        if mode == "md":
            new_mdp = "md.mdp"
        else:
            new_mdp = "eq.mdp"

        #shutil.copy(mdp, temp_mdp)
        self._modify_mdp(mdp, new_mdp, params)

        cmd1 = "%s -f %s -c %s -p %s -o %s -maxwarn 100" % (self.grompp, new_mdp, ingro, intop, out)
        self._run_suprocess(cmd1)

        if len(gpu_ids):
            cmd2 = "%s -deffnm %s -nt %d -v -gpu_id %s" % (self.mdrun, out, nt, gpu_ids)
        else:
            cmd2 = "%s -deffnm %s -nt %d -v" % (self.mdrun, out, nt)

        self._run_suprocess(cmd2)
        print("MD Simulation completed. mode = %s" % mode)
        return self


if __name__ == "__main__":
    d = """
    Run gromacs simulation in one-liner.

    Note: the performance of gromacs simulations is less than openmm simulations.

    Examples:

        md_gmx.py -h

        # run simulations with energy minimization followed by production run
        md_gmx.py -f input.pdb -o output -gbsa False -product True -ff amber14sb

    """
    parser = argparse.ArgumentParser(description=d, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument("-c", type=str, default="input.gro",
                        help="Input, str. The input gro or pdb file. Default is input.pdb")
    parser.add_argument("-p", type=str, default="topology.top",
                        help="Input, str, optional. The input topology file name. \n"
                             "Default is topology.top")
    parser.add_argument("-o", type=str, default="md_1",
                        help="Output file name, optional. \n")
    parser.add_argument("-nsteps", type=int, default=500,
                        help="Input, int. "
                             "The number of steps for simulations.")
    parser.add_argument("-t", type=float, default=300, help="Input, optional. Temperature, default is 300 (K).")
    parser.add_argument("-nt", default=4, type=int, help="Input, int, optional. Number of CPU cores to use.")
    parser.add_argument("-gpuids", default="", type=str,
                        help="Input, str, optional. The gpu ids to use in the simulation. \n"
                             "Example: 1,0 or 0,1,2,3 or 4")
    parser.add_argument("-task", default="md", type=str,
                        help="Input, task type. Options: em, eq, md. \n"
                             "em: energy minimization. eq: nvt equilibrium. md: product MD run.")
    parser.add_argument("-gmx_prefix", type=str, default="gmx ",
                        help="Input, optional. The prefix of mdrun and grompp cmd. ")
    parser.add_argument("-gmx_suffix", type=str, default="_mpi ",
                        help="Input, optional. The suffix of mdrun and grompp cmd. ")
    parser.add_argument("-mdp", default="xxx.mdp", type=str,
                        help="Input, optional. The mdp file for parameter control.")
    parser.add_argument("-c_run", default=0, type=int,
                        help="Input, optional. Continue previous run.")

    args = parser.parse_args()
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    app = AutoRunMD(args.p, args.c,
                    grompp_exe=args.gmx_prefix+"grompp"+args.gmx_suffix,
                    mdrun_exe=args.gmx_prefix+"mdrun"+args.gmx_suffix,)

    if args.task == "em":
        app.minimize(args.c, args.o, mdp=args.mdp, nt=args.nt, nsteps=args.nsteps)
    else:
        app.md(args.c, args.o, mdp=args.mdp, intop=args.p, nt=args.nt,
               gpu_ids=args.gpuids, restraints=False, continue_run=args.c_run,
               nsteps=args.nsteps, temperature=args.t, mode=args.task)

