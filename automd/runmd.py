# -*- coding: utf-8 -*-

import os, sys
from glob import glob
import subprocess as sp


class AutoRunMD :

    def __init__(self, topFile, taskName, grompp, mdrun, verbose=True, qsub=False):
        self.top = topFile
        self.task= taskName
        self.verbose = verbose
        self.grompp = grompp
        self.mdrun = mdrun

    def modifyMDP(self, inputMDPFile, outputMDPFile, parameters):
        """

        :param inputMDPFile:
        :param outputMDPFile:
        :param parameters:
        :return:
        """
        tofile = open(outputMDPFile, 'wb')
        with open(inputMDPFile) as lines :
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
        return(outputMDPFile)

    def energyMinimization(self, emMDPFile, groFile,
                           outTprFile, np=4,
                           qsub=False):
        """

        :param emMDPFile:
        :param groFile:
        :param outTprFile:
        :param np:
        :param qsub:
        :return:
        """
        if self.verbose :
            print("Start energy minimization for task %s."%self.task)

        # generate GMX tpr file
        job = sp.check_output("%s -f %s -c %s -o %s -p %s "%
                              (self.grompp, emMDPFile, groFile, outTprFile, self.top),
                              shell=True )
        if self.verbose :
            print(job)
            print("Generating TPR file for EM completed.")

        # run MD with qsub or on the terminal
        cmd = "mpirun -np %d %s -deffnm %s " % (np, self.mdrun, outTprFile)
        if qsub :
            pscript = PrepScripts('sample_qsub.sh')
            pscript.prepareQsub('em_qsub.sh', cmd)
            job = sp.Popen("qsub em_qsub.sh", shell=True)
            job.communicate()

        else :
            job = sp.Popen(cmd, shell=True)
            job.communicate()

        return(outTprFile+".gro")

    def runMD(self, MDPFile, groFile, outTprFile, mdpOptions,
              sampleScript='sample_qsub.sh',
              qsub=False, np=4,
              index="index.ndx",
              ):

        # prepare mdp file for MD
        current_mdp = MDPFile
        if len(mdpOptions.keys()) :
            self.modifyMDP(inputMDPFile=MDPFile,
                           outputMDPFile="modified_"+MDPFile,
                           parameters=mdpOptions)
            current_mdp = "modified_"+MDPFile

        # prepare tpr file for md
        job = sp.check_output("%s -f %s -c %s -o %s -p %s -n %s"
                              %(self.grompp, current_mdp, groFile, outTprFile, self.top, index),
                              shell=True)
        if self.verbose :
            print(job)

        # run md
        cmd = "mpirun -np %d %s -deffnm %s " % (np, self.mdrun, outTprFile[-4:])
        if qsub :
            pscript = PrepScripts(sampleScript)
            oscript = str(len(glob("./*_"+sampleScript)))+"_"+sampleScript
            pscript.prepareQsub(oscript, cmd)
            job = sp.Popen("qsub %s"%oscript, shell=True)
            job.communicate()
        else :
            job = sp.Popen(cmd+" &", shell=True)
            job.communicate()

        return(outTprFile[-4:]+".gro")


class PrepScripts :

    def __init__(self, qsubScriptSample):
        self.qsubMD = qsubScriptSample

    def prepareQsub(self,  outSH, cmds):
        tofile = open(outSH, 'wb')

        with open(self.qsubMD) as lines:
            for s in lines:
                if "#define_commends" in s:
                    for cmd in cmds :
                        tofile.write(cmd+" \n")
                else:
                    tofile.write(s)
        tofile.close()

    def RewriteAntechamber(self, nc, antechamber="sample_antechamber.sh"):
        tofile = open("antechamber.sh", 'wb')
        with open(antechamber) as lines:
            for s in lines:
                if len(s.split()) > 0 and s.split()[0] == "antechamber":
                    tofile.write(
                        "antechamber -i $1 -fi mol2 -o amberff.prep.$2 -fo prepi -at $2 -pf y -s 2 -c bcc -nc %d \n" % nc)
                else:
                    tofile.write(s)
        return 1
