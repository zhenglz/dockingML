#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess as sp


class Convert(object):

    """Convert structure files, using obabel (binary file)
    may change another way in future

    Parameters
    ----------
    obabel

    Attributes
    ----------
    obabel

    Methods
    -------

    """

    def __init__(self,  obabel="obabel"):
        self.obabel = obabel
        print("Make sure you have open babel installed and path exported.")

    def convert(self, input, output, verbose=True):
        """convert one format to another, keep hydrogen atoms

        Parameters
        ----------
        input: str,
            input file name
        output: str,
            output file name
        verbose: bool,
            verbose print option

        Returns
        -------

        """

        try:
            if os.path.exists(input):
                status = sp.check_output("%s %s -O %s " %
                                         (self.obabel, input, output),
                                         shell=True)
                if verbose:
                    print(status)
            else:
                print(input + " is not exist!!!\nExit Now!")

        except FileNotFoundError:
            print("Converting file %s failed! "% input)

        return self

    def removeLonePair(self, input, output):
        '''
        input a pdb file, remove the lone pair electron lines in the file
        only when you have a GOLD docking result pose file, you need this function
        :param inpdb: str, a pdb file
        :param outpdb: str, a pdb file
        :return:
        '''

        job = sp.Popen("awk \'$1 ~ /HETATM/ && $3 !~ /XX/ {print $0}\' %s > %s " % (input, output), shell=True)
        job.communicate()

        return 1

    def processHETATM(self, input, output,
                      hetatmsave=['WAT', 'HOH'],
                      dropWater=True,
                      headersave=True,
                      selectedChains=[],
                      ):

        '''

        :param input: str, input pdb file
        :param output: str, output pdb file
        :param hetatmsave: the resnames for saving
        :param dropWater: remove water molecules
        :param cleanedPDB:
        :param headersave: save the pdb header information
        :param selectedChains:
        :return: the pdbfile after removing unnecessary information
        '''
        tofile = open(output, 'wb')
        with open(input) as lines:
            for s in lines:
                if len(s.split()) > 0 and \
                                s.split()[0] in ['ATOM', 'HETATM', 'TER', 'END', 'ENDMDL']:
                    if dropWater:
                        if s[21] in selectedChains and s[17:20].strip() not in ['HOH', 'WAT'] \
                                and s[17:20].strip() in hetatmsave:
                            tofile.write(s)
                    else:
                        if s[21] in selectedChains and s[17:20].strip() in hetatmsave:
                            tofile.write(s)
                else:
                    if headersave:
                        tofile.write(s)

        return 1