# -*- coding: utf-8 -*-

"""
IO
handle matrix files, community output, network files
"""

class MatrixHandler :
    def __init__(self) :
        pass

class CommunityHandler :
    def __init__(self) :
        pass

    def readCommunityFile(self, filein, nres_cutoff=6):
        """
        read community output file, set a cutoff for least number of community
        """
        comm = []
        with open(filein) as lines :
            for s in lines :
                if "The residues in community" in s :
                    resi = s.split(":")[-1].split()
                    resi = [ (int(x)) for x in resi ]
                    if len(resi) >= nres_cutoff :
                        comm.append(resi)
        return comm
