# -*- coding: utf-8 -*-

import numpy as np

class NetworkPrepare :
    def __init__(self):
        pass

    def genNodeEdges(self, filein, community):
        """

        :param filein: str, betweeness data set, matrix format
        :param community: list, community compositions, list of lists
        :return:
        """

        dat = np.loadtxt(filein, comments="#")
        dat = list(dat)

        commbetw = []

        for comm1 in community :
            for comm2 in community :
                if comm1 == comm2 :
                    commbetw += [0.0]
                else :
                    btw = 0.0
                    for i in comm1 :
                        for j in comm2 :
                            btw += dat[i][j]
                    commbetw += [btw]

        return np.reshape(np.asarray(commbetw), (len(community), len(community)))
