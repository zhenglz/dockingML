# -*- coding: utf-8 -*-


class WorkingFlow(object):

    def __init__(self, calc_type):
        self.calc_type = calc_type

        self.permited_types_ = ["community", "pca", "cmpa"]

    def community_network(self):

        TYPE = "community"

        if self.calc_type in TYPE:
            print(workingflow(TYPE))

        return None


def workingflow(calc_type):
    working_flow = {}

    comm = """
    Working flow of the network generation

    1. Correlation analysis / ContactMap analysis
        Mutual information based correlation analysis could be performed using Wordom based
        on C-alpha atoms of the simulation system. The correlation coefficients should be
        within [0,1].

        ContactMap analysis could be performed using cmap class in dockingML module.

    2. Betweeness analysis
        The betweeness between different residues could be calculated from gncommunity software.
        The input for this analysis is the contact matrix or correlation matrix.
        The output of this analysis includes:
            a) the betweenness of the residue matrix
            b) the communities and their memebers

    3. Generate network
        Network based graph could be generated using networkx lib in python

    """
    working_flow["community"] = comm

    return working_flow[calc_type]
