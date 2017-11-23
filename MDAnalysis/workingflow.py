# -*- coding: utf-8 -*-

def workingflow() :
    wf = """
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
    print(wf)

    return 1
