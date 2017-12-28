# -*- coding: utf-8 -*-

from .io import *

import numpy as np
from matplotlib import pyplot as plt
import networkx as nx

import collections
import os, sys
import argparse
from argparse import RawTextHelpFormatter

from dockml import pdbIO

class NetworkPrepare :
    def __init__(self):
        pass

    def genNodeEdges(self, filein, community, output=""):
        """
        input a betweenness file, return community information
        :param filein: str, betweeness data set, matrix format
        :param community: list, community compositions, list of lists
        :return: a betweenness summed community matrix
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

        commu = np.reshape(np.asarray(commbetw), (len(community), len(community)))

        if output :
            np.savetxt(output, commu, delimiter=" ", fmt="%.4f")

        return commu

    def parseNodeEdges(self, filein):
        """
        input a matrix file, return the network information
        :param filein: str, community based matrix file
        :return: a list of sets, [ (edge_i, edge_j, connections), ()]
        """

        dat = np.loadtxt(filein, comments="#")
        nodes = range(dat.shape[0])

        node_edges = []
        for i in nodes:
            for j in nodes:
                if i <= j and dat[i][j] > 0:
                    node_edges.append((i, j, dat[i][j]))

        return node_edges

    def parseCommunities(self, filen):
        """
        input a community information file, return the community data
        :param filen: str, a output community information from gncommunities software
        :return: tuple, ( number_of_communities, { commu_id : residue list })
        """
        modularity = 0.0
        community = collections.defaultdict(list)
        with open(filen) as lines :
            for s in lines :
                if "The optimum number of communities" in s :
                    no_commu = int(s.split()[6])
                    modularity = float(s.split()[-1])

                elif "The residues in community" in s :
                    community[int(s.split()[4])] = [int(x) for x in s.split(":")[-1].split()]

        return (community, modularity)

    def resInDomains(self, domainf, residues):
        '''
        input a list of residues and domain information file
        output the ratio of residues in each domain
        :param domainf:
        :param residues: list, a list of residues from community analysis
        :return:
        '''

        pdb = pdbIO.parsePDB()
        dinfor = pdb.readDomainRes(domainf)

        # eg. { domain_name: [1, 3, 5]}
        domains = collections.defaultdict(list)
        d_count = {}

        for d in dinfor :
            domains[d[0]] = range(d[1], d[2]+1)
            d_count[d[0]] = 0

        for d in d_count.keys() :
            # calculate how many res (in parameter residues) in a domain
            d_count[d] = len(set(domains[d]).intersection(set(residues)))

        # ratio_outof means how much res in the list (residues) in different domains,
        # sum them up, you should get 1
        # eg, this list of residues is 100 residues, only 25 in domain HNH,
        # therefore, the ratio for HNH is 25%.
        ratio_outof    = {}

        # ratio_indomain means, for a specific domain,
        # some ratio of all res in this specific domain, is in the list of residues
        # eg. Domain HNH (have 250 residues), in this community residues list, 25 res in
        # domain HNH, thus the ratio for HNH is 10%
        ratio_indomain = {}

        for d in d_count.keys() :
            ratio_outof[d]    = d_count[d] / float(len(residues))
            ratio_indomain[d] = d_count[d] / float(len(domains[d]))

        return (ratio_indomain, ratio_outof)

class NetworkDraw :
    def __init__(self):
        pass

    def arguemnets(self):
        d = '''
        Descriptions of community network plot
        '''

        parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
        parser.add_argument('-node_edges', type=str, default='node-edges-x.dat',
                            help="File contains node-edges information. \n"
                                 "Default is node-edges.dat. If this file not exist, \n"
                                 "Another node-edge file will be generated from betweenness and community.\n")
        parser.add_argument('-betw', default="betweenness.dat", type=str,
                            help="A matrix file contain betweenness information. \n")
        parser.add_argument('-com', default='community.dat', type=str,
                            help="A result file from gncommunity analysis. \n")
        parser.add_argument('-nsf', type=float, default=100,
                            help="Node size factor. Multiple this number with number of \n"
                                 "Residues in a community to determine the node size.\n")
        parser.add_argument('-lwf', type=float, default=0.001,
                            help="Edge linewidth size factor. Default is 0.001 \n"
                                 "Multiple this number with the \n"
                                 "betweenness in a community to determine the edge size.\n")
        parser.add_argument('-fig', type=str, default='',
                            help="Output the plt figure as a file. Default is figure_1.pdf.\n")
        parser.add_argument('-dpi', type=int, default=2000,
                            help="Output file DPI. Default is 2000. \n")
        parser.add_argument('-label', default=False, type=bool,
                            help="Add labels to nodes. Default is False. \n")
        parser.add_argument('-fsize', default=14, type=int,
                            help="Font size of labels. Default is 16. \n")
        parser.add_argument('-cols', default=['red', 'blue', 'yellow', 'green', 'cyan', 'orange', 'navy',
                                              'pink', 'olive', 'purple', 'firebrick', 'brown'],
                            type=str, nargs="+",
                            help="Colors for the nodes. Default is r b y g c o navy pink, olive \n")
        parser.add_argument('-pos', default='pos.dat', type=str,
                            help="A file contains positions of the nodes. Default is pos.dat. \n"
                                 "If this file is not existed, default postions will be used. \n")
        parser.add_argument('-nres_cutoff', default=6, type=int,
                            help="If in a community, number of residues is less than this number,\n"
                                 "the community will be ignored. \n")

        args, unknown = parser.parse_known_args()

        if len(sys.argv) < 2:
            parser.print_help()
            print("\nYou chose non of the arguement!\nDo nothing and exit now!\n")
            sys.exit(1)

        return args

    def readPos(self, filein):
        """

        :param filein: str, a file contains the communities locations
        :return: list of floats, the positions of the communities
        """
        positions = []
        with open(filein) as lines:
            for s in lines:
                if "#" not in s and len(s.split()) >= 2:
                    positions.append((float(s.split()[0]), float(s.split()[1])))
        return positions

    def drawNetwork(node_edges, nodes, nodes_resnum,
                    nodefactor=300, lwfactor=0.001,
                    showlabel=False, fig=None, DPI=2000,
                    fsize=20, positions=[], colors=[],
                    ):


        '''
        Draw network plot based on weighted node edges
        :param node_edges: a list of sets, []
        :param nodes:
        :param nodes_resnum:
        :param nodefactor:
        :param lwfactor:
        :param showlabel:
        :param fig:
        :param DPI:
        :param fsize:
        :param positions:
        :param colors:
        :return:
        '''

        if len(colors) == 0:
            cols = ['red', 'blue', 'yellow', 'green', 'cyan', 'orange', 'gray', 'pink', 'megenta'] * 2
        elif len(colors) < len(nodes):
            cols = colors * 5
        elif len(colors) >= len(nodes):
            cols = colors
        else:
            cols = ['red', 'blue', 'yellow', 'green', 'cyan', 'orange', 'gray', 'pink', 'megenta'] * 2

        G = nx.Graph()
        G.add_nodes_from(nodes)
        # add edges
        G.add_weighted_edges_from(node_edges)

        node_sizes = np.asarray(nodes_resnum) * nodefactor
        node_colors = cols[:len(nodes)]
        node_labels = {}
        for x in nodes:
            node_labels[x] = "C" + str(x)
        edge_widths = []
        t = [edge_widths.append(x[2]) for x in node_edges]
        edge_widths = [x * lwfactor for x in edge_widths]

        node_pos = {}
        for i in nodes:
            node_pos[i] = positions[i]

        nx.draw_networkx_nodes(G, node_pos, nodelist=nodes, node_size=node_sizes, node_color=node_colors, alpha=1.0)
        nx.draw_networkx_edges(G, pos=node_pos, edge_color='black', width=edge_widths)
        if showlabel:
            nx.draw_networkx_labels(G, pos=node_pos, labels=node_labels, font_size=fsize)

        if fig:
            plt.savefig(fig, dpi=DPI)

        plt.show()

        return 1

def main() :
    #os.chdir(os.getcwd())

    nwd = NetworkDraw()
    nio = CommunityHandler()
    nwp = NetworkPrepare()

    args = nwd.arguemnets()

    comm = nio.readCommunityFile(args.com, nres_cutoff=args.nres_cutoff)
    nodes_resnum = [len(x) for x in comm if len(x) > args.nres_cutoff]

    if os.path.exists(args.node_edges) :
        node_edges = nwp.parseNodeEdges(args.node_edges)
    else :
        node_edges = []
        edges = nwp.genNodeEdges(args.betw, comm)
        nodes = range(edges.shape[0])
        for i in nodes :
            for j in nodes :
                if i < j :
                    node_edges.append((i, j, edges[i][j]))

    nodes = range(len(nodes_resnum))

    colors = args.cols

    if os.path.exists(args.pos) :
        positions = nwd.readPos(args.pos)
    else :
        positions = [
            (0.1, 0.1),
            (0.05, 0.3),
            (0.1, 0.6),
            (0.3, 0.2),
            (0.35, 0.4),
            (0.3, 0.6),
            (0.25, 0.7),
            (0.2, 0.1),
            (0.2, 0.2),
            (0.15, 0.4),
        ] * 2

    nwd.drawNetwork(node_edges, nodes, nodes_resnum, args.nsf, args.lwf, args.label,args.fig, args.dpi, args.fsize, positions, colors,)
