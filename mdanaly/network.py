# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import networkx as nx
import collections
import os
import sys
import argparse
from argparse import RawTextHelpFormatter
from dockml import pdbIO


class CommunityHandler(object):

    def __init__(self):
        pass

    def readCommunityFile(self, filein, nres_cutoff=6):
        """
        read community output file, set a cutoff for least number of community

        Parameters
        ----------
        filein: str,
            the commu file generated using gncommunities
        nres_cutoff: int,
            when a node with few residues (<nres_cutoff), ignore the node

        Returns
        -------
        comm: list,
            the communities and their related residues index (starting from zero)

        """
        comm = []
        with open(filein) as lines:
            for s in lines:
                if "The residues in community" in s:
                    resi = s.split(":")[-1].split()
                    resi = [(int(x)) for x in resi]
                    if len(resi) >= nres_cutoff:
                        comm.append(resi)
        return comm


class ParseCommunity(object):
    """
    Parse community information.

    Parameters
    ----------
    commu: str,
        the commu file generated using gncommunities

    Attributes
    ----------
    community: str,
        the commu file generated using gncommunities

    Methods
    -------

    """

    def __init__(self, commu):
        self.community = commu

        if not os.path.exists(self.community):
            print("File not exists: ", self.community)
            sys.exit(0)

    def parseCommunities(self):
        """
        input a community information file, return the community data

        Returns
        -------
        community, modularity: tuple, ( number_of_communities, { commu_id : residue list })
        """

        modularity = 0.0
        community = collections.defaultdict(list)
        with open(self.community) as lines:
            for s in lines:
                if "The optimum number of communities" in s:
                    #no_commu = int(s.split()[6])
                    modularity = float(s.split()[-1])

                elif "The residues in community" in s:
                    community[int(s.split()[4])] = [int(x) for x in s.split(":")[-1].split()]

        return community, modularity

    def genNodeEdges(self, filein, community, output=""):
        """
        input a betweenness file, return community information

        Parameters
        ----------
        filein: str,
            betweeness data set, matrix format
        community: list,
            community compositions, list of lists
        output: str,
            the output file name

        Returns
        -------
        commu:
            a betweenness summed community matrix
        """

        dat = np.loadtxt(filein, comments="#", skiprows=1)
        dat = list(dat)

        commbetw = []

        for comm1 in community:
            for comm2 in community:
                if comm1 == comm2:
                    commbetw += [0.0]
                else :
                    btw = 0.0
                    for i in comm1:
                        for j in comm2:
                            btw += dat[i][j]
                    commbetw += [btw]

        commu = np.reshape(np.asarray(commbetw), (len(community), len(community)))

        if output :
            np.savetxt(output, commu, delimiter=" ", fmt="%.4f")

        return commu


class NetworkPrepare(object):

    def __init__(self):
        pass

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

        ratio_outof = sorted(ratio_outof.items(), key=lambda x: x[1], reverse=True)
        ratio_indomain = sorted(ratio_indomain.items(), key=lambda x:x[1], reverse=True)

        return ratio_indomain, ratio_outof


class NetworkDraw(object):

    def __init__(self):
        pass

    def readPos(self, filein):
        """
        read community locations from a file
        the file should contain at least two columns,
        the last two columns should be the x y coordinates
        :param filein: str, a file contains the communities locations
        :return: list of floats, the positions of the communities
        """
        positions = []
        with open(filein) as lines:
            for s in lines:
                if "#" not in s and len(s.split()) >= 2:
                    positions.append((float(s.split()[-2]), float(s.split()[-1])))
        return positions

    def readColors(self, filein):
        """
        Read colors from a file
        this file contains at least 1 column
        :param filein: str, multiple lines file
        :return: list, list of strings
        """
        colors = []
        with open(filein) as lines:
            for s in [x for x in lines if "#" not in x]:
                colors += s.split()[-1]

        return colors

    def drawNetwork(self, node_edges, nodes, nodes_resnum,
                    nodefactor=300, lwfactor=0.001,
                    showlabel={}, fig=None, DPI=2000,
                    fsize=20, positions=[], colors=[],
                    ):
        """Draw network plot based on weighted node edges

        Parameters
        ----------
        node_edges: a list of sets, []
        nodes: list, int
            the nodes
        nodes_resnum: list, int
            the list of residues in each node
        nodefactor: float,
            the scaling factor for each node
        lwfactor: float,
            the scaling factor for each edge
        showlabel: dict,
            the node-label information
        fig: str,
            the output file name
        DPI: int,
            the output file resolution
        fsize: int,
            the fontsize for the labels
        positions: list

        colors: list

        Returns
        -------

        """

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
            node_labels[x] = str(x)
        edge_widths = []
        t = [edge_widths.append(x[2]) for x in node_edges]
        edge_widths = [x * lwfactor for x in edge_widths]

        node_pos = {}
        for i in nodes:
            node_pos[i] = positions[i]

        nx.draw_networkx_nodes(G, node_pos, nodelist=nodes, node_size=node_sizes,
                               node_color=node_colors, alpha=1.0, edgecolors='black')
        nx.draw_networkx_edges(G, pos=node_pos, edge_color='black', width=edge_widths)
        if len(showlabel.keys()):
            node_labels = showlabel
            nx.draw_networkx_labels(G, pos=node_pos, labels=node_labels, font_size=fsize)

        if fig:
            plt.savefig(fig, dpi=DPI)

        plt.show()

        return None

    def arguemnets(self):
        d = '''
        Community analysis and network plot.
        Calculate the communities from a Cmap of a protein or other biomolecules simulations.
        
        Example:
        Show help information
        network.py -h
        
        Generate community figure
        network.py -betw betweenness.dat -com commu -domf domain_information.dat -nsf 100 -lwf 0.0001 -fig 
        output_figure.pdf -pos postions.dat 
        
        Add labels in the plot:
        network.py -betw betweenness.dat -com commu -domf domain_information.dat -nsf 100 -lwf 0.0001 -fig 
        output_figure.pdf -pos postions.dat -label True
        network.py -betw betweenness.dat -com commu -domf domain_information.dat -nsf 100 -lwf 0.0001 -fig 
        output_figure.pdf -pos postions.dat -label 0 1 2 3 4 5 6 7 8 9 10 11
        '''

        parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
        parser.add_argument('-node_edges', type=str, default='node-edges-x.dat',
                            help="Input, optional. \n"
                                 "File contains node-edges information. \n"
                                 "Default is node-edges.dat. If this file not exist, \n"
                                 "Another node-edge file will be generated from betweenness and community.\n")
        parser.add_argument('-betw', default="betweenness.dat", type=str,
                            help="Input. Default is betweenness.dat \n"
                                 "A matrix file contain betweenness information. \n")
        parser.add_argument('-com', default='commun', type=str,
                            help="Input. Default is commu . \n"
                                 "A result file from gncommunity analysis. \n")
        parser.add_argument('-domf', type=str, default='domains.dat',
                            help="Input. Default is domains.dat. \n"
                                 "Domains and their residue information. \n")
        parser.add_argument('-nsf', type=float, default=100,
                            help="Input, optional. Default is 100. \n"
                                 "Node size factor. Multiple this number with number of \n"
                                 "Residues in a community to determine the node size.\n"
                                 "Default value is 100. \n")
        parser.add_argument('-lwf', type=float, default=0.001,
                            help="Input, optional. Default is 0.001. "
                                 "Edge linewidth size factor. Multiple this number with the \n"
                                 "betweenness in a community to determine the edge size.\n"
                                 "Default value is 0.001 \n")
        parser.add_argument('-fig', type=str, default='',
                            help="Output, default is empty. \n"
                                 "Output the plt figure as a file. Default is figure_1.pdf.\n")
        parser.add_argument('-dpi', type=int, default=2000,
                            help="Input, optional. \n"
                                 "Output file DPI. Default is 2000. \n")
        parser.add_argument('-label', default=[], type=str, nargs="+",
                            help="Input, optional. Default is empty. \n"
                                 "Add labels to nodes. Provide a list of labels for each nodes. \n"
                                 "Examples (4 communities): 0 1 2 3 \n"
                                 "Another exp. (3 communities): HNH RuvC CTD \n"
                                 "If you provide True, the default labeling method would be used.\n")
        parser.add_argument('-fsize', default=14, type=int,
                            help="Input, optional. \n"
                                 "Font size of labels. Default is 16. \n")
        parser.add_argument('-cols', default=['red', 'blue', 'yellow', 'green', 'cyan', 'orange', 'navy',
                                              'pink', 'olive', 'purple', 'firebrick', 'brown'],
                            type=str, nargs="+",
                            help="Input, optional. Default is red blue yellow green cyan orange navy pink, olive \n"
                                 "Colors for the nodes. \n")
        parser.add_argument('-pos', default='pos.dat', type=str,
                            help="Input, optional. Default is pos.dat \n"
                                 "A file contains positions of the nodes. Default is pos.dat. \n"
                                 "If this file is not existed, default postions will be used. \n"
                                 "Example position.dat file content: \n"
                                 "0  0.1 0.1 \n"
                                 "1  0.4 0.2 \n"
                                 "2  0.3 0.15 \n")
        parser.add_argument('-nres_cutoff', default=6, type=int,
                            help="Input, optional. Default value is 6. \n"
                                 "If in a community, number of residues is less than this number,\n"
                                 "the community will be ignored. \n")
        parser.add_argument('-seq_start', default=1, type=int,
                            help="Input, optional. Default is 1. "
                                 "For reporting the residue composition in each node, you need to provide\n"
                                 "a domain.dat file for domain information. Provide the first residue sequence \n"
                                 "number for domain residue composition interpretation. \n")

        args, unknown = parser.parse_known_args()

        if len(sys.argv) < 2:
            parser.print_help()
            print("\nYou chose non of the arguement!\nDo nothing and exit now!\n")
            sys.exit(1)

        return args


def workingflow():

    d = '''
    The working flow of drawing community network
    
    1. construct a cmap
    using distance, or LMI, DCC correlation, generating a residue-residue
    contact/correlation map
    
    2. process the cmap
    remove digonal (means setting digonal as zeroes, since they are tightly connected)
    set values (higher than a probability cutoff, say 0.8) as one, loose connecting values
    (less than the cutoff) as zero
    
    3. calculate community information
    using \'gncommunities\' generating communities, as well as their betweenness
    
    4. get community nodes and edges
    name the communities, get their inter-connection strengths (weighted edges)
    
    5. draw community network
    using networkx generating community network plot 
    
    '''
    print(d)


def main():

    nwd = NetworkDraw()
    nwp = NetworkPrepare()
    args = nwd.arguemnets()

    nio = ParseCommunity(args.com)

    comm, modu = nio.parseCommunities()

    comm_res = [comm[x] for x in comm.keys()
                if len(comm[x]) >= args.nres_cutoff]

    nodes_resnum = [len(x) for x in comm_res]

    # report node residue compositions
    for i in range(len(nodes_resnum)):
        shift_res = args.seq_start
        print("Community (node) %d  Number of residues %d "
              % (i, nodes_resnum[i]))

        info = nwp.resInDomains(args.domf, [x+shift_res for x in comm_res[i]])

        print("Ratio in domain : \n", info[0])
        print("Ratio in community: \n", info[1])

    # generate node-edge information [(node_i, node_j, connectivity)]
    if os.path.exists(args.node_edges):
        node_edges = nwp.parseNodeEdges(args.node_edges)
    else:
        node_edges = []
        edges = nio.genNodeEdges(args.betw, comm_res)
        nodes = range(edges.shape[0])

        for i in nodes:
            for j in nodes:
                if i < j:
                    node_edges.append((i, j, edges[i][j]))

    nodes = range(len(nodes_resnum))

    colors = args.cols
    if os.path.exists(colors[0]):
        colors = nwd.readColors(colors[0])

    print(colors)

    if os.path.exists(args.pos):
        positions = nwd.readPos(args.pos)
    else:
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

    labels = dict()
    if len(args.label):
        if args.label[0] in ["True", "true", "t", "T"]:
            labels = dict(zip(nodes, [str(x) for x in nodes]))
        else:
            try:
                labels = dict(zip(nodes, args.label))
            except IOError:
                print("Your input for labels is not correct. Ignore the labels now. ")

    nwd.drawNetwork(node_edges, nodes, nodes_resnum,
                    args.nsf, args.lwf, labels,
                    args.fig, args.dpi, args.fsize,
                    positions, colors,
                    )
