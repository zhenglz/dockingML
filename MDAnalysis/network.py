# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import networkx as nx

class NetworkPrepare :
    def __init__(self):
        pass

    def genNodeEdges(self, filein, community, output=""):
        """

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


class NetworkDraw :
    def __init__(self):
        pass

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

