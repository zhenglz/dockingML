#!/usr/bin/env python

from matplotlib.pyplot import cm
from matplotlib import pyplot as plt
import numpy as np
import sys, os
from argparse import RawTextHelpFormatter
import argparse

import dockml

def plot2dScatter(filenames,
                  xlim=[], ylim=[],
                  xcol=0, ycol=0,
                  colors = [], cmaptype='bwr',
                  title=None, label="FES (kJ/mol)",
                  dpi=2000, savefile=None,
                  fz=12,
                  xshift=0.0, yshift=0.0,
                  xscale=1.0, yscale=1.0, timescale=0.001,
                  xlab="X", ylab="Y",
                  marker='x', alpha=0.8,
                  pmf=False,
                  ) :
    """
    get X data from the first file and Y data from next file in filenames
    draw a time series scatter
    :param filenames:
    :param xlim:
    :param ylim:
    :param xcol:
    :param ycol:
    :param colors:
    :param cmaptype:
    :param title:
    :param label:
    :param dpi:
    :param savefile:
    :param fz:
    :param xshift:
    :param yshift:
    :param xscale:
    :param yscale:
    :param timescale:
    :param xlab:
    :param ylab:
    :param marker:
    :param alpha:
    :param pmf:
    :return:
    """
    X = np.loadtxt(filenames[0], comments=["#", "@"], usecols=[xcol, ycol], dtype=float)[:, 0]
    X = X * xscale + xshift
    Y = np.loadtxt(filenames[1], comments=["#", "@"], usecols=[xcol, ycol], dtype=float)[:, 0]
    Y = Y * yscale + yshift

    if len(colors) == 0 :
        colors = np.loadtxt(filenames[0], comments=["#", "@"],)[:, 0] * timescale

    # draw plot
    plt.scatter(X, Y, c=colors, marker=marker, alpha=alpha, cmap=cmaptype)

    plt.xlabel(xlab, fontsize=fz)
    plt.ylabel(ylab, fontsize=fz)
    if len(xlim) :
        plt.xlim(xlim)
    if len(ylim) :
        plt.ylim(ylim)
    if title :
        plt.title(title)

    # draw color bar if necessary
    cb = plt.colorbar()
    cb.set_label(label, fontsize=12)

    if savefile :
        plt.savefig(savefile, dpi=dpi)

    plt.show()


def plot2dFes(filename, dtype=[], zlim=[],
              xlab="X", ylab="Y", xyzcols=[0, 1, 2],
              title=None, label="FES (kJ/mol)",
              dpi=2000, savefile=None,
              fz=12, level=30,
              pmf=False, MIN=1.0, cmaptype='bwr',
              xshift=0.0, yshift=0.0,
              xlim=[], ylim=[],
              mesh=False,
              ) :
    # load file
    if len(dtype) :
        fes = np.loadtxt(filename, comments=["#","@"],
                         dtype={'names':('Rec', 'Lig', 'Cmap'), 'formats':(dtype[0], dtype[1], dtype[2])},
                         usecols=xyzcols)
        x_size = len(set(fes['Rec']))
        y_size = len(set(fes['Lig']))

        #z = np.reshape(fes['Cmap'], (y_size, x_size))
        z = np.array(fes[:, 2]).reshape((y_size, x_size))
        if mesh :
            x = sorted(list(set(fes['Rec'].astype(np.float))))
            y = sorted(list(set(fes['Lig'].astype(np.float))))
            x.append(x[1] - x[0] + x[-1])
            y.append(y[1] - y[0] + y[-1])

            x = np.asarray(x) + xshift
            y = np.asarray(y) + yshift
        else :
            x = [float(str(item).split("_")[0]) for item in list(fes['Rec'])]
            y = [float(str(item).split("_")[0]) for item in list(fes['Lig'])]
            # 1d arrary to 2d
            x = np.reshape(np.asarray(x) + xshift, (y_size, x_size))
            y = np.reshape(np.asarray(y) + yshift, (y_size, x_size))

    else :
        fes = np.loadtxt(filename, comments="#", usecols=xyzcols)

        # get x and y size
        x_size = len(set(fes[:, 0]))
        y_size = len(set(fes[:, 1]))
        print( "X Y size", x_size, y_size)

        if mesh :
            x = sorted(list(set(fes[:, 0].astype(np.float))))
            y = sorted(list(set(fes[:, 1].astype(np.float))))
            x.append(x[1] - x[0] + x[-1])
            y.append(y[1] - y[0] + y[-1])

            x = np.asarray(x) + xshift
            y = np.asarray(y) + yshift

        else :
            # 1d arrary to 2d
            x = np.reshape(fes[:,0] + xshift, (y_size,x_size))
            y = np.reshape(fes[:,1] + yshift, (y_size,x_size))

        z = np.reshape(fes[:,2], (y_size, x_size))

    if pmf :
        # PMF
        algo = dockml.BasicAlgorithm()
        # vectorize the PMF function, and try to apply it to all element of a list
        PMF = np.vectorize(algo.pmf)
        MAX = np.max(z.ravel())
        z = PMF(z, MIN, kt=2.5, max=MAX)

    # plot data
    cmaptypes = ['bwr', 'seismic', 'Spectral', 'hot', 'cool', 'binary', 'gray']
    cm_cmap = [cm.bwr, cm.seismic, cm.Spectral, cm.hot, cm.cool, cm.binary, cm.gray]
    if mesh :
        plt.pcolormesh(x.astype(np.float), y.astype(np.float), z, cmap=cm_cmap[cmaptypes.index(cmaptype)])
    else :
        plt.contourf(x, y, z, level, cmap=cm_cmap[cmaptypes.index(cmaptype)])

    plt.colorbar(label=label, )
    if len(zlim) :
        plt.clim(vmin=zlim[0], vmax=zlim[1])

    # labels
    plt.xlabel(xlab, fontsize=fz)
    plt.ylabel(ylab, fontsize=fz)

    if len(xlim) :
        plt.xlim(xlim)
    if len(ylim) :
        plt.ylim(ylim)

    if title :
        plt.title(title)

    if savefile :
        plt.savefig(savefile, dpi=dpi)

    plt.show()

    return 1

def plot1dTimeSeries(filename, color, xycol,
                     xlim, ylim,
                     xscale=1.0, yscale=1.0,
                     xstart=0,
                     xlab="", ylab="",
                     title=None, label="1",fz=12,
                     showlegend=True,
                     dpi=2000, savefile=None,
                     show=True,
                     linewidth=0.5,
                     linestyle='', marker='',
                     shiftX=0.0, shiftY=0.0,
                     legend_loc='',
                     legend_box=0,
                     alpha=1.0
                     ):
    # define x y data

    if xycol[0] < 0 :
        xy = np.loadtxt(filename, comments=["#", "@"], usecols=set(xycol))
        y = xy[xstart:, 0] * yscale + shiftY
        x = np.asarray(range(len(list(y)))) * xscale + shiftX

        plt.plot(x, y, color=color, label=label,
                 linestyle=linestyle, marker=marker,
                 linewidth=linewidth, alpha=alpha,
                 )
    else :
        xy = np.loadtxt(filename, comments=["#", "@"], usecols=set(xycol))
        x = xy[xstart:, 0] * xscale + shiftX
        y = xy[xstart:, -1] * yscale + shiftY

        plt.plot(x, y, color=color, label=label,
                 linestyle=linestyle, marker=marker,
                 linewidth=linewidth
                 )

    if title :
        plt.title(title)

    if len(xlab) :
        plt.xlabel(xlab, fontsize=fz)
    if len(ylab) :
        plt.ylabel(ylab, fontsize=fz)

    if len(xlim) :
        plt.xlim(xlim)
    if len(ylim) :
        plt.ylim(ylim)

    #print legend_box
    if len(legend_loc) > 0:
        plt.legend(loc=legend_loc, frameon=legend_box)
    else:
        plt.legend(frameon=legend_box)

    if savefile:
        plt.savefig(savefile, dpi=dpi)

    if show :
        plt.show()

    return 1

def plot1Dhistogram(filename, color,
                    xlim, ylim, xcol=1,
                    xscale=1.0,
                    xstart=0,
                    num_bins=20,
                    relative_prob=True,
                    xlab="", ylab="",
                    title=None,
                    label="1",fz=12,
                    showlegend=True,
                    dpi=2000,
                    savefile=None,
                    show=True,
                    linewidth=0.5,
                    linestyle='', marker='',
                    pmf=False, MIN=1.0,
                    legend_loc='', legend_box=0,
                    alpha=1.0,
                    ):
    data = np.loadtxt(filename, comments=["#","@"], usecols=[xcol])

    X = data[xstart:] * xscale

    hist, bin_edges = np.histogram(X, bins=num_bins)

    if relative_prob :
        prob = hist / float(np.max(hist))
    else :
        prob = hist / float(np.sum(hist))

    #print "Bins Frequency and Bin_Ranges"
    #print prob
    #print bin_edges

    if pmf :
        prob = list(hist)
        for i in range(len(prob)) :
            if prob[i] < MIN :
                prob[i] = MIN / 2.0
        prob = -2.5 * np.log( np.asarray(prob) / float(np.max(prob)))

    plt.plot(bin_edges[1:], prob, color=color, label=label,
             linestyle=linestyle, marker=marker,
             linewidth=linewidth, alpha=alpha,
             )

    if title:
        plt.title(title)

    if len(xlab):
        plt.xlabel(xlab, fontsize=fz)
    if len(ylab):
        plt.ylabel(ylab, fontsize=fz)

    if len(xlim):
        plt.xlim(xlim)
    if len(ylim):
        plt.ylim(ylim)

    if showlegend:
        if len(legend_loc) > 0 :
            plt.legend(loc=legend_loc, frameon=legend_box)
        else :
            plt.legend(frameon=legend_box)

    if savefile:
        plt.savefig(savefile, dpi=dpi)

    if show:
        plt.show()

    return 1

def histBins(files, num_bins=20, xcol=1, xscale=1.0, xstart=0, xshift=0) :
    X, bins = [], []
    for f in files :
        x = np.loadtxt(f, comments=["#", '@'], usecols=[int(xcol)])[xstart: ] * xscale + xshift
        X = X + list(x)

    for i in range(num_bins) :
        bins.append(np.min(X) + float(i) * (np.max(X) - np.min(X)) / float(num_bins))
    return bins

def main() :
    os.chdir(os.getcwd())

    d = '''
        Using plot.py to plot figures in one command line.

        Usage:
        plot.py fes -h
        plot.py xy -h
        plot.py hist -h
        plot.py 2d -h
        
        Notes: For \'FES\', if xcol > ycol, the data will be presented at 
               lower diagonal. If xcol < ycol, the data points will be shown
               at upper diagonal.
               
        '''

    if len(sys.argv) < 2:
        # no enough arguements, exit now
        print( d)
        print( "You chose non of the arguement!\nDo nothing and exit now!\n")
        sys.exit(1)

    try :
        plotType = sys.argv[1]
    except ValueError :
        plotType = 'xy'

    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)

    parser.add_argument('-data', type=str, nargs="+",
                        help="Input data filename. ")
    parser.add_argument('-dtype', type=str, nargs="+", default=[],
                        help="Data type of your input file. Such as S4 and f4 ")
    parser.add_argument('-xshift', default=[0.0, ] * 20, type=float, nargs="+",
                        help="Shift X axis. Default is 0.0. ")
    parser.add_argument('-yshift', default=[0.0, ] * 20, type=float, nargs="+",
                        help="Shift Y axis. Default is 0.0. ")
    parser.add_argument('-xlim', default=[], nargs='+', type=float,
                        help="X limits. Default is auto. \n")
    parser.add_argument('-ylim', default=[], nargs='+', type=float,
                        help="Y limits. Default is auto. \n")
    parser.add_argument('-legend_loc', default='', type=str,
                        help="Set legend location. Options: best, upper right, upper left, lower right, \n"
                             "lower left, center right, center left, lower center, upper center, center \n")
    parser.add_argument('-legend_box', default=0, type=int,
                        help="Set legend box. Default is 0 (off). \n"
                             "Options: 1 (on) and 0 (off). \n")
    parser.add_argument('-alpha', default=1.0, type=float,
                        help="Alpha value, the transparency of lines and markers. Default is 1.0\n")
    if plotType in ['FES','fes','XYZ','xyz', '2d', '2D'] :
        parser.add_argument('-xyzcol', type=int, nargs="+", default=[0, 1, 2], help="Column indeices for xyz data.")
        parser.add_argument('-zlim',type=float, nargs="+", default=[0.0,1.0],
                            help="The range of data in Z axis. \n"
                                 "Default is 0.0 to 1.0 \n")
        parser.add_argument('-label', default="FES", type=str,
                            help="Color bar label in FES plot. Default is FES.")
        parser.add_argument('-level', default=30, type=int,
                            help="Level of the contour figure. ")
        parser.add_argument('-cmaptype', default='bwr',
                            help="Color Map types. Default is bwr. \n"
                                 "Options: bwr, seismic, Spectral, hot, cool, binary, gray.")
        parser.add_argument('-mesh', default=False, type=bool,
                            help="Show Cmap using mesh grid type of plot. Default is False.")
    if plotType in ['timeseries', 'XY', 'xy'] or plotType in ["hist", 'histogram', 'Hist', 'H', '2d', '2D'] :
        parser.add_argument('-xstart', default=0, type=int,
                            help="Start from a specific row, discard the ahead rows.")
        parser.add_argument('-xycol',default=[0,1], type=int, nargs="+",
                            help="X and Y data col. If the first col index is -1,\n"
                                 "no X values defined, use system default x, from \n"
                                 "0 to N. Default is col 0 and col 1 (In python style index).\n")
        parser.add_argument('-colors', default=['red','blue','cyan','orange', 'gray', 'olive', 'yellow', 'black'], type=str, nargs="+",
                            help="Colors for different series. \n"
                                 "Default is red, blue, cyan, orange.\n")
        parser.add_argument('-xyscale', default=[1.0,1.0], type=float, nargs="+",
                            help="Rescale X and Y axis. \n"
                                 "Multiply the x and y axis by a factor.\n"
                                 "Default is [1.0, 1.0].\n")
        parser.add_argument('-labels', default=['1','2','3','4'], type=str, nargs="+",
                            help="Labels and legends for time series data. \n")
        parser.add_argument('-separated', default=0, type=int,
                            help="Show time series data separatedly. \n"
                                 "Default is False. \n")
        parser.add_argument('-linestyle', default=["-", "-.",  '--', ':',]*3, nargs="+", type=str,
                            help="Line styles. ")
        parser.add_argument('-linewidth', default=0.5, type=float,
                            help="Line width. Default is 0.5. ")
        parser.add_argument('-marker', default=[".","x","o","^", "v","s","*", ".", "x", "o", "^"], nargs="+", type=str,
                            help=r"Maker Type. ")
        if plotType in ["hist", 'histogram', 'Hist', 'H'] :
            parser.add_argument('-bins', default=20, type=int,
                                help="Number of bins for histogram.\n "
                                     "Default is 20. \n")
            parser.add_argument('-max2one', default=False, type=bool,
                                help="Maximium Probability to be 1. \n"
                                     "Default is False. \n")

    parser.add_argument('-xlab', type=str, default="XXX",
                        help="X-axis label.")
    parser.add_argument('-ylab', type=str, default="YYY",
                        help="Y-axis label.")

    parser.add_argument('-title', type=str, default=None,
                        help="Title of the figure. \n")
    parser.add_argument('-savefig', default=None, type=str,
                        help="Whether save figure to a file, \n"
                             "if yes, please provide the file name.\n"
                             "Default is None. ")
    parser.add_argument('-dpi', default=2000, type=int,
                        help="Save figure to file using this resolution."
                             "Default is 2000. \n")

    parser.add_argument('-fsize', default=12, type=int,
                        help="Font size of the labels. \n"
                             "Default is 12." )
    parser.add_argument('-pmf', default=False, type=bool,
                        help="From probability to generate PMF. \n"
                             "Default is False. ")
    parser.add_argument('-minX', default=1.0, type=float,
                        help="Min X value if calculating PMF.\n"
                             "Default is 1.0 .\n")
    args, unknown = parser.parse_known_args()
    plotType = sys.argv[1]

    colors = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

    if plotType in ['FES','fes','XYZ','xyz'] :
        for xyzfile in args.data :
            plot2dFes(filename=xyzfile, zlim=args.zlim,dtype=args.dtype,
                      xlab=args.xlab, ylab=args.ylab, xyzcols=args.xyzcol,
                      title=args.title, label= args.label,
                      dpi=args.dpi, savefile=args.savefig,
                      fz=args.fsize, level=args.level,
                      pmf=args.pmf, MIN=args.minX, cmaptype=args.cmaptype,
                      xshift=args.xshift[0], yshift=args.yshift[0],
                      xlim=args.xlim, ylim=args.ylim, mesh=args.mesh,
                      )
    elif plotType in ['timeseries', 'XY', 'xy'] or plotType in ["hist", 'histogram', 'Hist', 'H'] :
        if len(args.data) <= 1 :
            separated = True
        else :
            separated = args.separated

        if separated :
            for i in range(len(args.data)) :
                if len(args.colors) >= len(args.data) :
                    color = args.colors[i]
                else :
                    color = colors[i]

                if len(args.labels) >= len(args.data) :
                    label = args.labels[i]
                else :
                    label = ""
                if plotType in ["hist", 'histogram', 'Hist', 'H'] :
                    plot1Dhistogram(filename=args.data[i],
                                    xscale=args.xyscale[0],
                                    color=color,
                                    xcol=args.xycol[0],
                                    xstart=args.xstart,
                                    xlim=args.xlim,ylim=args.ylim,
                                    xlab=args.xlab, ylab=args.ylab,
                                    label=label,
                                    title=args.title, fz=args.fsize,
                                    showlegend=True, show=True,
                                    savefile=args.savefig, dpi=args.dpi,
                                    linewidth=args.linewidth,
                                    linestyle=args.linestyle[i], marker=args.marker[i],
                                    num_bins=args.bins, relative_prob=args.max2one,
                                    pmf=args.pmf, MIN=args.minX,
                                    legend_loc=args.legend_loc,
                                    legend_box=args.legend_box,
                                    alpha=args.alpha,
                                    )
                elif plotType in ['timeseries', 'XY', 'xy'] :
                    plot1dTimeSeries(filename=args.data[i],
                                     xscale=args.xyscale[0], yscale=args.xyscale[-1],
                                     xstart=args.xstart,
                                     color=color, xycol=args.xycol,
                                     xlim=args.xlim, ylim=args.ylim,
                                     xlab=args.xlab, ylab=args.ylab,
                                     label=label,
                                     title=args.title, fz=args.fsize,
                                     showlegend=True, show=True,
                                     savefile=args.savefig,dpi=args.dpi,
                                     linewidth=args.linewidth,
                                     linestyle=args.linestyle[i], marker=args.marker[i],
                                     shiftX=args.xshift[0], shiftY=args.yshift[0],
                                     legend_loc=args.legend_loc,
                                     legend_box=args.legend_box,
                                     alpha=args.alpha,
                                     )
        else :
            if plotType in ["hist", 'histogram', 'Hist', 'H'] :
                bins_range = histBins(args.data, num_bins=args.bins,
                                      xcol=args.xycol[0],
                                      xscale=args.xyscale[0],
                                      xstart=args.xstart,
                                      xshift=args.xshift[0] )

            for i in range(len(args.data)-1) :
                if len(args.colors) >= len(args.data) :
                    color = args.colors[i]
                else :
                    color = colors[i]

                if len(args.labels) >= len(args.data) :
                    label = args.labels[i]
                else :
                    label = str(i)

                if plotType in ["hist", 'histogram', 'Hist', 'H'] :
                    plot1Dhistogram(filename=args.data[i],
                                    color=color,
                                    xcol=args.xycol[0],
                                    xscale=args.xyscale[0],
                                    xlim=[], ylim=[],
                                    label=label,
                                    xstart=args.xstart,
                                    title=None, fz=args.fsize,
                                    showlegend=False, show=False,
                                    savefile=False, dpi=args.dpi,
                                    linestyle=args.linestyle[i],
                                    marker=args.marker[i],
                                    linewidth=args.linewidth,
                                    num_bins=bins_range,
                                    relative_prob=args.max2one,
                                    pmf=args.pmf, MIN=args.minX,
                                    legend_loc=args.legend_loc,
                                    legend_box=args.legend_box,
                                    alpha=args.alpha,
                                    )
                elif plotType in ['timeseries', 'XY', 'xy'] :
                    plot1dTimeSeries(filename=args.data[i],
                                     color=color, xycol=args.xycol,
                                     xstart=args.xstart,
                                     xscale=args.xyscale[0], yscale=args.xyscale[-1],
                                     xlim=[], ylim=[],
                                     label=label,
                                     title=None, fz=args.fsize,
                                     showlegend=False, show=False,
                                     savefile=False, dpi=args.dpi,
                                     linestyle=args.linestyle[i],
                                     marker=args.marker[i],
                                     linewidth=args.linewidth,
                                     shiftX=args.xshift[i], shiftY=args.yshift[i],
                                     legend_loc=args.legend_loc,
                                     legend_box=args.legend_box,
                                     alpha=args.alpha,
                                     )

            if len(args.labels) >= len(args.data):
                label = args.labels[len(args.data)-1]
            else:
                label = str(len(args.data)-1)
            if len(args.colors) >= len(args.data):
                color = args.colors[len(args.data)-1]
            else:
                color = colors[len(args.data)-1]

            if plotType in ["hist", 'histogram', 'Hist', 'H']:
                plot1Dhistogram(filename=args.data[len(args.data) - 1],
                                color=color,
                                xstart=args.xstart,
                                xcol=args.xycol[0],
                                xscale=args.xyscale[0],
                                xlim=args.xlim,ylim=args.ylim,
                                xlab=args.xlab, ylab=args.ylab,
                                label=label,
                                title=args.title, fz=args.fsize,
                                showlegend=True, show=True,
                                savefile=args.savefig, dpi=args.dpi,
                                linestyle=args.linestyle[len(args.data) - 1],
                                marker=args.marker[len(args.data) - 1],
                                linewidth=args.linewidth,
                                num_bins=bins_range,
                                relative_prob=args.max2one,
                                pmf=args.pmf, MIN=args.minX,
                                legend_loc=args.legend_loc,
                                legend_box=args.legend_box,
                                alpha=args.alpha,
                                )
            elif plotType in ['timeseries', 'XY', 'xy']:
                plot1dTimeSeries(filename=args.data[len(args.data)-1],
                                 color=color, xycol=args.xycol,
                                 xscale=args.xyscale[0], yscale=args.xyscale[-1],
                                 xstart=args.xstart,
                                 xlim=args.xlim, ylim=args.ylim,
                                 xlab=args.xlab, ylab=args.ylab,
                                 label=label,
                                 title=args.title, fz=args.fsize,
                                 showlegend=True, show=True,
                                 savefile=args.savefig, dpi=args.dpi,
                                 linestyle=args.linestyle[len(args.data)-1],
                                 marker=args.marker[len(args.data)-1],
                                 linewidth=args.linewidth,
                                 shiftX=args.xshift[-1], shiftY=args.yshift[-1],
                                 legend_loc=args.legend_loc,
                                 legend_box=args.legend_box,
                                 alpha=args.alpha,
                                 )

    elif plotType in ['2D', '2d', '2ds'] :
        plot2dScatter(args.data, xlim=args.xlim, ylim=args.ylim,
                      xcol=args.xycol[0], ycol=args.xycol[0],
                      colors=[], timescale=0.001, cmaptype=args.cmaptype,
                      xlab=args.xlab, ylab=args.ylab, title=args.title,
                      xscale=args.xyscale[0], yscale=args.xyscale[1],
                      xshift=args.xshift, yshift=args.yshift, marker=args.marker[0],
                      dpi=args.dpi, savefile=args.savefig, fz=args.fsize,
                      )
