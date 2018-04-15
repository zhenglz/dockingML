

#!/usr/bin/env python

from matplotlib import pyplot as plt
import numpy as np
import argparse
from argparse import RawTextHelpFormatter
import os
import sys

def loadData(datafile, cols=[0, 1], deli=" ") :
    df = np.loadtxt(datafile, comments=["#", "@"])

    if cols[0] >= 0 :
        X = df[:, cols[0]]
        Y = df[:, cols[1]]
    else :
        Y = df[:, cols[1]]
        X = np.asarray(range(Y.shape[0]))

    return X, Y

if __name__ == "__main__" :
    # PWD
    os.chdir(os.getcwd())

    d = '''
    Generate subplots with Python Matplotlib with one command.

    Examples :
    subplot.py -data num_hbonds_e723_3.xvg num_hbonds_e723_5.xvg num_hbonds_e723_6.xvg num_hbonds_e723_7.xvg num_hbonds_e723_8.xvg -xscale 0.001
    -xlab "Time (ns)" -ylab "Number of Hydrogen Bonds" -labels R899 L901 S902 V903 E904 -linewidth 2 -ylim 0 5.5
    '''

    if len(sys.argv) < 2:
        # no enough arguements, exit now
        print(d)
        print( "You chose non of the arguement!\nDo nothing and exit now!\n")
        sys.exit(1)

    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
    parser.add_argument("-data", type=str, nargs="+", default=[],
                        help="Data files to be ploted. ")
    parser.add_argument("-linestyle", default=["-"]*9, type=str, nargs="+",
                        help="The type of lines. Default is - .")
    parser.add_argument("-linewidth", default=1, type=float, help="Line width. Default is 1.")
    parser.add_argument("-color", default=["red", "blue", "orange", "gray", "cyan", "green"], type=str, nargs="+",
                        help="The color of plots. Default is red and blue .")
    parser.add_argument("-xscale", default=1.0, type=float,
                        help="Scale X axis. The X values will be multiplied with this value.")
    parser.add_argument("-xycols", default=[0, 1] * 9, type=int, nargs="+",
                        help="The column index of the data file. Default is [0, 1]")
    parser.add_argument("-xlab", default="XXX", type=str, help="X axis label. ")
    parser.add_argument("-ylab", default="YYY", type=str, help="Y axis label. ")
    parser.add_argument("-fsize", default=12, type=int, help="Font Size for axis.")
    parser.add_argument("-labels", default=[ str(x) for x in range(9) ], type=str, nargs="+",
                        help="Labels of subplots.")
    parser.add_argument("-legend_loc", default=["upper right"], type=str, nargs="+",
                        help="Set legend location. Default is upper right. Options:"
                             "upper left, upper right, middle.")
    parser.add_argument("-legend_frameon", default=False, type=bool,
                        help="Set boarder box for legends. Default False. ")
    parser.add_argument("-ylim", default=[], type=float, nargs="+",
                        help="Y axis lim. Number of float numbers should be twice "
                             "the number of total input files. Default is auto. ")
    parser.add_argument("-xlim", default=[], type=float, nargs="+",
                        help="X axis lim. ")
    parser.add_argument("-ytick", default=[], type=float, nargs="+",
                        help="Set ytick. Start and end, and step. "
                             "Example: 0 10 1, means from 0 to 9 step as 1.")
    parser.add_argument("-sharex", default=True, type=bool, help="Share X axis. Default is True.")
    parser.add_argument("-sharey", default=True, type=bool, help="Share Y axis. Default is True.")
    parser.add_argument("-savefig", default="", type=str, help="Save Figure. Default is not saving figure.")
    args, unknown = parser.parse_known_args()

    if len(args.legend_loc) == len(args.data) :
        legend_loc = args.legend_loc
    else :
        legend_loc = args.legend_loc[0] * len(args.data)

    if len(args.data) >= 1 :
        num = len(args.data)

        f, axes = plt.subplots(num, sharex=args.sharex, sharey=args.sharey)

        for i in range(num) :
            x1, y1 = loadData(args.data[i], args.xycols[i*2:i*2+2])
            axes[i].plot(x1 * args.xscale, y1,
                         linestyle=args.linestyle[i],
                         color=args.color[i],
                         label=args.labels[i],
                         linewidth=args.linewidth
                         )

            if num % 2 == 1 :
                if i == int(num / 2) :
                    axes[i].set_ylabel(args.ylab, fontsize=args.fsize)
            else :
                pass
                #axes[i].set_ylabel(args.ylab, fontsize=14)
            if i == num -1 :
                axes[i].set_xlabel(args.xlab, fontsize=args.fsize)

            axes[i].legend(frameon=False, loc=legend_loc[i])

            if len(args.ylim) >= 2 :
                axes[i].set_ylim(args.ylim[i*2 : i*2+2])
            if len(args.xlim) >= 2 :
                axes[i].set_ylim(args.xlim[i * 2: i * 2 + 2])

            if len(args.ytick) == 3 :
                ytick = np.arange(args.ytick[0], args.ytick[1], args.ytick[2])
                axes[i].yaxis.set_ticks(ytick)

        f.subplots_adjust(hspace=0)

        if num % 2 == 0 :
            f.text(0.06, 0.5, args.ylab, fontsize=args.fsize,
                   ha='center', va='center', rotation='vertical')

        if len(args.savefig) :
            plt.savefig(args.savefig, dpi=3000)
        plt.show()

    else :
        print("Number of data file is not larger than 2. Exit Now!")
        sys.exit(0)
