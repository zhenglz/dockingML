import numpy as np
import argparse
from argparse import RawTextHelpFormatter
import sys
from dockml import index

class TimeStamp :

    def __init__(self):
        pass

    def selectDataPoints(self, dataf, upbounds, lowbounds, dt=100, usecols=[0, 1]):
        """
        select the time points where values locate in up and low boundaries
        :param dataf: str, input data file name
        :param upbounds: list of float, up boundaries for the two cols
        :param lowbounds:list of float, low boundaries for the two cols
        :param dt: int, delta time step
        :param usecols: list of int, using which two cols values
        :return: ndarray, M*3 shape
        """

        df = np.loadtxt(dataf, comments="#", usecols=usecols)
        timestamp = np.arange(df.shape[0]) * dt

        timestamp = np.array([timestamp]).T

        # add time stamp information
        df = np.concatenate((timestamp, df), axis=1)

        selected = df

        for i in range(len(upbounds)) :
            selected = selected[ selected[:, i+1] < upbounds[i] ]
            selected = selected[ selected[:, i+1] > lowbounds[i]]

        return selected

    def outputIndex(self, output, groupname, indexes):
        """
        output atom index into a gmx index file
        :param output: str, output index file name
        :param groupname: str, group name in index file
        :param indexes: list of interger
        :return:
        """

        tofile = open(output, 'a')
        tofile.write("[ %s ] \n"%(groupname))
        i = 0
        for atom in indexes :
            i += 1
            tofile.write('%6d ' % atom)
            if i % 15 == 0:
                tofile.write('  \n')
        tofile.write(" \n")
        tofile.close()

        return 1

def arguments() :

    d = """
    extract specific frames from a local minimium 
    given a file contains two columns, select data points 
    within up and low boundarys
    
    examples: 
    
    timestampy -h
    
    timestamp.py -dat datafile.dat -up 1.0 1.0 -low -2 -3.0 -out output.ndx -gn groupname
    
    """

    parser = argparse.ArgumentParser(description=d)
    parser.add_argument('-dat', type=str, help="Input data set")
    parser.add_argument('-dt', type=int, default=10,
                        help="Time step used in the dataset file. Default is 10. \n")
    parser.add_argument('-up', type=float, default=[0, 1], nargs="+",
                        help="Upbounds for selection of data points. \n")
    parser.add_argument('-low', type=float, default=[0, 1], nargs="+",
                        help="Lowbounds for selection of data points. \n")
    parser.add_argument('-cols', type=int, default=[0, 1], nargs="+",
                        help="Used cols for datapoint selections. Default is 0 and 1. \n")
    parser.add_argument('-out', type=str, default="index.ndx",
                        help="Output the index of the data points to GMX format file. \n")
    parser.add_argument('-gn', type=str, default="",
                        help="Output the index group name. Default is empty. \n")

    args = parser.parse_args()

    if len(sys.argv) < 2 :
        parser.print_help()
        sys.exit(0)

    return args

def main() :

    args = arguments()

    ts = TimeStamp()

    print("Selecting data points ... ... ")

    dp = ts.selectDataPoints(args.dat, args.up, args.low, dt=args.dt, usecols=args.cols)

    indexes = dp[:, 0]

    if len(args.gn) :
        groupname = args.gn
    else :
        groupname = "+".join([str(x) for x in args.up ]) + "_" + "+".join([str(x) for x in args.low ])

    print("Writing data point time stamp ... ... ")
    ts.outputIndex(args.out, groupname, [int(x) for x in indexes])

    print("Completed ... ... ")
    sys.exit(1)