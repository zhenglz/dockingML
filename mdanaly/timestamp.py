import numpy as np
import argparse
from argparse import RawTextHelpFormatter
import sys


class TimeStamp(object):

    """

    Parameters
    ----------
    fn: str,
        the input dataset file.

    Attributes
    ----------


    """

    def __init__(self, fn, with_index=True):
        self.df = fn
        self.with_index = with_index

        self.selected_ = None

    def selectDataPoints(self, upbounds, lowbounds, dt=100, usecols=[0, 1]):
        """Select the time points where values locate in up and low boundaries

        Parameters
        ----------
        upbounds: list,
            the upper bound of the cols
        lowbounds: list,
            the lower bound of the cols
        dt: int,
            the gap between frames
        usecols: list
            the list of cols for input

        Returns
        -------
        selected: np.ndarray, shape=[N, M]
            the selected dataframe, N is number of frames or samples
            M is the number of dimensions.

        """

        df = np.loadtxt(self.df, comments=["#", "@"], usecols=set(usecols))

        if not self.with_index:
            timestamp = np.arange(df.shape[0]) * dt
            timestamp = np.array([timestamp]).T
            # add time stamp information
            df = np.concatenate((timestamp, df), axis=1)

            self.with_index = True

        selected = df.copy()

        for i in range(len(upbounds)):
            selected = selected[selected[:, i+1] < upbounds[i]]
            selected = selected[selected[:, i+1] > lowbounds[i]]

        self.selected_ = selected

        return selected

    def outputIndex(self, output, groupname):
        """Output atom index into a gromacs index file.

        Parameters
        ----------
        output: str,
            the output index file name
        groupname:
            the output group name

        Returns
        -------

        """

        tofile = open(output, 'a')
        tofile.write("[ %s ] \n" % groupname)
        i = 0

        if self.with_index:
            for atom in self.selected_[:, 0]:
                i += 1
                tofile.write('%6d ' % atom)
                if i % 15 == 0:
                    tofile.write('  \n')
            tofile.write(" \n")
            tofile.close()
        else:
            print("Index information not provided. Exit now!")

        return None


def arguments():

    d = """
    Extract specific frames from a local minimium 
    Given a file contains two columns, select data points 
    within up and low boundarys
    
    Examples: 
    
    timestampy -h
    
    timestamp.py -dat datafile.dat -up 1.0 1.0 -low -2 -3.0 -out output.ndx -gn groupname
    
    """

    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
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

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

    return args


def main():

    args = arguments()

    ts = TimeStamp(args.dat)

    print("Selecting data points ... ... ")

    ts.selectDataPoints(args.up, args.low, dt=args.dt, usecols=args.cols)

    if len(args.gn):
        group_name = args.gn
    else :
        group_name = "+".join([str(x) for x in args.up]) + "_" + "+".join([str(x) for x in args.low])

    print("Writing data point time stamp ... ... ")
    ts.outputIndex(args.out, group_name)

    print("Completed ... ... ")
    sys.exit(1)

