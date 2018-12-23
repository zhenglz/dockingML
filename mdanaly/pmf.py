#!/usr/bin/env python

import numpy as np
import sys
import argparse
from argparse import RawTextHelpFormatter


class PMF(object):
    """
    Perform PMF calculation
    pmf = -RT * ln p_i

    Parameters
    ----------

    Attributes
    ----------

    Methods
    -------

    """

    def __init__(self):
        pass

    def pmf2d(self, data, nbins=20, xcol=0, ycol=1, RT=-2.5):
        """
        give a 2d array, calculate the 2d distribution
        and then calculate its PMF fes

        Parameters
        ----------
        data
        nbins
        xcol
        ycol
        RT: float,
            Gas constant * Absolution temperature

        Returns
        -------
        pmf
        x
        y

        """
        data = np.asarray(data)

        # calculate 2d probability distributions
        hist2d, edges1, edges2 = np.histogram2d(data[:, xcol], data[:, ycol], bins=nbins)

        max_val = float(np.max(hist2d))
        min_val = 1.0 / 2.5

        prob = hist2d / max_val

        # Assign the half_of_the_min to the zero elements
        prob[prob == 0.0] = (min_val / (-1.0 * RT)) / max_val

        pmf = RT * np.log2(prob)

        x = np.repeat(edges1[:-1], pmf.shape[0])
        x = np.reshape(x, (pmf.shape[0], pmf.shape[1]))

        y = np.repeat(edges2[:-1], pmf.shape[0])
        y = np.reshape(y, pmf.shape)

        return pmf.T, x.T, y

    def pmf1d(self, data, nbins=20, RT=-2.5):
        """
        give a 1d array, calculate the 1d distribution
        and then calculate its PMF fes

        Parameters
        ----------
        data: np.ndarray,
            input data set
        nbins: int,
            number of bins for probability distribution
        RT: float,
            Gas constant * Absolution temperature

        Returns
        -------

        """

        # probability distribution
        hist, edges = np.histogram(data, bins=nbins)

        max_val = np.max(hist)
        min_val = np.sort(hist, axis=None)[1]

        prob = hist / max_val
        prob[prob == 0] = (min_val / (-1.0 * RT)) / max_val

        pmf = RT * np.log(prob)

        return pmf, edges


def arguments():
    d = '''
        Transform a 2D data (time series) file into a PMF matrix
        (Only numerical data files are accepted)
        Usage:
        python pmf2d.py -dat your2d.data -numbins numberofbins
    '''

    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)

    parser.add_argument('-dat', type=str, default='time_series.dat',
                        help="Input file name. Two columns or one column. ")
    parser.add_argument('-numbins', type=int, default=40,
                        help="Number of bins for histogram analysis. \n")
    parser.add_argument('-out', type=str, default='',
                        help="Output file name. \n")
    parser.add_argument('-cols', type=int, default=[0], nargs="+",
                        help="Use which cols for analysis. \n "
                             "Default is col 0. ")

    args = parser.parse_args()

    if len(sys.argv) < 2:
        print(d)
        parser.print_help()
        sys.exit(0)

    return args


def main():
    args = arguments()

    dataf = args.dat
    nb = args.numbins
    if len(args.out):
        out = args.out
    else:
        out = "pmf_"+dataf

    df = np.loadtxt(dataf, comments=["#", "@"], usecols= args.cols)

    pmf = PMF()

    if len(args.cols) == 2:
        matrix, edges1, edges2 = pmf.pmf2d(df, nb)
        print("X ticks ")
        print(",".join([str(x) for x in list(edges1)]))
        print("Y ticks ")
        print(",".join([str(x) for x in list(edges2)]))
    else :
        matrix, edges = pmf.pmf1d(df, nb)
        print("X ticks ")
        print(",".join([str(x) for x in list(edges)]))

    np.savetxt(out, matrix, fmt="%.3f")

