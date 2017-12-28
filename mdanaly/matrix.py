#!/usr/bin/env python

import numpy as np
import argparse
import os
import sys
from scipy import stats
from argparse import RawTextHelpFormatter
from matplotlib import pyplot as plt

class MatrixHandle :
    def __init__(self):
        pass

    def reshapeMtx(self, dataf, dtype, cols=[0, 1, 2], xyshift=[0, 0]):
        '''
        Load a matrix file, return a ndarray matrix shape object
        :param dataf:
        :param dtype:
        :param cols:
        :param xyshift:
        :return:
        '''

        if len(dtype):
            data = np.loadtxt(dataf, dtype=dtype, comments=['@', '#'], usecols=set(cols))
        else:
            data = np.loadtxt(dataf, comments=['@', '#'], usecols=set(cols))

        xsize = data.shape[0]
        ysize = data.shape[1]

        d = []
        data = list(data)

        for i in range(xsize):
            for j in range(ysize):
                d.append([i, j, data[i][j]])

        for c in [0, 1]:
            d[:, c] = d[:, c] + xyshift[c]

        return np.asarray(d)

    def xyz2matrix(self, data):

        dat = list(data)
        xyz = []

        for i in range(len(data)) :
            for j in range(len(data[0])) :
                xyz.append([i, j, data[i][j]])

        return xyz

    def loadxyz(self, dataf, dtype=[], cols=[0, 1, 2], xyshift=[0, 0]):
        '''
        Load xyz data, return ndarray
        :param dataf:
        :param dtype:
        :param cols:
        :param xyshift:
        :return:
        '''

        if os.path.exists(dataf):
            if len(dtype):

                if "S" in dtype[0]:
                    data = []
                    with open(dataf) as lines :
                        for s in lines :
                            if "#" != s[0] :
                                d = [s.split()[cols[0]].split("_")[0],
                                     s.split()[cols[1]].split("_")[0],
                                     float(s.split()[cols[2]]),
                                     ]
                                data.append(d)

                    data = np.asarray(data)

                if "S" not in dtype[0] and "S" not in dtype[1]:
                    data = np.loadtxt(dataf, dtype=dtype, comments=['@', '#'], usecols=set(cols))
                    for c in [0, 1]:
                        data[:, c] = data[:, c] + xyshift[c]
            else:
                data = np.loadtxt(dataf, comments=['@', '#'], usecols=set(cols))
                data[:, 0] = data[:, 0] + xyshift[0]
                data[:, 1] = data[:, 1] + xyshift[1]
        else:
            print("file %s not exist! Exit now!" % dataf)

            sys.exit(0)

        return data

    def extractDomainData(self, data, xrange, yrange):
        '''
        extract specific x y data based xy range
        :param data: float, 3D xyz array (matrix)
        :param xrange: list, [a, b], both a and b are integers
        :param yrange: list, [a, b]
        :return: 3D xyz array
        '''
        d = []
        data = np.asarray(data)
        num = data.shape[0]
        for i in range(num):
            if xrange[0] <= float(data[i][0]) < xrange[1] \
                    and yrange[0] <= float(data[i][1]) < yrange[1]:
                d.append(list(data[i]))

        return np.asarray(d)

    def neiborhood2zero(self, data, neiborsize=4, xyzshift=[0, 0, 0], zscale=1.0, outtype='xyz'):
        """
        treat diagonal elements as zero
        :param data:
        :param neiborsize:
        :param xyzshift:
        :param zscale:
        :param outtype:
        :return:
        """
        xsize = len(set(list(data[:, 0])))
        ysize = len(set(list(data[:, 1])))
        xysize = xsize * ysize
        newd = data + np.reshape(np.repeat(xyzshift, xysize), [xysize, 3])
        newd[:, 2] = newd[:, 2] * zscale
        for i in range(xysize):
            if abs(data[i, 0] - data[i, 1]) <= neiborsize:
                newd[i, 2] = 0.

        if outtype == 'xyz':
            return newd
        else:
            return np.reshape(newd[:, 2], [xsize, ysize])

    def zRangeSelect(self, data, zrange=[]):
        """
        select the data points whoes z values locate in the zrange
        :param data: a list, or matrix data
        :param zrange: list, the lower and upper boundary of the z values
        :return:
        """

        xyshape = np.asarray(data).shape

        newdata = list(data)

        for x in range(xyshape[0]) :
            for y in range(xyshape[1]) :
                if zrange[0] <= newdata[y][x] < zrange[1] :
                    pass
                else :
                    newdata[y][x] = 0.0

        return newdata

def arguments():
    d = '''
    ########################################################################
    #  Handling matrix files                                               #
    #  Author:  ZHENG Liangzhen & Mu Yuguang                               #
    #  Email:   LZHENG002@e.ntu.edu.sg                                     #
    #  Version: V1.2                                                       #
    #  Date:    27 Dec 2017                                                #
    ########################################################################
    
    Handling matrix files.

    Merge 2 matrix: File1 takes the upper diagonal, while File2 takes the lower diagonal.
    matrixHandle.py -opt merge -dat all_caCmap.dat.xyz ../../../repeat1/4un3_nrqq/cmap_acs9/all_caCmap.dat.xyz -dtype S6 S6 f4

    Extract part of a matrix
    matrixHandle.py -opt extract
    matrixHandle.py -opt extract -dat all_caCmap.dat.xyz -dtype S6 S6 f4 -ds xyz -out HNH_CTD_cmap.xyz -xyrange 765 925 1200 1369
    '''
    parser = argparse.ArgumentParser(description=d, formatter_class=RawTextHelpFormatter)
    parser.add_argument('-opt', default='merge', type=str,
                        help="Operations of matrix files. Default is merge.\n"
                             "Opts: transform, merge, extract, average, pair-t-test, ind-t-test\n"
                             "transform: transform file from xyz to mtx, or vice verse\n"
                             "merge: combine two matrix, 1st matrix takes the upper diagonal.\n"
                             "extract: extract specific x and y data points, using xyrange.\n"
                             "average: averaging matrix z values among the input matrices \n"
                             "pair-t-test and ind-t-test: perform t tests, paired and independednt\n"
                             "domain-aver: averaing domain wise data points\n"
                             "neib0: set diagonal neibor nodes value as 0 "
                        )
    parser.add_argument('-dat', default=["input.dat"], type=str, nargs="+",
                        help="Input data file. Default is input.dat. \n")
    parser.add_argument('-out', default="output.dat", type=str,
                        help="Output file name. Default is output.dat. \n")
    parser.add_argument('-dtype', default=[], type=str, nargs="+",
                        help="Data type. For example, S6 and f4 indicate \n"
                             "string and float data respectively. \n"
                             "Default data type is none, automated data type \n"
                             "will be determined.")
    parser.add_argument('-ds', default='xyz', type=str,
                        help="Data file structure, options: xyz, or matrix. Default is xyz.\n"
                             "xyz data file contains columns (>=3) only. \n"
                             "Matrix data file contains M*N data matrix points. \n")
    parser.add_argument('-xyrange', default=[], type=float, nargs="+",
                        help="Range for x and y. Select only part of the matrix, \n"
                             "thus choose the ranges for X and Y respectively. \n"
                             "Please give 4 numbers, two of them define a range for one axis.")
    parser.add_argument('-drange', default=[], type=str, nargs="+",
                        help="Ranges define domain-domain matrix data. \n"
                             "Data points are average data points. "
                        )
    parser.add_argument('-xyzcol', type=int, nargs="+", default=[0, 1, 2],
                        help="XYZ columns in data file. Default is 0 1 2.")
    parser.add_argument('-xyshift', default=[0, 0], type=float, nargs="+",
                        help="Shift X and Y axis. Add a value to the x or y. \n"
                             "Default is 0 and 0. ")
    parser.add_argument('-zscale', type=float, default=1.0,
                        help="Scale Z column, default is 0 \n")
    parser.add_argument('-dzero', type=bool, default=False,
                        help="Set digonal as zero. Default is False. ")
    parser.add_argument('-neibsize', default=4, type=int,
                        help="if -opt neib0, set neibor size. default is 4. \n")

    args, unknown = parser.parse_known_args()

    if len(sys.argv) < 2:
        # no enough arguements, exit now
        parser.print_help()
        print("\nYou chose non of the arguement!\nDo nothing and exit now!\n")
        sys.exit(1)

    if len(unknown) > 0:
        print("These arguments could not be understood! ")
        print(unknown)

    return args, unknown

def main():

    mtxh = MatrixHandle()

    args, unknown = arguments()

    if args.opt in ["merge", "pair-t-test", "ind-t-test"] :
        if args.ds in ['xyz', 'XYZ', '3d'] :
            data1 = mtxh.loadxyz(args.dat[0], args.dtype, args.xyzcol, xyshift=args.xyshift)
            data2 = mtxh.loadxyz(args.dat[1], args.dtype, args.xyzcol, xyshift=args.xyshift)

        elif args.ds in ['matrix', 'mtx'] :
            data1 = mtxh.reshapeMtx(args.dat[0], args.dtype, args.xyzcol, xyshift=args.xyshift)
            data2 = mtxh.reshapeMtx(args.dat[1], args.dtype, args.xyzcol, xyshift=args.xyshift)

        else :
            sys.exit(0)
        if args.opt == "merge" :
            d = np.zeros(data1.shape)
            for i in range(data1.shape[0]) :
                if data1[i][0] < data1[i][1] :
                    d[i] = data1[i]
                else :
                    d[i] = data2[i]
            np.savetxt(args.out, d, fmt="%.5f")
            print( "Merge matrix file completed!")

        elif args.opt == "pair-t-test" :
            t, p = stats.ttest_rel(data1[:, 2], data2[:, 2])
            print( "T statistics: %6.3f " % t)
            print( "P value     : %8.6f" % p)
        elif args.opt == "ind-t-test" :
            t, p = stats.ttest_ind(data1[:, 2], data2[:, 2])
            print( "T statistics: %6.3f " % t)
            print("P value     : %8.6f" % p)

        fit = np.polyfit(data1[:, 2], data2[:, 2], deg=1)
        x, y = data1[:, 2], data2[:, 2]
        plt.plot(x, fit[0] * x + fit[1], color='red', lw=2.5)
        plt.scatter(x, y, c='b')
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.show()

    elif args.opt == "transform" :
        if args.ds in ['xyz', 'XYZ', '3d'] :
            data = mtxh.loadxyz(args.dat[0], args.dtype, args.xyzcol, xyshift=args.xyshift)
            x_size = len(set(list(data[:, 0])))
            y_size = len(set(list(data[:, 1])))

            odata = np.reshape(data, (x_size, y_size))

            np.savetxt(args.out, odata, fmt="%.5f", delimiter=" ")
        else :
            data = mtxh.reshapeMtx(args.dat[0], args.dtype, args.xyzcol, xyshift=args.xyshift)

            np.savetxt(args.out, data, fmt="%.5f", delimiter=" ")
        print("Transform matrix file type completed!")

    elif args.opt == "extract" :
        if args.ds in ['xyz', 'XYZ', '3d']:
            data = mtxh.loadxyz(args.dat[0], args.dtype, args.xyzcol, xyshift=args.xyshift)
        elif args.ds in ['matrix', 'mtx'] :
            data = mtxh.reshapeMtx(args.dat[0], args.dtype, args.xyzcol, xyshift=args.xyshift)
        else :
            sys.exit(0)
        data = data.astype(np.float)
        d = mtxh.extractDomainData(data, args.xyrange[:2], args.xyrange[2:])
        np.savetxt(args.out, d, fmt="%.5f")
        print("Extract matrix file completed!")

    elif args.opt == "average" :
        aver_data = np.array([])
        for i in range(len(args.dat)) :
            if args.ds in ['xyz', 'XYZ', '3d']:
                data = mtxh.loadxyz(args.dat[0], args.dtype, args.xyzcol, xyshift=args.xyshift)
            elif args.ds in ['matrix', 'mtx'] :
                data = mtxh.reshapeMtx(args.dat[0], args.dtype, args.xyzcol, xyshift=args.xyshift)
            else :
                sys.exit(0)

            print( data.shape)

            if i == 0 :
                aver_data = data
            else :
                aver_data[:, 2] += data[:, 2]

        aver_data[:, 2] = aver_data[:, 2] / float(len(args.dat))

        np.savetxt(args.out, aver_data, fmt="%.5f")
        print( "Average matrix files completed!")

    elif args.opt == "domain-aver" :
        if args.ds in ['xyz', 'XYZ', '3d']:
            data = mtxh.loadxyz(args.dat[0], args.dtype, args.xyzcol, xyshift=args.xyshift)
        elif args.ds in ['matrix', 'mtx']:
            data = mtxh.reshapeMtx(args.dat[0], args.dtype, args.xyzcol, xyshift=args.xyshift)
        else:
            sys.exit(0)
        data = data.astype(np.float)

        drange = []
        if os.path.exists(args.drange[0]) :
            with open(args.drange[0]) as lines :
                for s in lines :
                    if "#" not in s :
                        drange += [ int(s.split()[1]), int(s.split()[2]) ]
        else :
            drange = [ float(x) for x in args.drange ]

        tofile = open(args.out, 'wb')
        for i in range(len(drange)/2) :
            tofile.write("# %d %.2f %.2f \n"%(i, drange[i*2], drange[i*2+1]))

        for i in range(len(drange)/2) :
            for j in range(len(drange) / 2):
                if args.dzero and i == j :
                    tofile.write("%3d %3d  0.0 \n" % (i, j))
                else :
                    d = mtxh.extractDomainData(data, xrange=drange[2*i:2*i+2], yrange=drange[2*j:2*j+2])
                    tofile.write("%3d %3d %12.5f \n"%(i, j, np.mean(d[:, 2])))
        tofile.close()

        print("Domain-wise matrix averaging completed!")

    elif args.opt == 'neib0' :
        if args.ds in ['xyz', 'XYZ', '3d']:
            data = mtxh.loadxyz(args.dat[0], args.dtype, args.xyzcol, xyshift=args.xyshift)
        elif args.ds in ['matrix', 'mtx']:
            data = mtxh.reshapeMtx(args.dat[0], args.dtype, args.xyzcol, xyshift=args.xyshift)
        else:
            sys.exit(0)
        data = data.astype(np.float)

        newd = mtxh.neiborhood2zero(data, neiborsize=args.neibsize, outtype='mtx', zscale=args.zscale)
        np.savetxt(args.out, newd, fmt='%.3f')
