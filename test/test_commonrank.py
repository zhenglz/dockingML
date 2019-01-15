#!/usr/bin/env python

import sys
import argparse
from dockml import gold

def get_solutions(output, lst) :

    with open(output) as lines :
        all_ligs = [ x.split()[0] for x in lines ]

    results = gold.GoldResults(lst).results

    for lig in all_ligs :
        pose = results[lig][-1].strip("\'")

        pose.replace("data", "scratch")

        gold.GoldResults(lst).copyLigandPoseFile(pose, lig.strip("\'")+".mol2")

if __name__ == "__main__" :

    d = '''
    Get common ranking ligands
    '''
    parser = argparse.ArgumentParser(description=d)

    parser.add_argument('-lst', default=['bestranking.lst'], nargs="+", type=str,
                        help="A list of bestranking files. \n")
    parser.add_argument('-topN', type=int, default=1000,
                        help="Top N molecules in bestrankings")
    parser.add_argument('-out', type=str, default='output.dat',
                        help="Output file name")

    args = parser.parse_args()

    topN = gold.RankResults(args.lst[0]).commonTopLigandsID(args.lst, topNum=args.topN)

    print(topN)

    with open(args.out, 'w') as tofile :
        tofile.write("\n".join(topN))
