#!/usr/bin/env python

import sys


def main():

    if len(sys.argv) < 2:
        print("usage: read_mmpbsa.py summary_full.dat")
        sys.exit(0)

    with open(sys.argv[1]) as lines:
        energies = {}
        for s in lines:
            if "van der Waal energy" in s:
                energies['vdw'] = [float(s.split()[5]), float(s.split()[7])]
            elif "Electrostattic energy" in s:
                energies['ele'] = [float(s.split()[3]), float(s.split()[5])]
            elif "Polar solvation energy"  in s:
                energies['polar'] = [float(s.split()[4]), float(s.split()[6])]
            elif "SASA energy" in s:
                energies['sasa'] = [float(s.split()[3]), float(s.split()[5])]
            elif "Binding energy" in s:
                energies['total'] = [float(s.split()[3]), float(s.split()[5])]

    print("# fn vdw vdw_std ele ele_std polar polar_std sasa sasa_std total total_std")

    string = sys.argv[1] + "  "
    for key in ['vdw', 'ele', 'polar', 'sasa', 'total']:
        string += "%8.3f %8.3f " % (energies[key][0], energies[key][1])

    print(string)

main()
