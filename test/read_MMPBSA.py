#!/usr/bin/env python

import sys

if __name__ == "__main__" :

    filen = sys.argv[1]

    all_data = []

    with open(filen) as lines :

        for s in lines :
            if " = " in s :
                p   = float(s.split("=")[1].split("+/-")[0])
                std = float(s.split("=")[1].split("+/-")[1].split()[0])

                all_data.append(p)
                all_data.append(std)

    print("  ".join([str(x) for x in all_data ]))
