import os, sys
import subprocess as sp

if __name__ == "__main__" :

    inp = sys.argv[1]
    info= []
    with open(sys.argv[2]) as lines :
        for s in lines :
            info.append(s.split())

    ligs = []
    with open(inp) as lines :
        for s in lines :
            ligs += [ x.strip("\'") for x in s.split() ]

    results = []
    id = 0
    for lig in ligs :
        id += 1
        for item in info :
            if item[-1] == lig :
                results.append([id, lig, item[0], item[1], item[2], item[3]])
                if os.path.exists(item[3]) :
                    sp.Popen("cp %s ./%d.mol2 "%(item[3], id), shell=True)
                else :
                    print("No such file %s"%item[3])

        print("Completed %d "%id + lig)

    tofile = open("ago_results.dat", 'w')
    for t in results :
        tofile.write(" ".join([str(id)] + [str(x) for x in t ] + [" \n"]))

    print(results)