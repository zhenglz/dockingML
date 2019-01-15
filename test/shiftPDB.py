def shiftDNA() :
    import os, sys
    from dockml import pdbIO
    from dockml import algorithms
    import numpy as np
    from automd import shiftpdb

    if len(sys.argv) < 3 :
        print("Number of arguments incorrect! Example: \npython shiftPDB.py input.pdb output.pdb 10.0 \n")
        sys.exit(0)

    dt = float(sys.argv[3])  # angstrom
    out = sys.argv[2]        # output pdb file name
    input = sys.argv[1]      # input pdb file name

    dchains = ['C', 'D']

    if not os.path.exists(input) :
        sys.exit(0)

    with open(input) as lines :
        dna_lines = [ x for x in lines if ("ATOM" in x  and x[21] in dchains) ]
        dna_crds = pdbIO.coordinatesPDB().getAtomCrdFromLines(dna_lines)

    vectors = algorithms.LineFit(np.array(dna_crds)).fit_line()

    vectors = vectors * dt

    tofiles = open(out, 'w')
    shift = shiftpdb.shiftPDB(input)

    with open(input) as lines :
        for atom in lines :
            if "ATOM" in atom and atom[21] in dchains :
                newline = shift.xyzChanger(atom, vectors)
                tofiles.write(newline)
            else :
                tofiles.write(atom)

    tofiles.close()

shiftDNA()
