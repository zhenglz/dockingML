
from automd.gmxtop import TopModifier, IndexModifier
import os, sys
from automd import gentop, cleanpdb
import subprocess as sp

def cleanMol2(pose, outname) :

    gtop = gentop.GenerateTop()

    gtop.runObabel('obabel', pose, outname+".pdb")

    cpdb = cleanpdb.CleanPDB(outname+".pdb")

    cpdb.removeLonePair(outname+".pdb", "cleaned_"+outname+".pdb")

    gtop.runObabel('obabel', "cleaned_"+outname+".pdb", 'noL_'+outname+'.mol2')

    return 'noL_'+outname+'.mol2'

def ligandTop(ligmol2, nc) :
    gtop = gentop.GenerateTop()

    script = gtop.runAntechamber(ligmol2, netCharge=nc )

    job = sp.Popen('sh ./{} {} {}'.format(script, ligmol2, 'AMBER'), shell=True)
    job.communicate()

def modifyTop() :

    mtop = TopModifier("cplx_GMX.top")
    head = "head_top"
    mtop.addFromFile("topol_1.top", newline_file=head, field="atomtypes", after=True)

    tail = '''
; Include water topology
#include "amber14sb.ff/tip3p.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include topology for ions
#include "amber14sb.ff/ions.itp"
    '''
    mtop = TopModifier("topol_1.top")
    mtop.addFromLines("topol.top", tail, "system", after=False)

    os.remove("topol_1.top")

def motifyIndex() :

    mndx = IndexModifier("index.ndx")

    mndx.appendFields(fields=['Protein', 'LIG'], output="output.ndx", field_name="")

def tt() :

    opt = sys.argv[1]

    if opt in ['top', 'Top', 'TOP'] :
        modifyTop()
    elif opt in ['index', 'ndx', 'Index', 'INDEX'] :
        motifyIndex()

pdb = cleanMol2(sys.argv[1], 'M1')
ligandTop(pdb, 0)
