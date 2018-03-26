
from automd.gmxtop import TopModifier, IndexModifier
import os, sys

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

if __name__ == "__main__" :

    opt = sys.argv[1]

    if opt in ['top', 'Top', 'TOP'] :
        modifyTop()
    elif opt in ['index', 'ndx', 'Index', 'INDEX'] :
        motifyIndex()