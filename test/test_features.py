
from dockml import features
import os

print(os.getcwd())

inp = './test/CplxM1.pdb'

ligf = features.BindingFeature()

params = ligf.getVdWParams()

print(params.keys())

receptorXYZ, ligandXYZ, recatomDetailInfor, ligatomDetailInfor = ligf.getAtomInfor(input=inp, ligCode='LIG')

alldistpairs = ligf.atomDistMatrix(receptorXYZ, ligandXYZ)

# for all the atomtypes combinations, what are the contacts counts given a cutoff as 6.0 angstrom?
# this function passes the test and works fine
#atomTypeCounts = ligf.contactsAtomtype(alldistpairs, recatomDetailInfor, ligatomDetailInfor, params, distanceCutoff=6.0, pdbqt=False)

#reslist2, backcount, sidecount = ligf.residueCounts(alldistpairs, recatomDetailInfor, distanceCutoff=4.0)

reslist2, backvan, sidevan = ligf.resVdWContribution(alldistpairs,
                                                     recatomDetailInfor,
                                                     ligatomDetailInfor,
                                                     params,
                                                     maxCutoff=12.0)

#print(atomTypeCounts)
print(sidevan)
print(backvan)
