#!/bin/bash


if [ $# -ne 2 ]; then
    echo "Usage: input RED III.1  output, charge & spin & residuNAME read from Mol2"
    echo "Format: file_in(Mol2); atom_type(gaff or amber)"
    exit
fi

$AMBERBIN/antechamber -i $1 -fi mol2 -o prep.$2 -fo prepi -at $2 -pf y -s 2
$AMBERBIN/parmchk -i prep.$2 -f prepi -o frcmod.$2 -a Y
grep need frcmod.$2

if test -f leap.in; then
   rm leap.in
fi

echo "============================================="
echo "dmfff = loadamberparams frcmod.$2" >> leap.in
echo "loadamberprep prep.$2" >> leap.in

echo "prepared run_parm.sh for tleap"
echo "tleap -s -f ${AMBERHOME}dat/leap/cmd/leaprc.ff99 -f leap.in"
echo "tleap -s -f ${AMBERHOME}dat/leap/cmd/leaprc.gaff -f leap.in"
echo "general AMBER : leaprc.gaff"
echo "AMBER : leaprc.ff90"
