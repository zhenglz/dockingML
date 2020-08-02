# Prepare molecule topology for gromacs simulations.

## 1. Generate topology examples
    
    cd examples/10gs/test
    # parameterize protein-ligand complex
    python ../../../generate_top.py -lig 10gs_ligand.mol2 -out complex -ion Na+ Cl- \
           -ff gaff 14SB -inp 10gs_protein.pdb

    # parameterize ligand only 
    python ../../../generate_top.py -lig 10gs_ligand.mol2 -out ligand \
           -ion Na+ Cl- -ff gaff -inp 10gs_protein.pdb

