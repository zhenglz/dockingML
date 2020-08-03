# Prepare molecule topology for gromacs simulations.

## 1. Generate topology examples
    
    cd examples/10gs/test
    # parameterize protein-ligand complex
    python ../../../generate_top.py -lig 10gs_ligand.mol2 -out complex -ion Na+ Cl- \
           -ff gaff 14SB -inp 10gs_protein.pdb

    # parameterize ligand only 
    python ../../../generate_top.py -lig 10gs_ligand.mol2 -out ligand \
           -ion Na+ Cl- -ff gaff -inp 10gs_protein.pdb

## 2. Run MD simulation 
    cd examples/10gs/test
    python ../../../md_gmx.py -c em_1.gro -p complex_GMX.top -o eq_1 -nsteps 10000 -nt 16 -gmx_suffix "" -gmx_prefix "gmx " -task
    python ../../../md_gmx.py -c eq_1.gro -p complex_GMX.top -o md_1 -nsteps 10000 -nt 16 -gmx_suffix "" -gmx_prefix "gmx " -task md 
