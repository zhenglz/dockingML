The Docking-MD-ML Project of PR LBD

    Progesterone Receptor (PR) ligand binding domain (LBD) is a 11 helical
    domain with around 250 amino acids. The identification of druggable molecules
    towards this domain is of great use.
    Here we try to combine docking, MD and ML methods to discover novel ligands.

******************************************************************************
* Part 1 Docking and docking results parse
    Docking with GOLD
    Post-processing the GOLD results 

* Part 2 Molecular Dynamics
    <p>Perform MD simulations where necessary</p>

* Part 3 Machine Learning
   
        From the results of md simulation or the results of molecular docking,
        we first extract the important features for ML.

        Features including, contacts, vdw energy contribution, electrostatic
        contributions. All the features are decomposed into single residues,
        where backbone and sidechain interactions are separated.

        After which, the features are selected based on the distribution,
        the difference between plus and minus groups. Only limited number of
        features are chosen for ML training.

        Step 1. From coordination files to obtain features.
            We prepare a input file, called "input.dat". The file contains two columns. First
            column is the filename of a complex file (receptor + ligand) in pdb file format.
            Second column is the ligand code, or residue name of the ligand involving the
            interactions.

