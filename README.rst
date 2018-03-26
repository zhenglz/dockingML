This project includes modules for docking prepare, post-processing,
MD simulation prepare and post-processing.

The original idea of the project is to combine machine learning with
docking to lower down the false positive rate.

Authors: Zheng L.Z. & Mu Y.G.
Contact: lzheng002_AT_e.ntu.edu.sg

* Docking and docking results parse (dockml)
  Docking with GOLD
  Post-processing the GOLD results
  

* Molecular Dynamics (autoMD & mdanaly)
 1 Protein initial structure preparing
    The general problem of a protein pdb file is that the pdb file generated 
    from an X-ray method would suffer from missing atoms, 
    sometimes even long missing loops.

    In order to preform molecular dynamics, we would firstly add missing atoms, 
    sometimes, also need to model the missing loops.

    So we need find a way, or build a pipline to process pdb files (of proteins) 
    to prepare input files for MD and molecular docking.

    Update 2017.3.25 Modeller-9.18 was applied to perform protein processing 
    to model the missing atoms. 
    PDB2PQR was used to add missing residues where necessary.

  2 Ligand topology building
    Using AmberTool and Acpype, charges and bonding, nonbonding parameters
    are calculated using AM1-BCC charge model for large set of ligands
    Amber format topology files are created and converted to GMX format

  3 Protein Ligand Complex simulation
    The complex is then subjected to gromacs for product simulation
  
  4 Trajectory analysis (mdanaly)
    Basic analysis such as contactmap, community analysis, lipid properties,
    time series analysis, PMF (1d, 2d), PCA and so on

* Machine Learning Aid rescoring and virtual screening (dockml)
   
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
    Step 2. Feature clean and selection

    Step 3. Construct machine learning models

    Step 4. Rescoring large screening poses and select high potent ligands


