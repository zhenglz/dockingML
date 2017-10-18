<h1>The Docking-MD-ML Project of PR LBD</h1>
<p>
    Progesterone Receptor (PR) ligand binding domain (LBD) is a 11 helical
    domain with around 250 amino acids. The identification of druggable molecules
    towards this domain is of great use.
    Here we try to combine docking, MD and ML methods to discover novel ligands.
</p>

<br>******************************************************************************<br>
<h2>Part 1 Docking</h2>
    <p>Docking with GOLD</p>
    <p>Post-processing the GOLD results </p>

<h2>Part 2 Molecular Dynamics</h2>
    <p>Perform MD simulations where necessary</p>

<h2>Part 3 Machine Learning</h2>
    <p>
        From the results of md simulation or the results of molecular docking,
        we first extract the important features for ML.
    </p>
    <p>
        Features including, contacts, vdw energy contribution, electrostatic
        contributions. All the features are decomposed into single residues,
        where backbone and sidechain interactions are separated.
    </p>
    <p>
        After which, the features are selected based on the distribution,
        the difference between plus and minus groups. Only limited number of
        features are chosen for ML training.
    </p>
    <h3>
        Step 1. From coordination files to obtain features.
    </h3>
        <p>
            We prepare a input file, called "input.dat". The file contains two columns. First
            column is the filename of a complex file (receptor + ligand) in pdb file format.
            Second column is the ligand code, or residue name of the ligand involving the
            interactions.
        </p>

        <p>We use the BindingFeature Module to extract the raw data features.</p>
        <p>Usage: python BindingFeature.py -h</p>
        <p>using the "input.dat" as an input file</p>
