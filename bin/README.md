## Command line tools for easy MD/Docking analysis

### Set up the environments

#### 1. How to setup the command line tools
<p> Firstly, locate where is the package. Generally, the package
could be found in $PYTHONPATH/lib/python2.7/site-packages/dockingML/bin
</p>
<p>Copy the files in the bin folder to your $HOME/bin/ directory. If 
$HOME/bin is not exisited, create one. Then add the following line in your 
$HOME/.bashrc
</p>
    
    export PATH=$HOME/bin:$PATH

#### 2. If you use git clone https://github.com/zhenglz/dockingML to install the package
    $ cd ~/applications
    $ git clone https://github.com/zhenglz/dockingML
    $ cd dockingML
    $ pip install -e . 
    $ export PATH=$HOME/applications/dockingML/bin:$PATH
    # you could add the above line into your ~/.bashrc 


<p> or add the following line in your ~/.bashrc file: </p>

    export PATH=$PYTHONPATH/lib/python2.7/site-packages/dockingML/bin:$PATH

### Use the tools 

#### 1. Gromacs-style command line
<p>The tools starting with gmx_ would be use to analyze 
trajectory files (xtc, multiple-frame pdb file), if reference
pdb file provided. The tools could also be used to process 
gromacs output xvg files.</p>

##### For example:
<p>gmx_angle.py: reading xtc file to calculate angles/dihedrals of 
selected atoms if a gromacs index file is provided.<p/>
    
    $ gmx_angle.py -h
    
<p>gmx_index.py: reading pdb file to generate gromacs index 
file. </p>
    
    $ gmx_index.py -h

#### 2. Plot dataset files
<p>Plot gromacs output dataset files, or python generated dataset files.
Plot types could be xy, xyz, histogram
</p>

    $ gmx_plot.py -h</p>
    $ gmx_plot.py xy -h</p>

#### 3. Process protein-ligand complex, extract features
<p>Extract binding interation features from the protein ligand 
complex.</p>

    $ genfeatures.py -h

