# dockingML v2.0 — Modern Docking, MD & ML Rescoring Pipeline

**dockingML** is a Python toolkit that integrates **molecular docking**,
**molecular dynamics simulation** and **machine-learning rescoring**
into a reproducible workflow for structure-based drug discovery /
virtual screening.

> **v2.0 highlights (2025):** modern docking engines (Vina 1.2+, GNINA,
> Smina/Vinardo, QVina-W, Vina-CUDA, AutoDock-GPU), updated GROMACS 2022+
> best-practice simulation parameters with **AMBER ff19SB/OPC** and
> **CHARMM36m/TIP3P-CHARMM** presets, vectorised interaction fingerprints
> (OnionNet multi-shell contacts, hydrogen-bond / π-π / salt-bridge /
> halogen-bond PLIFs), classical ML rescoring (XGBoost, LightGBM, RF) with
> consensus, and optional deep-learning (3D-CNN, GNN) rescoring backends.
> OpenFF Toolkit + ParmEd replace the legacy *tlap + acpype* workflow for
> topology generation (legacy path still preserved).

---

## 🖥️ GROMACS Graphical User Interface (GROMACS-GUI)

A Python & PyQt6 Graphical User Interface (GUI) tool for GROMACS molecular dynamics simulation setup, execution, parameter editing, and visualization.

### Quick Start

Launch the GUI on Windows or Linux:

**Windows:**
```cmd
run.bat
```

**Linux / macOS:**
```bash
python3 src/main.py
# or
./run.sh
```

### Features
1. **Ligand Simulator:** Complete protein-ligand complex setup, receptor processing (`pdb2gmx`), topology merge (`complex.gro` & `topol.top`), solvation, neutralization (`genion`), and workflow execution.
2. **Solution Simulator:** Solution system topology creation, water box setup, solvation, and step-by-step EM / NVT / NPT / Production MD.
3. **Visual MDP Parameter Editor:** Graphical parameter editor for `integrator`, `nsteps`, `emtol`, temperature & pressure coupling, and constraints.
4. **Trajectory Analysis & Charting:** Interactive Matplotlib plotting for RMSD, RMSF, Gyrate, and Energy, plus PBC removal (`trjconv`).
5. **External Tool Launcher:** Integration with VMD, PyMOL, and UCSF Chimera.
6. **Modules in Development (WIP):** Membrane Protein Simulator and Polymer Simulator.

---

## Installation

```bash
# 1. Clone the repository
git clone https://github.com/zhenglz/dockingML.git
cd dockingML

# 2. (Recommended) Create a fresh conda/mamba environment
mamba create -n dml python=3.11 -y
mamba activate dml

# 3. Install base dependencies and the package itself
pip install -U pip
pip install -e .

# 4. Optional: install deep-learning stack
pip install -e ".[dl]"

# 5. Optional: install OpenFF Toolkit for modern ligand parameterisation
mamba install -c conda-forge openff-toolkit openff-interchange openmmforcefields pdbfixer -y

# 6. External docking/MD engines (install separately & put on PATH):
#    - AutoDock Vina 1.2+   (https://github.com/ccsb-scripps/AutoDock-Vina)
#    - GNINA 1.0+           (https://github.com/gnina/gnina)
#    - Smina / QVina2-W / Vina-CUDA / AutoDock-GPU (optional)
#    - GROMACS 2022+        (https://manual.gromacs.org/)
#    - OpenBabel 3.x, meeko (pip install meeko) for PDBQT preparation
```

Legacy scripts still work; the original v1 API is preserved.

---

## Project layout

```
dockingML/
├── run.bat                 # Windows GUI launcher
├── run.sh                  # Linux / macOS GUI launcher
├── requirements.txt        # Python GUI dependencies
├── src/                    # GROMACS-GUI Application
│   ├── main.py             # Main GUI application entry point
│   ├── core/               # Backend GROMACS wrappers, MDP editor, analysis
│   └── gui/                # PyQt6 GUI windows, tabs, dialogs, and widgets
├── dockml/                 # core docking & ML code (legacy v1)
│   └── modern/             # v2 pipeline
│       ├── docking.py      # wrappers for Vina / Smina / GNINA / QVina / Vina-CUDA / AD-GPU
│       ├── pipeline.py     # end-to-end dock+featurise+rescore workflow
│       ├── features/       # vectorised contacts, OnionNet, PLIF (HB/pi/ionic/halogen)
│       ├── scoring/        # XGBoost/LightGBM/RF/LR + consensus
│       └── models/         # PyTorch 3D-CNN & GNN rescoring (optional)
├── automd/                 # MD preparation + execution
│   ├── mdrun.py            # modern GROMACS driver + OpenMM backend
│   ├── topology.py         # OpenFF/ParmEd + legacy AmberTools topology builder
│   ├── data/modern/*.mdp   # GROMACS 2022+ MDP files (ff19SB/OPC, CHARMM36m)
│   └── utils/              # legacy helpers
├── mdanaly/                # trajectory analysis (contact map, PCA, PMF, network, …)
├── bin/                    # gmx-style command-line scripts
└── data/                   # data files (atom types, VDW parameters, etc.)
```

---

## Quick start

### 1. Prepare receptor & ligands and dock with GNINA (deep-learning rescoring)

```bash
dml-prepare -r receptor.pdb -l ligand1.sdf ligand2.sdf --repair -o prepared/

dml-dock -r prepared/receptor_fixed.pdbqt \
         --ligand_dir prepared/ \
         --engine gnina \
         --center 15.19 53.90 16.92 --size 22 22 22 \
         --exhaustiveness 8 \
         --outdir docked/
```

The box can also be auto-derived from a co-crystallised ligand:

```python
from dockml.modern.docking import DockingBox
box = DockingBox.from_ligand("cocrystal.sdf", padding=8.0)
```

### 2. Build a protein-ligand system and run MD (GROMACS)

```bash
# Build topology (OpenFF2.1 for ligand, ff19SB+OPC for protein/water), solvate,
# add neutralising ions + 150 mM NaCl, write GROMACS .top/.gro and posre.itp
dml-mdtop -r receptor.pdb -l ligand.sdf \
          --preset amber19sb-opc --padding 10 --repair \
          -o topol

# Run EM → NVT (500 ps, restrained) → NPTeq (1 ns, restrained) → NPT production
dml-md -c topol.gro -p topol.top \
       --length-ns 100 --ncpu 16 --gpu-ids 0 \
       --ligand-resname LIG --prefix md
```

Recommended GROMACS/AMBER presets:

| Preset                | Protein FF   | Water  | Ligand FF  |
|-----------------------|-------------|--------|-----------|
| `amber19sb-opc` (**default/recommended**) | ff19SB | OPC | OpenFF 2.1 (fallback GAFF2) |
| `amber14sb-tip3p`     | ff14SB      | TIP3P | OpenFF 2.1 / GAFF2 |
| `charmm36m-tip3p`     | CHARMM36m   | CHARMM TIP3P | CGenFF/OpenFF |

### 3. Compute interaction features

```bash
# input.dat: two columns <complex_pdb> <lig_resname>
cat > input.dat <<EOF
pose1_complex.pdb  LIG
pose2_complex.pdb  LIG
EOF

dml-features -i input.dat -o features.csv --onion --plif
```

Produces a CSV with:

- per-residue backbone/side-chain **contact counts**, VdW, Coulomb terms
  (legacy-compatible);
- **OnionNet-style multi-shell** element–residue contacts (10 shells from 1–6 Å);
- protein–ligand interaction fingerprints (**hydrophobic, hydrogen-bond,
  π-π stacking, π-cation, salt-bridge, halogen-bond** counts per residue).

### 4. Train an ML rescorer and apply it

```bash
# Train an XGBoost classifier to discriminate actives from decoys
dml-train -X features.csv -y labels.csv -m xgb --task classify -o model.pkl

# Rescore new poses
dml-rescore -f features_new.csv -m model.pkl -o ranked.csv
```

`labels.csv` must contain a binary `label` column (1 = active/good pose,
0 = decoy) aligned row-by-row with `features.csv`.

### 5. End-to-end pipeline

```bash
dml-pipeline -r receptor.pdb -l library/*.sdf \
             --engine gnina \
             --center 15.19 53.90 16.92 --size 22 22 22 \
             --model model.pkl \
             --workdir vs_run/
```

---

## Python API (quick tour)

### Modern docking engines

```python
from dockml.modern.docking import (
    VinaDocking, GninaDocking, SminaDocking, get_engine, DockingBox,
    batch_dock,
)
engine = get_engine("gnina", ncpu=8, exhaustiveness=8, cnn_scoring="rescore")
box = DockingBox(center=(15.2, 53.9, 16.9), size=(22, 22, 22))
engine.dock("rec.pdbqt", "lig.pdbqt", out="out.pdbqt", logfile="out.log", box=box)
print(engine.parse_log("out.log").head())
```

### Vectorised features (fast!)

```python
from dockml.modern.features.contacts import (
    ResidueInteractionFeatures, OnionNetContacts,
)
from dockml.modern.features.interactions import InteractionFingerprints

rf = ResidueInteractionFeatures.from_pdb("complex.pdb", ligand_resname="LIG")
contacts = rf.contact_counts()       # per-residue backbone/side-chain contact counts
vdw      = rf.vdw_energies()         # LJ energies (vectorised)
vec, names = rf.flat_feature_vector()

onion = OnionNetContacts.from_pdb("complex.pdb", ligand_resname="LIG").compute()
plif  = InteractionFingerprints("complex.pdb", "LIG").compute()
```

The legacy `dockml.features.BindingFeature` is still available but the new
vectorised path is **50–200× faster** for typical complexes.

### ML rescoring

```python
from dockml.modern.scoring.classical import ClassicalRescorer, ConsensusRescorer
model = ClassicalRescorer(model="xgb", task="classify")
model.fit(X_train, y_train)
proba = model.predict_proba(X_test)
print("CV:", model.cv_evaluate(X, y))
model.save("model.pkl")
```

Consensus of engine score + ML + CNN affinity:

```python
cons = ConsensusRescorer()
cons.add(model_vina_score, weight=0.3)
cons.add(ml_model,        weight=0.5)
cons.add(gnina_cnn_model, weight=0.2)
final = cons.predict(X)
```

### MD with OpenMM (optional)

```python
from automd.mdrun import OpenMMRunner
from automd.topology import ProteinLigandSystem

sys = ProteinLigandSystem.from_files("receptor.pdb", "ligand.sdf",
                                     preset="amber19sb-opc")
sys.solvate_box(edge=10.0, ion_strength=0.15)
sys.write_gromacs("topol.top", "conf.gro")  # optional GROMACS output

runner = OpenMMRunner(platform="CUDA")
runner.setup_from_amber("out.prmtop", "out.inpcrd")  # or load directly from openmm Modeller
runner.minimize()
runner.equilibrate(nvt_ps=500, npt_ps=1000)
runner.production(length_ns=100, traj="md.dcd", report_ps=20)
```

---

## Trajectory analysis (mdanaly)

Legacy scripts remain available:

| Script               | Purpose                                          |
|----------------------|-------------------------------------------------|
| `gmx_cmap.py`        | residue-residue contact maps                    |
| `gmx_pca.py`         | PCA (xyz / contact-map / dihedral)              |
| `gmx_angle.py`       | angle / dihedral timeseries                     |
| `gmx_network.py`     | community / betweenness network analysis        |
| `gmx_dssp.py`        | secondary-structure time series                 |
| `gmx_pmf2d.py`       | 2D potentials of mean force                     |
| `pmf.py`, `plot.py`, `matrix.py`, … | core analysis utilities        |

---

## Bugs fixed in v2.0 (summary of critical issues in v1.x)

* `BindingFeature.getVdWParams` previously misinterpreted the AMBER
  AtomType.dat columns and never converted units correctly; now uses the
  AMBER `V = ε[(R*/r)^12 − 2(R*/r)^6]` form with proper kcal→kJ conversion.
* `atomicVdWEnergy` was called with an extra `vdwParams` positional
  argument (TypeError at runtime). Fixed.
* Several `atomtype1 = "Cl"` typos inside the `atomtype2` branch have been
  corrected.
* `contactsAtomtype` mis-indexed PDB columns (single-character `line[13]`)
  causing incorrect atom-type counts.
* Legacy MDP files used Berendsen pressure coupling for production; modern
  MDPs use `C-rescale` / `Parrinello-Rahman`, `V-rescale` thermostat,
  Verlet+PME, LINCS h-bond constraints, DispCorr=EnerPres.
* The old `_modify_mdp` over-wrote parameter lines without checking value
  types; the new runner uses an explicit key→value rewriter.

---

## References

If you use dockingML in published work, please cite:

* Zheng L. et al. (manuscript in preparation) — original dockingML.
* For recommended presets:
  - **AMBER ff19SB**: Tian C. et al., *J. Chem. Theory Comput.* **2020**.
  - **OPC water**: Izadi S. et al., *J. Phys. Chem. Lett.* **2014**.
  - **CHARMM36m**: Huang J. et al., *Nat. Methods* **2017**.
  - **OpenFF 2.x**: https://openforcefield.org/
  - **GNINA** (CNN scoring): McNutt A. et al., *J. Cheminform.* **2021**.
  - **Smina/Vinardo**: Koes D.R. et al., *J. Chem. Inf. Model.* **2013**;
    Quiroga R. & Villarreal M., *PLoS ONE* **2016**.
  - **OnionNet**: Liang H. et al., *J. Chem. Inf. Model.* **2020**; OnionNet-2,
    Wang Z. et al. **2022**.

---

## Authors & license

Original code by **Liangzhen Zheng** & **Yuguang Mu** (Nanyang Technological
University, Singapore). Modernised refactor (v2.0) by the dockingML
contributors. Licensed under **GPL-3.0**.
