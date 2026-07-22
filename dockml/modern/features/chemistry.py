"""
Chemical constants & SMARTS definitions for PLIF calculation.
============================================================

- VDW_RADII / ELEMENT_NEGATIVITY
- Residue / element groups used for simple interaction types
- Hydrogen-bond / aromatic / hydrophobic / ionic SMARTS (RDKit-based)

References:
    - OnionNet / OnionNet-2 (Liang et al. 2020; Wang et al. 2022)
    - RDKit-contrib IFP examples
    - MDAnalysis hydrogen-bond criteria
"""
from __future__ import annotations

# Van der Waals radii (Å) from Bondi (1964) / RDKit defaults
VDW_RADII = {
    "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
    "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98,
    "Se": 1.90, "Zn": 1.39, "Fe": 1.24, "Mg": 1.73, "Ca": 1.74,
    "Mn": 1.39, "Cu": 1.40,
}

ELEMENT_NEGATIVITY = {  # Pauling scale
    "H": 2.20, "C": 2.55, "N": 3.04, "O": 3.44, "F": 3.98,
    "P": 2.19, "S": 2.58, "Cl": 3.16, "Br": 2.96, "I": 2.66,
    "Se": 2.55, "Zn": 1.65, "Fe": 1.83, "Mg": 1.31, "Ca": 1.00,
    "Mn": 1.55, "Cu": 1.90,
}

# elements that form hydrophobic contacts (any element is carbon;
# halogens also contribute in halogen-bond sense, here ignored for simplicity)
HYDROPHOBIC_ELEMENTS = {"C"}

# Aromatic residues
AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS"}
# Positively charged (cation) residues
CATION_RES = {"LYS", "ARG", "HIP"}
# Negatively charged (anion) residues
ANION_RES = {"ASP", "GLU"}

# Hydrogen bond donor / acceptor heavy-atom element sets (coarse).
HBOND_DONOR_HEAVY = {"N", "O"}   # N-H, O-H donors
HBOND_ACCEPTOR = {"N", "O", "F"}  # acceptors

# Typical heavy-atom distance cutoffs (Å)
CUTOFF_CONTACT = 4.5
CUTOFF_HBOND = 3.5
CUTOFF_HYDROPHOBIC = 4.5
CUTOFF_PISTACK = 7.0        # aromatic ring centres within this
CUTOFF_SALTBRIDGE = 4.0
CUTOFF_PICATION = 6.0
CUTOFF_HALOGEN = 3.5

# Angle cutoffs (degrees)
ANGLE_HBOND_MIN = 120.0
ANGLE_PISTACK_MIN = 0.0     # angle between ring normals
ANGLE_PISTACK_MAX = 60.0    # face-to-face ~0°, edge-to-face ~90° (T-shape)
