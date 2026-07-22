"""autoMD -- modernised GROMACS + OpenMM MD preparation & running.

Highlights (v2.0):
    - Default force fields updated to AMBER ff19SB/OPC and CHARMM36m/TIP3P-CHARMM
      with OpenFF 2.x for ligands.
    - Uses ParMEd / Interchange for topology conversion; legacy ACPYPE retained.
    - Provides ready-to-use MDP files for modern GROMACS (2022-2026) best
      practice (Verlet scheme, PME, Parrinello-Rahman pressure coupling after
      equilibration, v-rescale thermostat, h-bond constraints, etc.).
    - OpenMM backend scriptable from Python (GPU, mixed precision, PME).
"""
__version__ = "2.0.0"
