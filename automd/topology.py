#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Modern protein-ligand topology builder
======================================

Supports two workflows:

1. **Preferred**: OpenFF Toolkit + Interchange -> GROMACS (or AMBER/OpenMM)
   inputs. Works for any SMILES/SDF/MOL2 ligand with OpenFF 2.x (e.g.
   ``openff-2.1.0``, the successor to GAFF2) and any standard protein FF
   supported by OpenFF (AMBER ff14SB / ff19SB via ``openmmforcefields``).
   No external AmberTools installation required (uses RDKit + AmberTools
   conda package when available).

2. **Legacy**: AmberTools tleap + ACPYPE for GROMACS conversion
   (preserved for back-compatibility, see ``automd.utils.gentop``).

Example
-------
>>> from automd.topology import ProteinLigandSystem, ForceFieldPreset
>>> sys = ProteinLigandSystem.from_files("receptor.pdb", "ligand.sdf",
...                                      preset="amber19sb-opc")
>>> sys.solvate_box(edge=10.0, ion_strength=0.15)
>>> sys.write_gromacs("topol.top", "conf.gro", posre=True)
"""
from __future__ import annotations

import os
import logging
import tempfile
import shutil
from pathlib import Path
from typing import Optional, List, Tuple

log = logging.getLogger("automd.topology")
log.addHandler(logging.NullHandler())


# ---------------------------------------------------------------------------
class ProteinLigandSystem:
    """High-level topology builder.

    The class tries the OpenFF path first and falls back to AmberTools/ACPYPE
    when OpenFF unavailable.
    """

    def __init__(self,
                 protein_ff: str = "ff14sb",
                 ligand_ff: str = "openff-2.1.0",
                 water_model: str = "opc3",
                 charge_method: str = "am1bcc",
                 ):
        self.protein_ff = protein_ff
        self.ligand_ff = ligand_ff
        self.water = water_model
        self.charge_method = charge_method
        self._interchange = None
        self._pdb = None
        self._ligand_mol = None

    # ------------------------------------------------------------------
    @classmethod
    def from_files(cls, protein_pdb: str, ligand_sdf_or_smi: Optional[str],
                   preset: str = "amber19sb-opc",
                   **kw) -> "ProteinLigandSystem":
        from automd.mdrun import ForceFieldPreset
        p = ForceFieldPreset.lookup(preset)
        ff_map = {
            "amber19sb-opc": ("amber/ff19SB", "openff-2.1.0", "opc3"),
            "amber14sb-tip3p": ("amber/ff14SB", "openff-2.1.0", "tip3p"),
            "charmm36m-tip3p": ("charmm36m", "openff-2.1.0", "charmm36/tip3p-pme-b"),
        }
        prot, lig, water = ff_map[p.name]
        self = cls(protein_ff=prot, ligand_ff=lig, water_model=water, **kw)
        self._protein_pdb = protein_pdb
        self._ligand_file = ligand_sdf_or_smi
        self._build()
        return self

    # ------------------------------------------------------------------
    def _build(self):
        """Dispatch to the available backend."""
        try:
            self._build_openff()
            return
        except ImportError as e:
            log.warning("OpenFF stack unavailable (%s); falling back to legacy "
                        "AmberTools+ACPYPE path.", e)
        self._build_legacy()

    # ------------------------------------------------------------------
    def _build_openff(self):
        """Build topology with the OpenFF toolkit + Interchange."""
        from openff.toolkit import Molecule, ForceField, Topology
        from openff.units import unit
        from openff.interchange import Interchange
        from openmm import app
        import numpy as np

        # protein
        pdb = app.PDBFile(self._protein_pdb)
        protein_top = pdb.topology
        protein_pos = pdb.positions

        # ligand
        lig_mol = None
        if self._ligand_file:
            fname = self._ligand_file
            if fname.lower().endswith(".smi") or fname.lower().endswith(".smiles"):
                with open(fname) as fh:
                    smi = fh.read().split()[0]
                lig_mol = Molecule.from_smiles(smi, allow_undefined_stereo=True)
                lig_mol.generate_conformers(n_conformers=1)
            elif fname.lower().endswith(".sdf"):
                lig_mol = Molecule.from_file(fname, allow_undefined_stereo=True)
            elif fname.lower().endswith(".mol2"):
                lig_mol = Molecule.from_file(fname, allow_undefined_stereo=True)
            else:
                lig_mol = Molecule.from_file(fname, allow_undefined_stereo=True)
            # assign partial charges if missing
            if lig_mol.partial_charges is None:
                lig_mol.assign_partial_charges(self.charge_method,
                                              use_conformers=lig_mol.conformers)

        # load FF; if amber/... requested, use openmmforcefields ForceField for water+protein
        try:
            from openmmforcefields.generators import SystemGenerator, SMIRNOFFTemplateGenerator
            ff_names = []
            if self.protein_ff.startswith("amber/"):
                ff_names += ["amber14/protein.ff14SB.xml" if "ff14" in self.protein_ff
                             else "amber14/protein.ff19SB.xml",
                             "amber14/tip3p_standard.xml" if "tip3p" in self.water
                             else "amber14/opc3.xml"]
            elif "charmm" in self.protein_ff:
                ff_names += ["charmm36.xml", "charmm36/water.xml"]
            smirnoff = SMIRNOFFTemplateGenerator(forcefield=self.ligand_ff,
                                                 molecules=[lig_mol] if lig_mol else [])
            forcefield_omm = app.ForceField(*ff_names)
            forcefield_omm.registerTemplateGenerator(smirnoff.generator)
            modeller = app.Modeller(protein_top, protein_pos)
            self._modeller = modeller
            self._omm_ff = forcefield_omm
            self._lig_mol = lig_mol
            self._backend = "openmm"
            log.info("OpenMM/openff topology builder ready (protein=%s, ligand=%s, water=%s)",
                     self.protein_ff, self.ligand_ff, self.water)
        except Exception as e:
            log.warning("openmmforcefields path failed (%s); trying native Interchange.", e)
            # native openff/interchange path (limited protein support)
            raise

    # ------------------------------------------------------------------
    def solvate_box(self, edge: float = 10.0,
                    ion_strength: float = 0.15,
                    cation: str = "Na+", anion: str = "Cl-",
                    neutralize: bool = True):
        from openmm import app, unit
        from openmm.app import Modeller
        from pdbfixer import PDBFixer  # optional dependency

        modeller = self._modeller
        # Add hydrogens/pH with modeller? We assume receptor is prepped.
        # Solvate
        if "opc3" in self.water:
            water_model = app.OPC3
        elif "tip3p" in self.water:
            water_model = app.TIP3P
        elif "tip4p" in self.water:
            water_model = app.TIP4PEw
        else:
            water_model = app.TIP3P
        modeller.addSolvent(self._omm_ff, model=water_model,
                            padding=edge * unit.angstroms,
                            ionicStrength=ion_strength * unit.molar,
                            positiveIon=cation, negativeIon=anion,
                            neutralize=neutralize)
        self._modeller = modeller
        log.info("System solvated with model=%s padding=%.1f A, ionicStrength=%.3f M",
                 self.water, edge, ion_strength)

    # ------------------------------------------------------------------
    def write_gromacs(self, top_out: str, gro_out: str, posre: bool = True):
        """Write GROMACS top + gro via ParmEd from the OpenMM system."""
        import parmed as pmd
        from openmm import app
        system = self._omm_ff.createSystem(self._modeller.topology,
                                           nonbondedMethod=app.PME,
                                           nonbondedCutoff=1.0 * pmd.unit.nanometer,
                                           constraints=app.HBonds)
        struct = pmd.openmm.load_topology(self._modeller.topology, system,
                                          xyz=self._modeller.positions)
        struct.save(gro_out, overwrite=True)
        struct.save(top_out, overwrite=True)
        if posre:
            # write posre.itp for heavy atoms of protein + ligand
            self._write_posre(struct, "posre.itp")
        log.info("Wrote %s, %s (and posre.itp)", top_out, gro_out)

    # ------------------------------------------------------------------
    @staticmethod
    def _write_posre(struct, outname: str, k: float = 1000.0):
        """Write a GROMACS position-restraints itp for all non-H atoms."""
        with open(outname, "w") as f:
            f.write("[ position_restraints ]\n")
            f.write(";  i funct       fcx        fcy        fcz\n")
            for atom in struct.atoms:
                if atom.atomic_number != 1:
                    f.write(f"{atom.idx+1:5d}    1   {k:.1f}   {k:.1f}   {k:.1f}\n")

    # ------------------------------------------------------------------
    def _build_legacy(self):
        """Fallback to AmberTools tleap + ACPYPE. See automd.utils.gentop."""
        from automd.utils.gentop import GenerateTop, AcpypeGenTop
        log.info("Using legacy AmberTools+ACPYPE topology workflow.")
        self._legacy = (GenerateTop(), AcpypeGenTop(self._ligand_file)
                        if self._ligand_file else None)
        self._backend = "legacy"

    # ------------------------------------------------------------------
    def write_gromacs_legacy(self, out_prefix: str, amberhome: str = os.environ.get("AMBERHOME", ""),
                             box_edge: float = 10.0, ff=None):
        """Run legacy tleap+acpype path."""
        from automd.utils.gentop import GenerateTop, AcpypeGenTop
        gt = GenerateTop()
        ff = ff or ["gaff2", "ff14SB"]
        if self._ligand_file:
            lig_top = AcpypeGenTop(self._ligand_file)
            lig_top.run_acpype("ligand")
            frcmod = "ligand.acpype/ligand_AC.frcmod"
            off = "ligand.acpype/ligand_AC.lib"
            gt.gmxTopBuilder([self._protein_pdb], out_prefix,
                             frcmodFile=frcmod, offFile=off,
                             amberhome=amberhome,
                             solveBox="TIP3PBOX", boxEdge=box_edge,
                             FField=ff)
        else:
            gt.gmxTopBuilder([self._protein_pdb], out_prefix,
                             amberhome=amberhome, solveBox="TIP3PBOX",
                             boxEdge=box_edge, FField=ff)


# ---------------------------------------------------------------------------
def repair_protein_pdb(pdb_in: str, pdb_out: str, pH: float = 7.4,
                       add_missing_residues: bool = False):
    """Repair a PDB file with PDBFixer (add hydrogens, missing atoms).

    Parameters
    ----------
    pdb_in, pdb_out : str
    pH : float
        Protonation state pH.
    add_missing_residues : bool
        If True, model missing loops/residues (PDBFixer can add terminal
        residues; full loop modelling requires MODELLER).
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except ImportError:
        log.warning("pdbfixer not available; skipping repair (copying as-is).")
        shutil.copy(pdb_in, pdb_out)
        return
    fixer = PDBFixer(filename=pdb_in)
    fixer.findMissingResidues()
    if not add_missing_residues:
        # only keep missing residues in the middle; drop termini
        chains = list(fixer.topology.chains())
        keys = list(fixer.missingResidues.keys())
        for key in keys:
            chain = chains[key[0]]
            if key[1] == 0 or key[1] == len(list(chain.residues())):
                del fixer.missingResidues[key]
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH)
    with open(pdb_out, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    log.info("Repaired PDB written to %s", pdb_out)
