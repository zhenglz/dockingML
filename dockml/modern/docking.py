#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Modern docking engine wrappers
==============================

Supports:
    - AutoDock Vina 1.2.x  (autodock-vina or vina  on PATH or custom path)
    - AutoDock Vina GPU / Vina-CUDA / QuickVina-W
    - Smina (fork of Vina with Vinardo, AutoDock4, and vina scoring,
             plus flexible residue support)
    - GNINA  1.0+  (CNN scoring / CNN pose scoring, built on Smina)
    - AutoDock-GPU (OpenCL/CUDA AD4 GPU docking)

All wrappers accept input in PDBQT, expose consistent ``dock(receptor, ligand,
box)`` and ``parse_log(logfile)`` entry points, and write results in PDBQT.

Author: dockingML contributors, 2024+
"""
from __future__ import annotations

import os
import glob
import shutil
import logging
import subprocess as sp
from dataclasses import dataclass, field, asdict
from typing import List, Tuple, Dict, Optional, Union
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger("dockml.modern.docking")
log.addHandler(logging.NullHandler())


# ---------------------------------------------------------------------------
# Utility: locate executable
# ---------------------------------------------------------------------------
def _which(exe: str, user_path: Optional[str] = None) -> str:
    """Return the full path to *exe*. Prefer ``user_path`` if given."""
    if user_path:
        if os.path.isfile(user_path) and os.access(user_path, os.X_OK):
            return user_path
        # treat as directory
        p = os.path.join(user_path, exe)
        if os.path.isfile(p) and os.access(p, os.X_OK):
            return p
    found = shutil.which(exe)
    if not found:
        raise FileNotFoundError(
            f"Executable `{exe}` not found on PATH. Please install it or "
            f"pass the full path explicitly."
        )
    return found


# ---------------------------------------------------------------------------
# Docking box dataclass
# ---------------------------------------------------------------------------
@dataclass
class DockingBox:
    """Binding-pocket definition used by all docking engines.

    Parameters
    ----------
    center : (x, y, z)
    size   : (sx, sy, sz) in Angstrom
    """
    center: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    size: Tuple[float, float, float] = (22.0, 22.0, 22.0)

    @classmethod
    def from_ligand(cls, ligand_pdb: str, padding: float = 8.0) -> "DockingBox":
        """Build a box that encloses the ligand coordinates plus padding."""
        try:
            import mdtraj as mt
            t = mt.load(ligand_pdb)
        except Exception as e:
            raise RuntimeError(f"Cannot load {ligand_pdb} to build box: {e}")
        xyz = t.xyz[0] * 10.0  # nm -> A
        center = xyz.mean(axis=0).tolist()
        extent = (xyz.max(axis=0) - xyz.min(axis=0)).tolist()
        size = [extent[i] + 2 * padding for i in range(3)]
        return cls(center=tuple(center), size=tuple(size))


# ---------------------------------------------------------------------------
# Base class
# ---------------------------------------------------------------------------
class BaseDocking:
    """Base docking engine wrapper. Subclasses set ``self.name`` and implement
    ``_build_cmd`` / ``parse_log``."""
    name: str = "base"

    def __init__(self, exe: str, ncpu: int = 8, seed: int = 12345,
                 exhaustiveness: int = 32, num_modes: int = 20,
                 energy_range: float = 5.0, scoring: Optional[str] = None,
                 **extra):
        self.exe = _which(exe)
        self.ncpu = ncpu
        self.seed = seed
        self.exhaustiveness = exhaustiveness
        self.num_modes = num_modes
        self.energy_range = energy_range
        self.scoring = scoring
        self.extra = extra

    # ------------------------------------------------------------------
    def dock(self,
             receptor: str,
             ligand: str,
             out: str = "out.pdbqt",
             logfile: str = "out.log",
             box: DockingBox = DockingBox(),
             flex: Optional[str] = None,
             cwd: Optional[str] = None,
             extra_args: Optional[List[str]] = None) -> int:
        """Run docking. Returns returncode."""
        cmd = self._build_cmd(receptor, ligand, out, logfile, box, flex,
                              extra_args or [])
        log.info("[%s] running: %s", self.name, " ".join(cmd))
        proc = sp.Popen(cmd, cwd=cwd, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = proc.communicate()
        if proc.returncode != 0:
            log.error("[%s] docking failed. stderr:\n%s", self.name,
                      stderr.decode(errors="ignore"))
        return proc.returncode

    # ------------------------------------------------------------------
    def _build_cmd(self, receptor, ligand, out, logfile, box: DockingBox,
                   flex: Optional[str], extra_args: List[str]) -> List[str]:
        raise NotImplementedError

    # ------------------------------------------------------------------
    @staticmethod
    def parse_log(logfile: str) -> pd.DataFrame:
        """Parse engine log → DataFrame with columns ``mode, affinity, rmsd_lb, rmsd_ub``.
        Default works for Vina / Smina / QVina / GNINA."""
        rows = []
        in_table = False
        separator_seen = False
        with open(logfile) as fh:
            for line in fh:
                if "mode |" in line and "affinity" in line:
                    in_table = True
                    separator_seen = False
                    continue
                if in_table and ("-----+" in line or line.strip().startswith("-----")):
                    separator_seen = True
                    continue
                if in_table and separator_seen:
                    s = line.split()
                    if len(s) >= 4:
                        try:
                            mode = int(s[0])
                            aff = float(s[1])
                            lb = float(s[2])
                            ub = float(s[3])
                            rows.append((mode, aff, lb, ub))
                        except ValueError:
                            break
                    elif line.strip() == "":
                        break
        return pd.DataFrame(rows, columns=["mode", "affinity_kcal_mol",
                                           "rmsd_lb", "rmsd_ub"])


# ---------------------------------------------------------------------------
# Vina 1.2.x / AutoDock Vina
# ---------------------------------------------------------------------------
class VinaDocking(BaseDocking):
    """Wrapper for AutoDock Vina (>= 1.2).

    Notes
    -----
    Default parameters reflect 2024 best practice:
        exhaustiveness=32 (higher reproducibility, cost ~4x of default 8)
        num_modes=20
        scoring="vina" is the default.  Vina 1.2 also supports
        ``ad4`` and (in Vina 1.2.5+) ``vinardo`` scoring.
    """
    name = "vina"

    def __init__(self, exe: str = "vina", scoring: str = "vina", **kw):
        super().__init__(exe=exe, scoring=scoring, **kw)

    def _build_cmd(self, receptor, ligand, out, logfile, box, flex, extra):
        cx, cy, cz = box.center
        sx, sy, sz = box.size
        cmd = [self.exe,
               "--receptor", receptor,
               "--ligand", ligand,
               "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}",
               "--center_z", f"{cz:.3f}",
               "--size_x", f"{sx:.3f}", "--size_y", f"{sy:.3f}",
               "--size_z", f"{sz:.3f}",
               "--out", out,
               "--log", logfile,
               "--cpu", str(self.ncpu),
               "--exhaustiveness", str(self.exhaustiveness),
               "--num_modes", str(self.num_modes),
               "--energy_range", str(self.energy_range),
               "--seed", str(self.seed),
               "--scoring", self.scoring or "vina",
               ]
        if flex:
            cmd += ["--flex", flex]
        cmd += list(extra)
        return cmd


# ---------------------------------------------------------------------------
# Smina  (fork with Vinardo scoring and flexible residues)
# ---------------------------------------------------------------------------
class SminaDocking(BaseDocking):
    """Smina/Vinardo docking.  The ``vinardo`` scoring function outperforms
    default Vina on several benchmarks (Quiroga & Villarreal, 2016; McNutt et al.)."""
    name = "smina"

    def __init__(self, exe: str = "smina", scoring: str = "vinardo", **kw):
        super().__init__(exe=exe, scoring=scoring, **kw)
        # Smina uses ``--exhaustiveness`` like Vina but also supports
        # ``--autobox_ligand`` (when box is not explicitly given).
        self.autobox_add = kw.get("autobox_add", 4)

    def _build_cmd(self, receptor, ligand, out, logfile, box, flex, extra):
        cx, cy, cz = box.center
        sx, sy, sz = box.size
        cmd = [self.exe,
               "-r", receptor,
               "-l", ligand,
               "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}",
               "--center_z", f"{cz:.3f}",
               "--size_x", f"{sx:.3f}", "--size_y", f"{sy:.3f}",
               "--size_z", f"{sz:.3f}",
               "-o", out,
               "--log", logfile,
               "--cpu", str(self.ncpu),
               "--exhaustiveness", str(self.exhaustiveness),
               "--num_modes", str(self.num_modes),
               "--energy_range", str(self.energy_range),
               "--seed", str(self.seed),
               "--scoring", self.scoring,
               "--autobox_add", str(self.autobox_add),
               ]
        if flex:
            cmd += ["--flex", flex]
        cmd += list(extra)
        return cmd


# ---------------------------------------------------------------------------
# GNINA  (deep-learning augmented docking/rescoring)
# ---------------------------------------------------------------------------
class GninaDocking(BaseDocking):
    """GNINA 1.0+ docking with CNN pose/affinity rescoring.

    Key parameters (defaults from GNINA paper McNutt et al., J Cheminform 2021):
        exhaustiveness=8    (GNINA default; CNN rescoring compensates for lower sampling)
        cnn_scoring="rescore"   poses sampled with Vina SF, then CNN rescores
        cnn="default2018"       DenseNet ensemble
        num_modes=9
    """
    name = "gnina"

    def __init__(self,
                 exe: str = "gnina",
                 cnn_scoring: str = "rescore",
                 cnn_model: str = "default2018",
                 exhaustiveness: int = 8,
                 num_modes: int = 9,
                 cnn_rotation: int = 0,
                 min_rmsd_filter: float = 1.0,
                 num_mc_saved: int = 50,
                 **kw):
        super().__init__(exe=exe, exhaustiveness=exhaustiveness,
                         num_modes=num_modes, **kw)
        self.cnn_scoring = cnn_scoring
        self.cnn_model = cnn_model
        self.cnn_rotation = cnn_rotation
        self.min_rmsd_filter = min_rmsd_filter
        self.num_mc_saved = num_mc_saved

    def _build_cmd(self, receptor, ligand, out, logfile, box, flex, extra):
        cx, cy, cz = box.center
        sx, sy, sz = box.size
        cmd = [self.exe,
               "-r", receptor,
               "-l", ligand,
               "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}",
               "--center_z", f"{cz:.3f}",
               "--size_x", f"{sx:.3f}", "--size_y", f"{sy:.3f}",
               "--size_z", f"{sz:.3f}",
               "-o", out,
               "--log", logfile,
               "--cpu", str(self.ncpu),
               "--exhaustiveness", str(self.exhaustiveness),
               "--num_modes", str(self.num_modes),
               "--energy_range", str(self.energy_range),
               "--seed", str(self.seed),
               "--cnn_scoring", self.cnn_scoring,
               "--cnn", self.cnn_model,
               "--cnn_rotation", str(self.cnn_rotation),
               "--min_rmsd_filter", str(self.min_rmsd_filter),
               "--num_mc_saved", str(self.num_mc_saved),
               ]
        if flex:
            cmd += ["--flex", flex]
        cmd += list(extra)
        return cmd

    @staticmethod
    def parse_log(logfile: str) -> pd.DataFrame:
        """Parse GNINA log (mode, affinity, CNN pose score, CNN affinity)."""
        rows = []
        start = False
        with open(logfile) as fh:
            for line in fh:
                if "mode |" in line and "CNN" in line:
                    start = True
                    next(fh, None); next(fh, None)
                    continue
                if start:
                    s = line.split()
                    if len(s) >= 5:
                        try:
                            rows.append((int(s[0]),
                                         float(s[1]),    # Vina/minimized affinity
                                         float(s[2]),    # CNN pose score
                                         float(s[3])))   # CNN affinity
                        except ValueError:
                            break
        return pd.DataFrame(rows, columns=["mode", "affinity_kcal_mol",
                                           "cnn_pose_score", "cnn_affinity"])


# ---------------------------------------------------------------------------
# QVina-W / Vina-CUDA / AutoDock-GPU
# ---------------------------------------------------------------------------
class QVinaDocking(VinaDocking):
    """QuickVina 2 / QuickVina-W. Drop-in Vina replacement, faster search."""
    name = "qvina"
    def __init__(self, exe: str = "qvina-w", **kw):
        super().__init__(exe=exe, **kw)


class VinaCudaDocking(VinaDocking):
    """Vina-CUDA / Vina-GPU. Use ``gpu_ids`` for multi-GPU."""
    name = "vina-cuda"
    def __init__(self, exe: str = "vina_cuda", gpu_ids: str = "0", **kw):
        super().__init__(exe=exe, **kw)
        self.gpu_ids = gpu_ids

    def _build_cmd(self, receptor, ligand, out, logfile, box, flex, extra):
        cmd = super()._build_cmd(receptor, ligand, out, logfile, box, flex, extra)
        cmd += ["--gpu_id", self.gpu_ids]
        return cmd


class AutoDockGPUDocking(BaseDocking):
    """AutoDock-GPU (AD4 GPU docking). Requires grid maps (.fld / .map) prepared
    by ``AutoGrid`` or by MGLTools ``prepare_gpf.py``."""
    name = "adgpu"

    def __init__(self, exe: str = "autodock_gpu", **kw):
        super().__init__(exe=exe, **kw)

    def _build_cmd(self, receptor, ligand, out, logfile, box, flex, extra):
        # AutoDock-GPU takes a .fld file and a ligand pdbqt
        cmd = [self.exe,
               "--ffile", receptor,    # here receptor must be .fld
               "--lfile", ligand,
               "--nrun", str(self.exhaustiveness),
               "--devnum", str(self.extra.get("devnum", 0)),
               "--gbest",
               "--resnam", os.path.splitext(out)[0],
               ]
        cmd += list(extra)
        return cmd


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------
ENGINES = {
    "vina": VinaDocking,
    "smina": SminaDocking,
    "gnina": GninaDocking,
    "qvina": QVinaDocking,
    "vina-cuda": VinaCudaDocking,
    "adgpu": AutoDockGPUDocking,
}


def get_engine(name: str, **kw) -> BaseDocking:
    name = name.lower().strip()
    if name not in ENGINES:
        raise ValueError(f"Unknown docking engine {name}. Choose: {list(ENGINES)}")
    return ENGINES[name](**kw)


# ---------------------------------------------------------------------------
# Batch docking helper
# ---------------------------------------------------------------------------
def batch_dock(engine: BaseDocking,
               receptor: str,
               ligands: List[str],
               box: DockingBox,
               outdir: str = "docked",
               joblib_njobs: int = 1,
               **kw) -> pd.DataFrame:
    """Run docking for a list of ligands; return a concatenated DataFrame."""
    os.makedirs(outdir, exist_ok=True)
    from joblib import Parallel, delayed

    def _one(lig):
        base = os.path.splitext(os.path.basename(lig))[0]
        out = os.path.join(outdir, f"{base}_out.pdbqt")
        log = os.path.join(outdir, f"{base}.log")
        rc = engine.dock(receptor, lig, out=out, logfile=log, box=box, **kw)
        df = engine.parse_log(log) if os.path.exists(log) else pd.DataFrame()
        df["ligand"] = base
        df["out_pdbqt"] = out
        return df

    results = Parallel(n_jobs=joblib_njobs, backend="threading")(
        delayed(_one)(lig) for lig in ligands
    )
    return pd.concat([r for r in results if not r.empty], ignore_index=True)


# ---------------------------------------------------------------------------
# Preparation: PDB -> PDBQT via Meeko / obabel fallback
# ---------------------------------------------------------------------------
def prepare_receptor_pdbqt(pdb: str, out_pdbqt: Optional[str] = None,
                           method: str = "meeko") -> str:
    """Prepare a receptor PDBQT from a PDB file.

    method="meeko"  uses the modern `meeko` package (RDKit-based, MGLTools-free).
    method="obabel" falls back to OpenBabel.
    """
    out_pdbqt = out_pdbqt or os.path.splitext(pdb)[0] + ".pdbqt"
    if method == "meeko":
        try:
            from meeko import MoleculePreparation, PDBQTMolecule
            from meeko import ReceptorPreparation
            # meeko>=0.5
            ReceptorPreparation.from_pdb(pdb, output_pdbqt_filename=out_pdbqt)
            return out_pdbqt
        except Exception as e:
            log.warning("meeko receptor prep failed (%s), falling back to obabel", e)
    # openbabel fallback
    cmd = ["obabel", pdb, "-O", out_pdbqt, "-xr", "-h"]
    sp.check_call(cmd)
    return out_pdbqt


def prepare_ligand_pdbqt(ligand: str, out_pdbqt: Optional[str] = None,
                         pH: float = 7.4, method: str = "meeko") -> str:
    """Prepare ligand PDBQT (add hydrogens at *pH*, assign Gasteiger charges).

    method="meeko"     preferred: RDKit protonation + Meeko PDBQT writer
    method="obabel"    OpenBabel (obabel -h --partialcharge gasteiger)
    """
    out_pdbqt = out_pdbqt or os.path.splitext(ligand)[0] + ".pdbqt"
    if method == "meeko":
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from meeko import MoleculePreparation, PDBQTWriterLegacy
            # Protonate at given pH using RDKit Chem.MolFromSmiles/MolFromXYZ fallback
            if ligand.endswith((".sdf", ".mol2")):
                mol = Chem.MolFromMolFile(ligand, removeHs=False) if ligand.endswith(".mol") \
                    else Chem.MolFromMol2File(ligand, removeHs=False)
            else:
                mol = Chem.MolFromPDBFile(ligand, removeHs=False)
            if mol is None:
                raise RuntimeError("RDKit failed to read ligand")
            mol = Chem.AddHs(mol, addCoords=True)
            # basic 3D if missing
            if mol.GetNumConformers() == 0:
                AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                AllChem.MMFFOptimizeMolecule(mol)
            preparator = MoleculePreparation()
            mol_setup = preparator.prepare(mol)[0]
            pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(mol_setup)
            if not is_ok:
                raise RuntimeError(f"meeko pdbqt write errors: {err}")
            with open(out_pdbqt, "w") as f:
                f.write(pdbqt_str)
            return out_pdbqt
        except Exception as e:
            log.warning("meeko ligand prep failed (%s), falling back to obabel", e)
    cmd = ["obabel", ligand, "-O", out_pdbqt, "-h", "-p", str(pH),
           "--partialcharge", "gasteiger"]
    sp.check_call(cmd)
    return out_pdbqt
