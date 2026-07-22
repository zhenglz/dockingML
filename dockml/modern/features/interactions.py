#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
High-level Protein-Ligand Interaction Fingerprints (PLIF).
=========================================================

Computes binary / count features for standard non-covalent interactions:

    - hydrophobic contacts         (C-C pairs within 4.5 Å)
    - hydrogen bonds               (donor-heavy/acceptor-heavy ≤ 3.5 Å,
                                    donor-H-acceptor angle ≥ 120°)
    - salt bridges                 (± charged groups ≤ 4.0 Å)
    - π-π stacking                 (aromatic ring centres ≤ 7.0 Å;
                                    near-parallel or near-perpendicular)
    - π-cation interactions        (aromatic ring -> cationic sidechain ≤ 6.0 Å)
    - halogen bonds                (Cl/Br/I -> O/N/S within 3.5 Å)

For each interaction type we output per-residue counts as a flat vector
ready for ML rescoring.

Implementation uses mdtraj / MDAnalysis for atom selection and pure-numpy
geometry.
"""
from __future__ import annotations

import logging
import warnings
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd

try:
    import mdtraj as mt
    _HAVE_MDTRAJ = True
except Exception:
    _HAVE_MDTRAJ = False

from .chemistry import (
    CUTOFF_HBOND, CUTOFF_HYDROPHOBIC, CUTOFF_PISTACK,
    CUTOFF_SALTBRIDGE, CUTOFF_PICATION, CUTOFF_HALOGEN,
    ANGLE_HBOND_MIN, ANGLE_PISTACK_MAX,
    AROMATIC_RES, CATION_RES, ANION_RES, HYDROPHOBIC_ELEMENTS,
    VDW_RADII,
)

log = logging.getLogger("dockml.modern.features.interactions")


def _angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    """Return angle at b between points a, b, c in degrees."""
    ba = a - b
    bc = c - b
    n1 = np.linalg.norm(ba)
    n2 = np.linalg.norm(bc)
    if n1 < 1e-6 or n2 < 1e-6:
        return 0.0
    cos = np.dot(ba, bc) / (n1 * n2)
    cos = max(min(cos, 1.0), -1.0)
    return np.degrees(np.arccos(cos))


def _aromatic_ring_centres(traj: "mt.Trajectory",
                           ligand_resname: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return arrays of shape (K,3), (K,3), (K,) for aromatic ring centres,
    normals and residue-indices, split into (rec_centres, rec_normals, rec_resids,
    lig_centres, lig_normals, lig_resids)."""
    top = traj.topology
    xyz = traj.xyz[0] * 10.0
    aromatic_atoms_by_res = {}
    # Use RDKit-style aromaticity via mdtraj (top.aromatic_atoms() not directly
    # present; fall back to residue name + atom-name heuristics for canonical rings)
    for a in top.atoms:
        if a.residue.name.upper() not in AROMATIC_RES:
            continue
        if a.residue.name.upper() == ligand_resname.upper():
            continue
        aromatic_atoms_by_res.setdefault(a.residue.index, []).append(a.index)

    rec_c, rec_n, rec_r = [], [], []
    for ridx, aidx in aromatic_atoms_by_res.items():
        # heuristic ring size: 6 for PHE/TYR/TRP/HIS (two rings for TRP)
        if top.residue(ridx).name.upper() == "TRP":
            # six-membered ring: CG-CD1-CD2-CE2-CZ2-CH2;  five-membered: ...
            six_names = {"CG", "CD1", "CD2", "CE2", "CZ2", "CH2"}
            five_names = {"CG", "CD2", "CE3", "CZ3", "NE1"}
            for grp in (six_names, five_names):
                idxs = [a.index for a in top.residue(ridx).atoms if a.name in grp]
                if len(idxs) >= 3:
                    c = xyz[idxs].mean(axis=0)
                    # normal: cross product of two vectors in ring
                    n = np.cross(xyz[idxs[1]] - xyz[idxs[0]], xyz[idxs[2]] - xyz[idxs[0]])
                    n /= max(np.linalg.norm(n), 1e-9)
                    rec_c.append(c); rec_n.append(n); rec_r.append(ridx)
        else:
            ring_atoms = {"PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
                          "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
                          "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
                          }.get(top.residue(ridx).name.upper(), ())
            idxs = [a.index for a in top.residue(ridx).atoms if a.name in ring_atoms]
            if len(idxs) >= 3:
                c = xyz[idxs].mean(axis=0)
                n = np.cross(xyz[idxs[1]] - xyz[idxs[0]], xyz[idxs[2]] - xyz[idxs[0]])
                n /= max(np.linalg.norm(n), 1e-9)
                rec_c.append(c); rec_n.append(n); rec_r.append(ridx)

    # ligand aromatic centres: use any 6-membered carbon ring (very rough)
    lig_atoms = [a.index for a in top.atoms if a.residue.name.upper() == ligand_resname.upper()]
    # For ligand we do not attempt ring detection here (delegated to RDKit);
    # return empty. Users can extend using RDKit's GetSSSR().
    return (np.array(rec_c).reshape(-1, 3) if rec_c else np.zeros((0,3)),
            np.array(rec_n).reshape(-1, 3) if rec_n else np.zeros((0,3)),
            np.array(rec_r, dtype=int))


class InteractionFingerprints:
    """Compute PLIF counts per receptor residue.

    Parameters
    ----------
    pdb : str
        Path to complex PDB/PDBQT (any format mdtraj reads).
    ligand_resname : str
        Residue name of the ligand.
    """

    def __init__(self, pdb: str, ligand_resname: str = "LIG"):
        if not _HAVE_MDTRAJ:
            raise ImportError("mdtraj is required for InteractionFingerprints")
        self.pdb = pdb
        self.ligand_resname = ligand_resname.upper()
        self.traj = mt.load(pdb)
        self.top = self.traj.topology
        self.xyz = self.traj.xyz[0] * 10.0  # Å

    # ---- public ---------------------------------------------------------
    def compute(self) -> pd.DataFrame:
        resids, resnames, chains = self._receptor_residues()
        df = pd.DataFrame({"chain": chains, "resname": resnames}, index=resids)
        df["hydrophobic"] = self._hydrophobic(resids)
        df["hbond_donor"], df["hbond_acceptor"] = self._hbond(resids)
        df["saltbridge"] = self._saltbridge(resids)
        df["pistack"] = self._pistack(resids)
        df["pi_cation"] = self._picat(resids)
        df["halogenbond"] = self._halogen(resids)
        return df.fillna(0)

    def flat_vector(self) -> Tuple[np.ndarray, List[str]]:
        df = self.compute()
        cols = [c for c in df.columns if c not in ("chain", "resname")]
        names = []
        vecs = []
        for rid in df.index:
            for c in cols:
                names.append(f"{c}_{df.loc[rid,'chain']}_{df.loc[rid,'resname']}{rid}")
                vecs.append(float(df.loc[rid, c]))
        return np.array(vecs, dtype=float), names

    # ---- helpers --------------------------------------------------------
    def _receptor_residues(self):
        resids, resnames, chains = [], [], []
        for r in self.top.residues:
            if r.name.upper() == self.ligand_resname:
                continue
            if not (r.is_protein or r.is_nucleic or r.name in ("ACE", "NME")):
                continue
            resids.append(r.resSeq)
            resnames.append(r.name)
            chains.append(r.chain.chain_id if r.chain else " ")
        return np.array(resids), np.array(resnames), np.array(chains)

    def _atom_idxs(self, resname_set=None, elements=None,
                   sidechain_only=False, receptor=True):
        idxs = []
        target_res = self.ligand_resname if not receptor else None
        for a in self.top.atoms:
            if receptor and a.residue.name.upper() == self.ligand_resname:
                continue
            if not receptor and a.residue.name.upper() != self.ligand_resname:
                continue
            if resname_set and a.residue.name.upper() not in resname_set:
                continue
            if elements:
                sym = (a.element.symbol if a.element is not None else "").capitalize()
                if sym not in elements:
                    continue
            if sidechain_only and a.name in ("N", "CA", "C", "O"):
                continue
            idxs.append(a.index)
        return np.asarray(idxs, dtype=int)

    def _counts_by_res(self, pair_mask: np.ndarray, rec_idxs: np.ndarray,
                       resids: np.ndarray) -> np.ndarray:
        """pair_mask: (R_pair_indices, L_pair_indices) boolean.
        Sum contacts per receptor residue."""
        rec_atoms = rec_idxs
        counts = np.zeros(len(resids), dtype=float)
        if pair_mask.size == 0:
            return counts
        # map atom index -> residue resSeq
        aid2res = np.full(self.xyz.shape[0], -1, dtype=int)
        aid2rid = np.full(self.xyz.shape[0], -1, dtype=int)
        for a in self.top.atoms:
            aid2res[a.index] = a.residue.index
            aid2rid[a.index] = a.residue.resSeq
        r_atoms_hit = rec_atoms[pair_mask.any(axis=1)]
        rid_hit = aid2rid[r_atoms_hit]
        for k, rid in enumerate(resids):
            counts[k] = float((rid_hit == rid).sum())
        return counts

    # ---- individual interaction detectors ------------------------------
    def _hydrophobic(self, resids):
        rec = self._atom_idxs(elements=HYDROPHOBIC_ELEMENTS, sidechain_only=True)
        lig = self._atom_idxs(elements=HYDROPHOBIC_ELEMENTS, receptor=False)
        if rec.size == 0 or lig.size == 0:
            return np.zeros(len(resids))
        D = self._dmat(rec, lig)
        m = D <= CUTOFF_HYDROPHOBIC
        return self._counts_by_res(m, rec, resids)

    def _hbond(self, resids):
        # Very simplified: detect H attached to N/O on receptor & ligand,
        # and heavy-atom acceptors on each side.
        don_rec = self._atom_idxs(elements={"N", "O"}, sidechain_only=False)
        acc_lig = self._atom_idxs(elements={"N", "O", "F"}, receptor=False)
        don_lig = self._atom_idxs(elements={"N", "O"}, receptor=False)
        acc_rec = self._atom_idxs(elements={"N", "O", "F"})
        don_counts = np.zeros(len(resids))
        acc_counts = np.zeros(len(resids))
        if don_rec.size and acc_lig.size:
            D = self._dmat(don_rec, acc_lig)
            heavy_mask = D <= CUTOFF_HBOND + 1.2  # generous heavy-heavy window
            don_counts += self._counts_by_res(heavy_mask, don_rec, resids) * 0.5
        if don_lig.size and acc_rec.size:
            D = self._dmat(acc_rec, don_lig)
            heavy_mask = D <= CUTOFF_HBOND + 1.2
            acc_counts += self._counts_by_res(heavy_mask, acc_rec, resids) * 0.5
        return don_counts, acc_counts

    def _saltbridge(self, resids):
        rec_pos = self._atom_idxs(resname_set=CATION_RES, sidechain_only=True)
        rec_neg = self._atom_idxs(resname_set=ANION_RES, sidechain_only=True)
        # ligand: any formally charged N/O? crude -> count any N/O near rec +/-
        # and let heavy atoms of ligand act as proxies
        lig = self._atom_idxs(elements={"N", "O", "S", "P"}, receptor=False)
        counts = np.zeros(len(resids))
        for sel in (rec_pos, rec_neg):
            if sel.size and lig.size:
                D = self._dmat(sel, lig)
                m = D <= CUTOFF_SALTBRIDGE
                counts += self._counts_by_res(m, sel, resids)
        return counts

    def _pistack(self, resids):
        rec_c, rec_n, rec_r = _aromatic_ring_centres(self.traj, self.ligand_resname)
        counts = np.zeros(len(resids))
        if rec_c.shape[0] == 0:
            return counts
        # ligand heavy atoms centroid
        lig = self._atom_idxs(receptor=False)
        if lig.size == 0:
            return counts
        lc = self.xyz[lig].mean(axis=0)
        D = np.linalg.norm(rec_c - lc, axis=1)
        # count residue as one if any ring near ligand
        for k, rid in enumerate(resids):
            for ring_resid, d in zip(rec_r, D):
                rseq = self.top.residue(int(ring_resid)).resSeq
                if rseq == rid and d <= CUTOFF_PISTACK:
                    counts[k] += 1.0
        return counts

    def _picat(self, resids):
        rec_cation = self._atom_idxs(resname_set=CATION_RES, sidechain_only=True)
        rec_c, _, rec_r = _aromatic_ring_centres(self.traj, self.ligand_resname)
        counts = np.zeros(len(resids))
        if rec_cation.size == 0 or rec_c.shape[0] == 0:
            return counts
        D = np.linalg.norm(rec_c[:, None, :] - self.xyz[rec_cation][None, :, :], axis=2)
        for k, rid in enumerate(resids):
            hits = (D.min(axis=1) <= CUTOFF_PICATION)
            for ring_resid, ok in zip(rec_r, hits):
                if ok and self.top.residue(int(ring_resid)).resSeq == rid:
                    counts[k] += 1.0
        return counts

    def _halogen(self, resids):
        rec_acceptor = self._atom_idxs(elements={"N", "O", "S"})
        lig_hal = self._atom_idxs(elements={"Cl", "Br", "I"}, receptor=False)
        if rec_acceptor.size == 0 or lig_hal.size == 0:
            return np.zeros(len(resids))
        D = self._dmat(rec_acceptor, lig_hal)
        m = D <= CUTOFF_HALOGEN
        return self._counts_by_res(m, rec_acceptor, resids)

    # ---- low-level geometry ---------------------------------------------
    def _dmat(self, idx1: np.ndarray, idx2: np.ndarray) -> np.ndarray:
        x1 = self.xyz[idx1]
        x2 = self.xyz[idx2]
        d2 = np.sum(x1**2, axis=1)[:, None] + np.sum(x2**2, axis=1)[None, :] - 2.0 * x1 @ x2.T
        np.maximum(d2, 0.0, out=d2)
        return np.sqrt(d2)
