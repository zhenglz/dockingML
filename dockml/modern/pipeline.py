#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
End-to-end docking + rescoring pipeline.
========================================

Steps:
    1. **Prepare**: protonate receptor (pdbfixer), prepare PDBQT.
       Prepare ligands (RDKit protonation at pH 7.4, 3D embedding, PDBQT).
    2. **Dock**: run the selected engine (vina/smina/gnina/...) for each ligand.
    3. **Pose extraction**: split multi-pose PDBQT into single-pose PDBs (mdtraj),
       optionally top-N per ligand.
    4. **Featurise**: compute per-residue contacts/OnionNet/PLIF features for
       each pose.
    5. **Rescore**: apply trained ML/DL model or consensus.
    6. **Report**: ranked ligand list with final score.

Example
-------
>>> from dockml.modern.pipeline import DockingPipeline
>>> p = DockingPipeline(engine="gnina", receptor="rec.pdb", ligand_resname="LIG")
>>> p.prepare_ligands(["lig1.sdf", "lig2.sdf"], pH=7.4)
>>> p.run_docking(box=DockingBox.from_ligand("cocrystallized.sdf"))
>>> df = p.featurise_and_rescore(model_path="model.pkl")
>>> df.to_csv("ranked.csv")
"""
from __future__ import annotations

import os
import logging
import tempfile
import shutil
from pathlib import Path
from typing import List, Optional, Dict, Union

import numpy as np
import pandas as pd

from .docking import (
    BaseDocking, VinaDocking, GninaDocking, SminaDocking,
    get_engine, DockingBox,
    prepare_receptor_pdbqt, prepare_ligand_pdbqt,
)
from .features.contacts import ResidueInteractionFeatures, OnionNetContacts
from .features.interactions import InteractionFingerprints
from .scoring.classical import ClassicalRescorer

log = logging.getLogger("dockml.modern.pipeline")


class DockingPipeline:
    """One-stop-shop for virtual-screening rescoring."""

    def __init__(self,
                 engine: Union[str, BaseDocking] = "vina",
                 receptor_pdb: str = "receptor.pdb",
                 ligand_resname: str = "LIG",
                 workdir: str = "dockingML_run",
                 exhaustiveness: int = 32,
                 ncpu: int = 8,
                 engine_kwargs: Optional[Dict] = None,
                 ):
        self.workdir = Path(workdir)
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.receptor_pdb = Path(receptor_pdb)
        self.ligand_resname = ligand_resname
        # prepare receptor
        self.receptor_pdbqt = self.workdir / "receptor.pdbqt"
        prepare_receptor_pdbqt(str(self.receptor_pdb), str(self.receptor_pdbqt))
        if isinstance(engine, str):
            kw = dict(ncpu=ncpu, exhaustiveness=exhaustiveness)
            if engine_kwargs:
                kw.update(engine_kwargs)
            self.engine = get_engine(engine, **kw)
        else:
            self.engine = engine
        self.ligand_pdbqt: Dict[str, Path] = {}
        self.docked: Dict[str, Path] = {}
        self.poses_df: Optional[pd.DataFrame] = None

    # ------------------------------------------------------------------
    def prepare_ligands(self, ligands: List[str], pH: float = 7.4):
        outdir = self.workdir / "ligands_pdbqt"
        outdir.mkdir(exist_ok=True)
        for lig in ligands:
            name = Path(lig).stem
            out = outdir / f"{name}.pdbqt"
            prepare_ligand_pdbqt(lig, str(out), pH=pH)
            self.ligand_pdbqt[name] = out
        log.info("Prepared %d ligands", len(self.ligand_pdbqt))

    # ------------------------------------------------------------------
    def run_docking(self, box: DockingBox, topn: int = 10,
                    flex_residues: Optional[str] = None):
        outdir = self.workdir / "docked"
        outdir.mkdir(exist_ok=True)
        all_rows = []
        for name, lig_pdbqt in self.ligand_pdbqt.items():
            out = outdir / f"{name}_out.pdbqt"
            logf = outdir / f"{name}.log"
            rc = self.engine.dock(str(self.receptor_pdbqt), str(lig_pdbqt),
                                  out=str(out), logfile=str(logf),
                                  box=box, flex=flex_residues)
            if rc != 0:
                log.warning("Docking failed for %s", name)
                continue
            self.docked[name] = out
            df = self.engine.parse_log(str(logf))
            df["ligand"] = name
            df = df.head(topn)
            all_rows.append(df)
        self.poses_df = pd.concat(all_rows, ignore_index=True) if all_rows else pd.DataFrame()
        log.info("Docked %d ligands, total %d poses", len(self.docked), len(self.poses_df))
        return self.poses_df

    # ------------------------------------------------------------------
    def featurise_poses(self,
                        use_plif: bool = True,
                        use_onion: bool = True,
                        split_pdb_dir: Optional[str] = None,
                        ) -> pd.DataFrame:
        """For each docked pose build a receptor-ligand complex PDB and compute
        feature vectors. Returns DataFrame indexed by (ligand, mode)."""
        import mdtraj as mt
        split_dir = Path(split_pdb_dir or (self.workdir / "poses_complexes"))
        split_dir.mkdir(exist_ok=True)

        rec_traj = mt.load(str(self.receptor_pdb))
        vecs, names_meta = [], []
        feature_names = None

        # NOTE: mdtraj cannot load PDBQT directly; convert ligand poses to PDB
        # with obabel first.
        from dockml.convert import Convert
        conv = Convert(obabel="obabel")
        rows = []

        for _, row in self.poses_df.iterrows():
            lig_name = row["ligand"]
            mode = int(row["mode"])
            docked_pdbqt = self.docked[lig_name]
            lig_pdb = split_dir / f"{lig_name}_m{mode}.pdb"
            # split poses and select this mode (obabel can split pdbqt w/ -m)
            conv.convert(str(docked_pdbqt), str(lig_pdb), verbose=False)
            # Build complex PDB by concatenating receptor + ligand
            complex_pdb = split_dir / f"{lig_name}_m{mode}_complex.pdb"
            self._concat_pdb(self.receptor_pdb, lig_pdb, complex_pdb,
                             lig_resname=self.ligand_resname)
            # Features
            vec_parts = []
            col_parts = []
            try:
                # OnionNet multi-shell contacts
                if use_onion:
                    onion = OnionNetContacts.from_pdb(str(complex_pdb),
                                                      ligand_resname=self.ligand_resname)
                    df_onion = onion.compute()
                    vec_parts.append(df_onion.values.sum(axis=0).ravel())
                    col_parts += list(df_onion.columns)
                # PLIF
                if use_plif:
                    plif = InteractionFingerprints(str(complex_pdb),
                                                   ligand_resname=self.ligand_resname)
                    v, n = plif.flat_vector()
                    vec_parts.append(v)
                    col_parts += n
                # Legacy backbone/sidechain contacts+vdw+coul
                rf = ResidueInteractionFeatures.from_pdb(str(complex_pdb),
                                                         ligand_resname=self.ligand_resname)
                v_leg, n_leg = rf.flat_feature_vector()
                vec_parts.append(v_leg)
                col_parts += n_leg
            except Exception as e:
                log.warning("Featurisation failed for %s pose %d: %s", lig_name, mode, e)
                continue
            vec = np.concatenate(vec_parts)
            rows.append((lig_name, mode, row.get("affinity_kcal_mol", np.nan),
                         row.get("cnn_affinity", np.nan), vec))
            feature_names = col_parts

        if not rows:
            return pd.DataFrame()
        X = np.vstack([r[4] for r in rows])
        df = pd.DataFrame(X, columns=feature_names)
        df.insert(0, "mode", [r[1] for r in rows])
        df.insert(0, "ligand", [r[0] for r in rows])
        df["affinity_vina"] = [r[2] for r in rows]
        df["cnn_affinity"] = [r[3] for r in rows]
        return df

    # ------------------------------------------------------------------
    def rescore(self, features: pd.DataFrame, model_path: str,
                include_engine_score: bool = True) -> pd.DataFrame:
        model = ClassicalRescorer.load(model_path)
        feat_cols = [c for c in features.columns if c not in
                     ("ligand", "mode", "affinity_vina", "cnn_affinity")]
        X = features[feat_cols].fillna(0.0).values
        proba = model.predict_proba(X)
        features = features.copy()
        features["ml_score"] = proba
        if include_engine_score and "affinity_vina" in features.columns:
            # rank combine
            features["rank_ml"] = features["ml_score"].rank(ascending=False)
            features["rank_vina"] = features["affinity_vina"].rank(ascending=True)
            features["final_score"] = 0.5 * features["ml_score"] + \
                0.5 * (1.0 - features["rank_vina"] / len(features))
        return features.sort_values("final_score" if "final_score" in features else "ml_score",
                                    ascending=False)

    # ------------------------------------------------------------------
    @staticmethod
    def _concat_pdb(rec_pdb: Path, lig_pdb: Path, out_pdb: Path,
                    lig_resname: str = "LIG"):
        """Concatenate ATOM/HETATM lines of receptor (PDB) and ligand (PDB)
        into a single complex PDB."""
        with open(out_pdb, "w") as out:
            for f in (rec_pdb, lig_pdb):
                with open(f) as fh:
                    for line in fh:
                        if line.startswith(("ATOM", "HETATM", "TER", "END")):
                            out.write(line)
            out.write("END\n")
