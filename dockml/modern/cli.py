#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Command-line interface for the modern dockingML pipeline.

Subcommands
-----------
    prepare   -- prepare receptor/ligand PDBQT files
    dock      -- run docking with selected engine
    features  -- compute interaction features
    train     -- train a classical ML rescorer
    rescore   -- rescore existing poses with a trained model
    pipeline  -- run the full prepare->dock->featurise->rescore workflow
"""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger("dockml.modern.cli")


def _setup_logging(level=logging.INFO):
    logging.basicConfig(level=level,
                        format="%(asctime)s %(levelname)s %(name)s :: %(message)s")


# ---------------------------------------------------------------------------
def prepare_main(argv=None):
    p = argparse.ArgumentParser(prog="dml-prepare")
    p.add_argument("-r", "--receptor", required=True, help="Input receptor PDB")
    p.add_argument("-l", "--ligands", nargs="*", default=[],
                   help="Ligand files (sdf/mol2/pdb/smi)")
    p.add_argument("-o", "--outdir", default="prepared")
    p.add_argument("--pH", type=float, default=7.4)
    p.add_argument("--repair", action="store_true",
                   help="Repair receptor with pdbfixer (add missing atoms/H)")
    args = p.parse_args(argv)
    _setup_logging()

    from .docking import prepare_receptor_pdbqt, prepare_ligand_pdbqt
    from automd.topology import repair_protein_pdb

    outdir = Path(args.outdir); outdir.mkdir(exist_ok=True)
    rec = Path(args.receptor)
    if args.repair:
        fixed = outdir / (rec.stem + "_fixed.pdb")
        repair_protein_pdb(str(rec), str(fixed), pH=args.pH)
        rec = fixed
    rec_out = outdir / (rec.stem + ".pdbqt")
    prepare_receptor_pdbqt(str(rec), str(rec_out))
    print(f"[prepare] receptor PDBQT -> {rec_out}")
    for lig in args.ligands:
        out = outdir / (Path(lig).stem + ".pdbqt")
        prepare_ligand_pdbqt(lig, str(out), pH=args.pH)
        print(f"[prepare] ligand PDBQT  -> {out}")


# ---------------------------------------------------------------------------
def dock_main(argv=None):
    p = argparse.ArgumentParser(prog="dml-dock")
    p.add_argument("-r", "--receptor", required=True)
    p.add_argument("-l", "--ligand_dir", default=None,
                   help="Directory of PDBQT ligands; if omitted, ligands must be "
                        "listed positionally.")
    p.add_argument("ligands", nargs="*", help="Ligand PDBQT files")
    p.add_argument("--engine", default="vina",
                   choices=["vina", "smina", "gnina", "qvina", "vina-cuda"])
    p.add_argument("--center", nargs=3, type=float, required=True,
                   metavar=("CX", "CY", "CZ"), help="Binding box center (Å)")
    p.add_argument("--size", nargs=3, type=float, default=[22, 22, 22])
    p.add_argument("--exhaustiveness", type=int, default=32)
    p.add_argument("--ncpu", type=int, default=8)
    p.add_argument("--outdir", default="docked")
    p.add_argument("--scoring", default=None, help="Override scoring function (vina/vinardo/ad4/...)")
    args = p.parse_args(argv)
    _setup_logging()

    from .docking import get_engine, DockingBox, batch_dock
    ligs = args.ligands
    if args.ligand_dir:
        ligs += [str(p) for p in Path(args.ligand_dir).glob("*.pdbqt")]
    engine = get_engine(args.engine, ncpu=args.ncpu,
                        exhaustiveness=args.exhaustiveness,
                        **({"scoring": args.scoring} if args.scoring else {}))
    box = DockingBox(center=tuple(args.center), size=tuple(args.size))
    df = batch_dock(engine, args.receptor, ligs, box, outdir=args.outdir,
                    joblib_njobs=1)
    df.to_csv(Path(args.outdir) / "docking_summary.csv", index=False)
    print(df.head(20).to_string())


# ---------------------------------------------------------------------------
def features_main(argv=None):
    p = argparse.ArgumentParser(prog="dml-features")
    p.add_argument("-i", "--input", required=True,
                   help="Two-column text file: <complex_pdb> <lig_resname>")
    p.add_argument("-o", "--output", default="features.csv")
    p.add_argument("--onion", action="store_true", help="OnionNet multi-shell contacts")
    p.add_argument("--plif", action="store_true", help="Interaction fingerprints")
    args = p.parse_args(argv)
    _setup_logging()

    from .features.contacts import ResidueInteractionFeatures, OnionNetContacts
    from .features.interactions import InteractionFingerprints

    records = []
    vecs = []
    cols = None
    with open(args.input) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            pdb, lig = parts[0], parts[1]
            print(f"[features] processing {pdb} ({lig})")
            parts_v = []
            parts_n = []
            rf = ResidueInteractionFeatures.from_pdb(pdb, lig)
            v, n = rf.flat_feature_vector()
            parts_v.append(v); parts_n += n
            if args.onion:
                oc = OnionNetContacts.from_pdb(pdb, lig).compute()
                parts_v.append(oc.values.sum(0).ravel())
                parts_n += list(oc.columns)
            if args.plif:
                v2, n2 = InteractionFingerprints(pdb, lig).flat_vector()
                parts_v.append(v2); parts_n += n2
            vec = np.concatenate(parts_v)
            vecs.append(vec)
            records.append((pdb, lig))
            cols = parts_n
    X = np.vstack(vecs)
    df = pd.DataFrame(X, columns=cols)
    df.insert(0, "lig_resname", [r[1] for r in records])
    df.insert(0, "pdb", [r[0] for r in records])
    df.to_csv(args.output, index=False)
    print(f"[features] wrote {args.output} with shape {df.shape}")


# ---------------------------------------------------------------------------
def train_main(argv=None):
    p = argparse.ArgumentParser(prog="dml-train")
    p.add_argument("-X", "--features", required=True,
                   help="CSV of features (rows = samples)")
    p.add_argument("-y", "--labels", required=True,
                   help="CSV containing binary/regression labels; must include a "
                        "`label` column aligned 1:1 with -X rows.")
    p.add_argument("-m", "--model", default="xgb",
                   choices=["rf", "xgb", "lgbm", "lr", "ridge", "gbt"])
    p.add_argument("--task", default="classify", choices=["classify", "regress"])
    p.add_argument("-o", "--output", default="model.pkl")
    p.add_argument("--cv", type=int, default=5, help="N-fold CV reporting")
    args = p.parse_args(argv)
    _setup_logging()

    from .scoring.classical import ClassicalRescorer
    X = pd.read_csv(args.features)
    y = pd.read_csv(args.labels)
    drop = [c for c in X.columns if c in ("pdb", "lig_resname", "label")]
    Xf = X.drop(columns=drop)
    if "label" in y.columns:
        yv = y["label"].values
    else:
        yv = y.iloc[:, -1].values
    model = ClassicalRescorer(model=args.model, task=args.task, scale=(args.model in ("lr", "ridge")))
    metrics = model.cv_evaluate(Xf.values, yv, n_fold=args.cv)
    print("[train] CV metrics:", metrics)
    model.fit(Xf.values, yv)
    model.save(args.output)
    print(f"[train] model saved to {args.output}")


# ---------------------------------------------------------------------------
def rescore_main(argv=None):
    p = argparse.ArgumentParser(prog="dml-rescore")
    p.add_argument("-f", "--features", required=True, help="Features CSV (from dml-features)")
    p.add_argument("-m", "--model", required=True, help="Trained model .pkl")
    p.add_argument("-o", "--output", default="rescored.csv")
    args = p.parse_args(argv)
    _setup_logging()

    from .scoring.classical import ClassicalRescorer
    df = pd.read_csv(args.features)
    model = ClassicalRescorer.load(args.model)
    drop = [c for c in df.columns if c in ("pdb", "lig_resname", "label")]
    X = df.drop(columns=drop).fillna(0).values
    if model.task == "classify":
        df["score"] = model.predict_proba(X)
    else:
        df["score"] = model.predict(X)
    df = df.sort_values("score", ascending=False)
    df.to_csv(args.output, index=False)
    print(f"[rescore] wrote ranked results to {args.output}")


# ---------------------------------------------------------------------------
def pipeline_main(argv=None):
    p = argparse.ArgumentParser(prog="dml-pipeline")
    p.add_argument("-r", "--receptor", required=True)
    p.add_argument("-l", "--ligands", nargs="+", required=True)
    p.add_argument("--engine", default="vina")
    p.add_argument("--center", nargs=3, type=float, required=True)
    p.add_argument("--size", nargs=3, type=float, default=[22, 22, 22])
    p.add_argument("--ligand_resname", default="LIG")
    p.add_argument("--model", default=None, help="Pre-trained model pkl; "
                   "if omitted, only docking is performed.")
    p.add_argument("--workdir", default="dml_workdir")
    args = p.parse_args(argv)
    _setup_logging()
    from .pipeline import DockingPipeline
    from .docking import DockingBox
    pipe = DockingPipeline(engine=args.engine, receptor_pdb=args.receptor,
                           ligand_resname=args.ligand_resname, workdir=args.workdir)
    pipe.prepare_ligands(args.ligands)
    box = DockingBox(center=tuple(args.center), size=tuple(args.size))
    pipe.run_docking(box)
    feats = pipe.featurise_poses()
    if args.model:
        ranked = pipe.rescore(feats, args.model)
        ranked.to_csv(Path(args.workdir) / "ranked.csv", index=False)
        print(ranked[["ligand", "mode", "ml_score", "final_score"]].head(20).to_string())
    else:
        feats.to_csv(Path(args.workdir) / "features.csv", index=False)


if __name__ == "__main__":
    pipeline_main()
