#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
CLI for topology building and MD running (modern).
"""
from __future__ import annotations
import argparse
import logging
from pathlib import Path


def _setup():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s :: %(message)s")


def gentop_main(argv=None):
    p = argparse.ArgumentParser(prog="dml-mdtop",
                                 description="Build modern protein-ligand topology "
                                             "(OpenFF + Parmed or AmberTools legacy).")
    p.add_argument("-r", "--receptor", required=True, help="Receptor PDB file")
    p.add_argument("-l", "--ligand", default=None, help="Ligand file (sdf/mol2/smi)")
    p.add_argument("--preset", default="amber19sb-opc",
                   choices=["amber19sb-opc", "amber14sb-tip3p", "charmm36m-tip3p"])
    p.add_argument("--padding", type=float, default=10.0,
                   help="Solvation box padding (Å)")
    p.add_argument("--ionic-strength", type=float, default=0.15)
    p.add_argument("-o", "--prefix", default="topol", help="Output file prefix")
    p.add_argument("--repair", action="store_true", help="Repair receptor with pdbfixer")
    p.add_argument("--legacy", action="store_true", help="Use legacy AmberTools path")
    p.add_argument("--amberhome", default=None, help="AMBERHOME (legacy path)")
    args = p.parse_args(argv)
    _setup()

    from automd.topology import ProteinLigandSystem, repair_protein_pdb
    rec = Path(args.receptor)
    if args.repair:
        fixed = rec.with_name(rec.stem + "_fixed.pdb")
        repair_protein_pdb(str(rec), str(fixed))
        rec = fixed
    if args.legacy:
        sys = ProteinLigandSystem(preset=args.preset)
        sys._protein_pdb = str(rec)
        sys._ligand_file = args.ligand
        sys._build_legacy()
        sys.write_gromacs_legacy(args.prefix, amberhome=args.amberhome or "",
                                  box_edge=args.padding)
    else:
        sys = ProteinLigandSystem.from_files(str(rec), args.ligand, preset=args.preset)
        sys.solvate_box(edge=args.padding, ion_strength=args.ionic_strength)
        sys.write_gromacs(f"{args.prefix}.top", f"{args.prefix}.gro", posre=True)


def md_main(argv=None):
    p = argparse.ArgumentParser(prog="dml-md",
                                 description="Run the modern GROMACS EM -> NVT -> NPT "
                                             "equilibration and production pipeline.")
    p.add_argument("-c", "--gro", required=True, help="Starting coordinates (.gro)")
    p.add_argument("-p", "--top", required=True, help="Topology (.top)")
    p.add_argument("--preset", default="amber19sb-opc")
    p.add_argument("--length-ns", type=float, default=100.0)
    p.add_argument("--dt-ps", type=float, default=0.002)
    p.add_argument("--ncpu", type=int, default=8)
    p.add_argument("--gpu-ids", default="", help="GPU ids, e.g. '0' or '0,1'")
    p.add_argument("--ntmpi", type=int, default=1)
    p.add_argument("--ligand-resname", default="LIG")
    p.add_argument("--prefix", default="md")
    p.add_argument("--gmx", default="gmx")
    args = p.parse_args(argv)
    _setup()

    from automd.mdrun import GromacsRunner
    runner = GromacsRunner(preset=args.preset, gmx=args.gmx, ncpu=args.ncpu,
                           gpu_ids=args.gpu_ids, ntmpi=args.ntmpi)
    runner.equilibrate_produce(args.gro, args.top, length_ns=args.length_ns,
                               dt_ps=args.dt_ps,
                               ligand_resname=args.ligand_resname,
                               out_prefix=args.prefix)


if __name__ == "__main__":
    md_main()
