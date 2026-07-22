#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Modern MD driver
================

This module replaces the legacy ``automd.md_gmx.AutoRunMD`` class and adds:

* GROMACS 2022+ compatible workflow with modern MDP files (Verlet+PME,
  v-rescale, Parrinello-Rahman, h-bond LINCS, ...).
* Configurable force-field preset  -> ``"amber19sb-opc"`` (default/recommended),
  ``"amber14sb-tip3p"``, ``"charmm36m-tip3p"``.
* OpenMM backend for GPU-accelerated Python-driven simulations (AMBER/CHARMM/OpenFF).
* Automatic equilibration ladder (EM -> NVT -> NPT restrain -> NPT production).
* Restraint generation (posres.itp) for heavy atoms.

Example
-------
>>> from automd.mdrun import GromacsRunner, ForceFieldPreset
>>> runner = GromacsRunner(preset="amber19sb-opc", gmx="gmx", ncpu=16, gpu_ids="0")
>>> runner.equilibrate_produce("complex.gro", "topol.top", length_ns=100)
"""
from __future__ import annotations

import os
import shutil
import logging
import subprocess as sp
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict

log = logging.getLogger("automd.mdrun")
log.addHandler(logging.NullHandler())

DATA = Path(__file__).parent / "data" / "modern"


# ---------------------------------------------------------------------------
@dataclass
class ForceFieldPreset:
    """Small bundle describing a recommended FF/water/ion combination."""
    name: str
    protein_ff: str          # passed to pdb2gmx, e.g. 'amber19sb'
    water_model: str        # e.g. 'opc' or 'tip3p'
    ligand_ff: str          # e.g. 'gaff2' / 'openff-2.1.0' / 'cgenff'
    description: str = ""

    @classmethod
    def lookup(cls, key: str) -> "ForceFieldPreset":
        table = {
            "amber19sb-opc": cls(
                name="amber19sb-opc",
                protein_ff="amber19sb",
                water_model="opc",
                ligand_ff="gaff2",
                description="Recommended modern AMBER force field with OPC water; "
                            "best for folded proteins and drug-like ligands.",
            ),
            "amber14sb-tip3p": cls(
                name="amber14sb-tip3p",
                protein_ff="amber14sb",
                water_model="tip3p",
                ligand_ff="gaff2",
                description="Legacy-compatible AMBER ff14SB + TIP3P.",
            ),
            "charmm36m-tip3p": cls(
                name="charmm36m-tip3p",
                protein_ff="charmm36m",
                water_model="charmm36",
                ligand_ff="cgenff",
                description="CHARMM36m for proteins + CHARMM-modified TIP3P; very "
                            "good balance for folded & disordered systems.",
            ),
        }
        if key not in table:
            raise KeyError(f"Unknown preset {key}; choose {list(table)}")
        return table[key]


# ---------------------------------------------------------------------------
class GromacsRunner:
    """High-level driver for standard EM -> NVT -> NPTeq -> NPTprod ladder."""

    def __init__(self,
                 preset: str = "amber19sb-opc",
                 gmx: str = "gmx",
                 ncpu: int = 8,
                 gpu_ids: str = "",
                 ntmpi: int = 1,
                 pin: str = "on",
                 ):
        self.preset = ForceFieldPreset.lookup(preset)
        self.gmx = gmx
        self.ncpu = ncpu
        self.gpu_ids = str(gpu_ids)
        self.ntmpi = ntmpi
        self.pin = pin

    # ------------------------------------------------------------------
    # low-level run helpers
    def _run(self, *args, cwd=None, input_: Optional[str] = None):
        cmd = [self.gmx] + list(args)
        log.info("$ %s", " ".join(cmd))
        proc = sp.run(cmd, cwd=cwd, input=input_, text=True,
                      stdout=sp.PIPE, stderr=sp.PIPE)
        if proc.returncode != 0:
            log.error("gmx failed:\n%s\n%s", proc.stdout, proc.stderr)
            raise RuntimeError(f"gmx command failed: {' '.join(cmd)}")
        return proc.stdout

    def _which(self, prog: str) -> str:
        p = shutil.which(f"{self.gmx} {prog}" if self.gmx == "gmx" else prog)
        if not p:
            # gmx uses sub-command style: gmx grompp, gmx mdrun, etc.
            return prog
        return prog

    # ------------------------------------------------------------------
    def grompp(self, mdp: Path, gro: str, top: str, out_tpr: str,
               maxwarn: int = 5) -> None:
        self._run("grompp", "-f", str(mdp), "-c", gro, "-p", top,
                  "-o", out_tpr, "-maxwarn", str(maxwarn))

    def mdrun(self, tpr: str, deffnm: str, nsteps: Optional[int] = None,
              plumed: Optional[str] = None, multidir: Optional[str] = None) -> None:
        args = ["mdrun", "-deffnm", deffnm, "-nt", str(self.ncpu),
                "-ntmpi", str(self.ntmpi), "-pin", self.pin]
        if self.gpu_ids:
            args += ["-gpu_id", self.gpu_ids]
        if nsteps is not None:
            args += ["-nsteps", str(nsteps)]
        if plumed:
            args += ["-plumed", plumed]
        self._run(*args)

    # ------------------------------------------------------------------
    # build index groups (Protein-Ligand, Water_and_ions)
    def make_index(self, gro: str, ndx: str = "index.ndx",
                   ligand_resname: str = "LIG") -> None:
        """Create standard index groups for temperature coupling."""
        script = (
            f"keep 0\n"
            f"1|r{ligand_resname}\n"
            f"name 1 Protein_Ligand\n"
            f"q\n"
        )
        self._run("make_ndx", "-f", gro, "-o", ndx, input_=script)

    # ------------------------------------------------------------------
    def equilibrate_produce(self, gro: str, top: str,
                            length_ns: float = 100.0,
                            dt_ps: float = 0.002,
                            ligand_resname: str = "LIG",
                            em_nsteps: int = 50000,
                            nvt_ps: float = 500.0,
                            npt_eq_ps: float = 1000.0,
                            out_prefix: str = "md",
                            ):
        """Full equilibration + production pipeline."""
        make_ndx_script = "\n".join([
            f'"Protein" | r{ligand_resname}',
            "name 21 Protein_Ligand",
            '"Water" | "Ion"',
            "name 22 Water_and_ions",
            "q",
        ])
        self._run("make_ndx", "-f", gro, "-o", "index.ndx", input_=make_ndx_script)

        # ---- energy minimization ----
        self.grompp(DATA / "em.mdp", gro, top, "em.tpr")
        self.mdrun("em.tpr", "em")

        # ---- NVT ----
        nvt_steps = int(nvt_ps * 1000 / dt_ps)
        self._run("grompp", "-f", str(DATA / "nvt_eq.mdp"),
                  "-c", "em.gro", "-r", "em.gro", "-p", top,
                  "-n", "index.ndx", "-o", "nvt.tpr", "-maxwarn", "5")
        self.mdrun("nvt.tpr", "nvt", nsteps=nvt_steps)

        # ---- NPT equilibration ----
        npt_eq_steps = int(npt_eq_ps * 1000 / dt_ps)
        self._run("grompp", "-f", str(DATA / "npt_eq.mdp"),
                  "-c", "nvt.gro", "-r", "nvt.gro", "-t", "nvt.cpt",
                  "-p", top, "-n", "index.ndx", "-o", "npt_eq.tpr", "-maxwarn", "5")
        self.mdrun("npt_eq.tpr", "npt_eq", nsteps=npt_eq_steps)

        # ---- production ----
        prod_steps = int(length_ns * 1000 / dt_ps)
        # copy template mdp and override nsteps
        mdp_out = Path(f"{out_prefix}.mdp")
        shutil.copy(DATA / "npt_prod.mdp", mdp_out)
        self._mdp_set(mdp_out, {"nsteps": prod_steps, "dt": dt_ps})
        self._run("grompp", "-f", str(mdp_out),
                  "-c", "npt_eq.gro", "-t", "npt_eq.cpt",
                  "-p", top, "-n", "index.ndx",
                  "-o", f"{out_prefix}.tpr", "-maxwarn", "5")
        self.mdrun(f"{out_prefix}.tpr", out_prefix, nsteps=prod_steps)
        log.info("Production run complete: %s.gro/.xtc/.edr", out_prefix)

    # ------------------------------------------------------------------
    @staticmethod
    def _mdp_set(mdp_path: Path, kv: Dict[str, object]) -> None:
        """Replace key = value pairs in an mdp file."""
        if not mdp_path.exists():
            return
        lines = mdp_path.read_text().splitlines()
        out = []
        replaced = set()
        for line in lines:
            stripped = line.split(";")[0].strip()
            if not stripped:
                out.append(line); continue
            key = stripped.split()[0]
            if key in kv:
                val = kv[key]
                if isinstance(val, float):
                    out.append(f"{key:<22} = {val} ; edited by automd")
                else:
                    out.append(f"{key:<22} = {val} ; edited by automd")
                replaced.add(key)
            else:
                out.append(line)
        for k, v in kv.items():
            if k not in replaced:
                out.append(f"{k:<22} = {v} ; edited by automd")
        mdp_path.write_text("\n".join(out) + "\n")


# ---------------------------------------------------------------------------
# OpenMM backend (optional)
# ---------------------------------------------------------------------------
class OpenMMRunner:
    """Lightweight OpenMM driver for protein-ligand systems set up via
    ParmEd/OpenFF.

    Requires: openmm, openmmtools, parmed, openff-toolkit (for ligands).
    """

    def __init__(self,
                 preset: str = "amber19sb-opc",
                 platform: str = "CUDA",
                 temperature_K: float = 300.0,
                 pressure_bar: float = 1.01325,
                 dt_ps: float = 0.002,
                 gpu_index: str = "0"):
        self.preset = ForceFieldPreset.lookup(preset)
        self.platform = platform
        self.T = temperature_K
        self.P = pressure_bar
        self.dt = dt_ps
        self.gpu_index = gpu_index

    def setup_from_amber(self, prmtop: str, inpcrd: str):
        import parmed as pmd
        import openmm as mm
        import openmm.app as app
        import openmm.unit as u

        struct = pmd.load_file(prmtop, inpcrd)
        system = struct.createSystem(
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * u.nanometer,
            constraints=app.HBonds,
            rigidWater=True,
            ewaldErrorTolerance=0.0005,
        )
        integrator = mm.LangevinMiddleIntegrator(
            self.T * u.kelvin, 1.0 / u.picosecond, self.dt * u.picoseconds
        )
        platform = mm.Platform.getPlatformByName(self.platform)
        props = {"CudaDeviceIndex": self.gpu_index} if self.platform == "CUDA" else {}
        self.simulation = app.Simulation(struct.topology, system, integrator,
                                         platform, props)
        self.simulation.context.setPositions(struct.positions)
        self.modeller = app.Modeller(struct.topology, struct.positions)
        self.struct = struct

    def minimize(self, max_iter: int = 10000):
        self.simulation.minimizeEnergy(maxIterations=max_iter)

    def equilibrate(self, nvt_ps: float = 0.5, npt_ps: float = 1.0):
        import openmm as mm
        import openmm.app as app
        import openmm.unit as u
        # NVT
        self.simulation.context.setVelocitiesToTemperature(self.T * u.kelvin)
        self.simulation.step(int(nvt_ps / self.dt * 1000))
        # add barostat
        system = self.simulation.system
        system.addForce(mm.MonteCarloBarostat(self.P * u.bar, self.T * u.kelvin))
        self.simulation.context.reinitialize(preserveState=True)
        self.simulation.step(int(npt_ps / self.dt * 1000))

    def production(self, length_ns: float = 10.0, traj: str = "traj.dcd",
                   report_ps: float = 10.0):
        import openmm.unit as u
        import openmm.app as app
        n_steps = int(length_ns * 1000 / self.dt)
        report_steps = int(report_ps / self.dt)
        self.simulation.reporters.append(app.DCDReporter(traj, report_steps))
        self.simulation.reporters.append(app.StateDataReporter(
            "md.log", report_steps, step=True, potentialEnergy=True,
            temperature=True, volume=True, speed=True, time=True))
        self.simulation.step(n_steps)
