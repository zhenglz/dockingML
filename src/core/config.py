"""
GROMACS GUI Configuration Module
Manages application configuration, GROMACS paths, forcefields, and MDP default templates.
"""

import os
import json
import shutil
import subprocess

# Default GROMACS path - can be updated by user or configuration file
GMX_PATH = r"gmx"

CONFIG_FILE = os.path.expanduser("~/.gromacs_gui.json")

# Forcefields options for pdb2gmx
FORCEFIELDS = [
    {"id": "amber19sb-opc", "name": "AMBER ff19SB with OPC water (Recommended 2024+)", "ff": "amber19sb"},
    {"id": "amber99sb-ildn", "name": "AMBER99SB-ILDN protein, nucleic AMBER94", "ff": "amber99sb-ildn"},
    {"id": "amber14sb", "name": "AMBER14SB protein, nucleic AMBER94", "ff": "amber14sb"},
    {"id": "charmm36m", "name": "CHARMM36 all-atom force field (jul 2021)", "ff": "charmm36-jul2021"},
    {"id": "opls-aa", "name": "OPLS-AA/L all-atom force field (2001 aminoacid) ", "ff": "oplsaa"},
    {"id": "gromos54a7", "name": "GROMOS54a7 force field", "ff": "gromos54a7"},
]

# Water models
WATER_MODELS = [
    {"id": "opc", "name": "OPC (Modern 4-point, optimal for AMBER ff19SB)"},
    {"id": "tip3p", "name": "TIP3P (Standard 3-point water)"},
    {"id": "spce", "name": "SPC/E (Extended Simple Point Charge)"},
    {"id": "tip4p", "name": "TIP4P (4-point water)"},
    {"id": "none", "name": "None (For gas-phase or non-solvated)"},
]

# Box shapes
BOX_SHAPES = [
    {"id": "cubic", "name": "Cubic"},
    {"id": "dodecahedron", "name": "Rhombic Dodecahedron (Optimal volume)"},
    {"id": "triclinic", "name": "Triclinic"},
    {"id": "octahedron", "name": "Truncated Octahedron"},
]

DEFAULT_CONFIG = {
    "gmx_path": GMX_PATH,
    "vmd_path": "vmd",
    "pymol_path": "pymol",
    "chimera_path": "chimera",
    "cpu_threads": 0,
    "theme": "dark",
    "last_workspace": "",
}


def load_config():
    """Load configuration from JSON config file or auto-detect settings."""
    config = DEFAULT_CONFIG.copy()
    if os.path.exists(CONFIG_FILE):
        try:
            with open(CONFIG_FILE, "r", encoding="utf-8") as f:
                saved = json.load(f)
                config.update(saved)
        except Exception as e:
            print(f"[Warning] Failed to load config from {CONFIG_FILE}: {e}")

    # Auto-detect GMX path if default isn't valid
    if not is_valid_gmx(config["gmx_path"]):
        detected = detect_gmx_path()
        if detected:
            config["gmx_path"] = detected

    global GMX_PATH
    GMX_PATH = config["gmx_path"]
    return config


def save_config(config_data):
    """Save configuration dictionary to file."""
    try:
        with open(CONFIG_FILE, "w", encoding="utf-8") as f:
            json.dump(config_data, f, indent=4)
        global GMX_PATH
        if "gmx_path" in config_data:
            GMX_PATH = config_data["gmx_path"]
        return True
    except Exception as e:
        print(f"[Error] Failed to save config to {CONFIG_FILE}: {e}")
        return False


def detect_gmx_path():
    """Try to detect local GROMACS executable path."""
    for cmd in ["gmx", "gmx_mpi", r"C:\Program Files\gromacs\bin\gmx.exe"]:
        path = shutil.which(cmd)
        if path and is_valid_gmx(path):
            return path
        if os.path.exists(cmd) and is_valid_gmx(cmd):
            return cmd
    return ""


def is_valid_gmx(gmx_executable):
    """Test if GROMACS executable exists and responds to version call."""
    if not gmx_executable:
        return False
    try:
        res = subprocess.run(
            [gmx_executable, "version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=5,
        )
        return res.returncode == 0 and "GROMACS" in (res.stdout + res.stderr)
    except Exception:
        return False


def get_gmx_version(gmx_executable=None):
    """Get string version info from GROMACS."""
    path = gmx_executable or GMX_PATH
    if not is_valid_gmx(path):
        return "Not Found / Invalid"
    try:
        res = subprocess.run(
            [path, "version"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=5,
        )
        for line in (res.stdout + res.stderr).splitlines():
            if "GROMACS version" in line:
                return line.strip()
        return "GROMACS detected"
    except Exception as e:
        return f"Error: {e}"
