"""
GROMACS Trajectory Analysis & Visualization Core
Parses XVG files and executes GROMACS analysis commands (trjconv, rms, rmspf, gyrate, energy).
"""

import os
import re
import pandas as pd
import numpy as np


def parse_xvg(xvg_file_path):
    """
    Parse a GROMACS XVG file.
    Returns a dict with metadata:
      - 'title': plot title
      - 'xaxis': x-axis label
      - 'yaxis': y-axis label
      - 'legends': list of series names
      - 'data': pandas DataFrame containing time/residue and dataset values
    """
    if not os.path.exists(xvg_file_path):
        raise FileNotFoundError(f"XVG file not found: {xvg_file_path}")

    title = os.path.basename(xvg_file_path)
    xaxis = "X"
    yaxis = "Y"
    legends = []

    raw_data_lines = []

    with open(xvg_file_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue

            if line_str.startswith("#"):
                continue

            if line_str.startswith("@"):
                # Header metadata parsing
                if 'title "' in line_str:
                    m = re.search(r'title\s+"([^"]+)"', line_str)
                    if m:
                        title = m.group(1)
                if 'xaxis' in line_str and 'label "' in line_str:
                    m = re.search(r'xaxis\s+label\s+"([^"]+)"', line_str)
                    if m:
                        xaxis = m.group(1)
                if 'yaxis' in line_str and 'label "' in line_str:
                    m = re.search(r'yaxis\s+label\s+"([^"]+)"', line_str)
                    if m:
                        yaxis = m.group(1)
                if line_str.startswith("@ s") and 'legend "' in line_str:
                    m = re.search(r'legend\s+"([^"]+)"', line_str)
                    if m:
                        legends.append(m.group(1))
                continue

            # Data numerical line
            parts = line_str.split()
            try:
                vals = [float(p) for p in parts]
                raw_data_lines.append(vals)
            except ValueError:
                continue

    if not raw_data_lines:
        df = pd.DataFrame()
    else:
        num_cols = len(raw_data_lines[0])
        col_names = [xaxis if xaxis else "X"]

        if legends and len(legends) == num_cols - 1:
            col_names.extend(legends)
        else:
            for i in range(1, num_cols):
                col_names.append(f"Series {i}")

        df = pd.DataFrame(raw_data_lines, columns=col_names[:num_cols])

    return {
        "title": title,
        "xaxis": xaxis,
        "yaxis": yaxis,
        "legends": legends,
        "df": df,
    }


def build_trjconv_cmd(gmx_path, tpr_file, xtc_file, out_file, pbc_mode="mol", center=True):
    """Build command line for trjconv PBC removal."""
    cmd = [
        gmx_path,
        "trjconv",
        "-s", tpr_file,
        "-f", xtc_file,
        "-o", out_file,
        "-pbc", pbc_mode,
    ]
    if center:
        cmd.extend(["-center"])
    return cmd


def build_analysis_cmd(gmx_path, analysis_type, tpr_file, xtc_file, out_xvg):
    """Build GROMACS analysis command for RMSD, RMSF, or Gyrate."""
    atype = analysis_type.lower()
    if atype == "rmsd":
        return [gmx_path, "rms", "-s", tpr_file, "-f", xtc_file, "-o", out_xvg]
    elif atype == "rmsf":
        return [gmx_path, "rmsf", "-s", tpr_file, "-f", xtc_file, "-o", out_xvg, "-res"]
    elif atype in ["gyrate", "rg"]:
        return [gmx_path, "gyrate", "-s", tpr_file, "-f", xtc_file, "-o", out_xvg]
    elif atype == "hbond":
        return [gmx_path, "hbond", "-s", tpr_file, "-f", xtc_file, "-num", out_xvg]
    else:
        raise ValueError(f"Unknown analysis type: {analysis_type}")
