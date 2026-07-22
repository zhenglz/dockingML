"""
Automated Test Suite for GROMACS GUI Application
Verifies backend modules, GUI widgets, dialogs, tabs, main window, and headless execution.
"""

import os
import sys

# Add project root to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import unittest
import tempfile
import pandas as pd

# Headless env bootstrap
from src.core import env_setup
env_setup.prepare_headless_env()

from PyQt6.QtWidgets import QApplication

from src.core.config import (
    load_config, save_config, is_valid_gmx,
    detect_gmx_path, GMX_PATH, FORCEFIELDS, WATER_MODELS
)
from src.core.mdp_editor import MDPEditor
from src.core.topology_builder import merge_gro_files, update_topology
from src.core.analysis import parse_xvg, build_trjconv_cmd, build_analysis_cmd
from src.core.gmx_runner import GMXWorker

from src.gui.widgets.log_console import LogConsole
from src.gui.widgets.plot_canvas import PlotCanvas
from src.gui.mdp_dialog import MDPDialog
from src.gui.tabs.ligand_tab import LigandSimulatorTab
from src.gui.tabs.solution_tab import SolutionSimulatorTab
from src.gui.tabs.membrane_tab import MembraneSimulatorTab
from src.gui.tabs.polymer_tab import PolymerSimulatorTab
from src.gui.tabs.analysis_tab import AnalysisTab
from src.gui.tabs.settings_tab import SettingsTab
from src.gui.main_window import MainWindow

app = None


def setUpModule():
    global app
    app = QApplication.instance()
    if app is None:
        app = QApplication([])


class TestBackendCore(unittest.TestCase):

    def test_config(self):
        cfg = load_config()
        self.assertIn("gmx_path", cfg)
        self.assertTrue(len(FORCEFIELDS) > 0)
        self.assertTrue(len(WATER_MODELS) > 0)

    def test_mdp_editor(self):
        preset_em = MDPEditor.get_preset("em")
        self.assertEqual(preset_em.get("integrator"), "steep")

        preset_nvt = MDPEditor.get_preset("nvt")
        self.assertEqual(preset_nvt.get("integrator"), "md")

        preset_npt = MDPEditor.get_preset("npt")
        self.assertIn("Pcoupl", preset_npt)

        preset_md = MDPEditor.get_preset("md")
        self.assertEqual(preset_md.get("integrator"), "md")

        with tempfile.TemporaryDirectory() as tmpdir:
            mdp_path = os.path.join(tmpdir, "test.mdp")
            MDPEditor.save_mdp(preset_em, mdp_path, title="Test EM MDP")
            self.assertTrue(os.path.exists(mdp_path))

            parsed = MDPEditor.parse_mdp(mdp_path)
            self.assertEqual(parsed.get("integrator"), "steep")

    def test_topology_builder(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            rec_gro = os.path.join(tmpdir, "rec.gro")
            lig_gro = os.path.join(tmpdir, "lig.gro")
            cpx_gro = os.path.join(tmpdir, "complex.gro")
            topol_top = os.path.join(tmpdir, "topol.top")
            lig_itp = os.path.join(tmpdir, "ligand.itp")

            # Create dummy rec.gro
            with open(rec_gro, "w") as f:
                f.write("Receptor\n 2\n    1ALA     CA    1   0.100   0.200   0.300\n    1ALA     CB    2   0.150   0.250   0.350\n   2.00000   2.00000   2.00000\n")

            # Create dummy lig.gro
            with open(lig_gro, "w") as f:
                f.write("Ligand\n 1\n    1LIG    C01    1   0.500   0.500   0.500\n   1.00000   1.00000   1.00000\n")

            # Create dummy topol.top
            with open(topol_top, "w") as f:
                f.write('; Topology\n#include "amber99sb-ildn.ff/forcefield.itp"\n[ system ]\nTest System\n[ molecules ]\nProtein_chain_A 1\n')

            # Create dummy ligand.itp
            with open(lig_itp, "w") as f:
                f.write("; Ligand itp\n[ moleculetype ]\nLIG 3\n")

            # Test merging GRO
            merge_gro_files(rec_gro, lig_gro, cpx_gro)
            self.assertTrue(os.path.exists(cpx_gro))
            with open(cpx_gro, "r") as f:
                content = f.read()
                self.assertIn("3", content)  # Total 2 + 1 = 3 atoms
                self.assertIn("LIG", content)

            # Test updating topol.top
            update_topology(topol_top, lig_itp, ligand_resname="LIG", count=1)
            with open(topol_top, "r") as f:
                top_content = f.read()
                self.assertIn('#include "ligand.itp"', top_content)
                self.assertIn("LIG", top_content)

    def test_analysis(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            xvg_path = os.path.join(tmpdir, "test.xvg")
            with open(xvg_path, "w") as f:
                f.write('# Comments\n@ title "RMSD Analysis"\n@ xaxis label "Time (ps)"\n@ yaxis label "RMSD (nm)"\n0.000 0.120\n10.000 0.180\n20.000 0.210\n')

            data = parse_xvg(xvg_path)
            self.assertEqual(data["title"], "RMSD Analysis")
            self.assertEqual(data["xaxis"], "Time (ps)")
            self.assertEqual(data["yaxis"], "RMSD (nm)")
            self.assertFalse(data["df"].empty)
            self.assertEqual(len(data["df"]), 3)


class TestGUIComponents(unittest.TestCase):

    def test_log_console(self):
        console = LogConsole()
        console.append_log("🚀 Launching Test")
        console.append_log("❌ Error Test")
        self.assertIn("Error Test", console.log_text.toPlainText())

    def test_plot_canvas(self):
        canvas = PlotCanvas()
        df = pd.DataFrame({"Time (ps)": [0, 10, 20], "RMSD": [0.1, 0.2, 0.25]})
        canvas.plot_data({"title": "Test Plot", "xaxis": "Time (ps)", "yaxis": "RMSD (nm)", "df": df})
        self.assertIsNotNone(canvas.current_data)

    def test_mdp_dialog(self):
        dialog = MDPDialog(stage="em")
        self.assertEqual(dialog.stage, "em")

    def test_tabs(self):
        lig_tab = LigandSimulatorTab()
        sol_tab = SolutionSimulatorTab()
        mem_tab = MembraneSimulatorTab()
        pol_tab = PolymerSimulatorTab()
        ana_tab = AnalysisTab()
        set_tab = SettingsTab()

        self.assertIsNotNone(lig_tab)
        self.assertIsNotNone(sol_tab)
        self.assertIsNotNone(mem_tab)
        self.assertIsNotNone(pol_tab)
        self.assertIsNotNone(ana_tab)
        self.assertIsNotNone(set_tab)

    def test_main_window(self):
        win = MainWindow()
        win.show()
        self.assertTrue(win.isVisible())
        self.assertEqual(win.tabs.count(), 6)


if __name__ == "__main__":
    unittest.main()
