"""
Main Application Window for GROMACS GUI
Assembles tabs, header title bar, status bar, menu bar, and theme styling.
"""

import os
from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QTabWidget,
    QLabel, QStatusBar, QMenuBar, QMenu, QMessageBox, QFileDialog
)
from PyQt6.QtGui import QFont, QIcon, QAction
from PyQt6.QtCore import Qt

from src.core.config import load_config, get_gmx_version, is_valid_gmx
from src.gui.tabs.ligand_tab import LigandSimulatorTab
from src.gui.tabs.solution_tab import SolutionSimulatorTab
from src.gui.tabs.membrane_tab import MembraneSimulatorTab
from src.gui.tabs.polymer_tab import PolymerSimulatorTab
from src.gui.tabs.analysis_tab import AnalysisTab
from src.gui.tabs.settings_tab import SettingsTab


class MainWindow(QMainWindow):
    """GROMACS GUI Primary Main Window Application UI Container."""

    def __init__(self):
        super().__init__()
        self.config = load_config()
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle("GROMACS Graphical User Interface (GROMACS-GUI)")
        self.resize(1150, 800)

        # Central Widget & Main Layout
        central_widget = QWidget()
        main_layout = QVBoxLayout(central_widget)
        main_layout.setContentsMargins(10, 10, 10, 10)

        # Top Title Header Banner
        header_layout = QHBoxLayout()

        title_lbl = QLabel("🧪 GROMACS-GUI Workbench")
        title_lbl.setFont(QFont("Segoe UI", 16, QFont.Weight.Bold))
        title_lbl.setStyleSheet("color: #1976d2;")

        subtitle_lbl = QLabel("Graphical Configuration & Task Manager for Molecular Dynamics")
        subtitle_lbl.setFont(QFont("Segoe UI", 9))
        subtitle_lbl.setStyleSheet("color: #777;")

        header_info_layout = QVBoxLayout()
        header_info_layout.addWidget(title_lbl)
        header_info_layout.addWidget(subtitle_lbl)

        header_layout.addLayout(header_info_layout)
        header_layout.addStretch()

        main_layout.addLayout(header_layout)

        # Main Tab Widget
        self.tabs = QTabWidget()
        self.tabs.setFont(QFont("Segoe UI", 10))

        # Add Tab Modules
        self.ligand_tab = LigandSimulatorTab()
        self.solution_tab = SolutionSimulatorTab()
        self.membrane_tab = MembraneSimulatorTab()
        self.polymer_tab = PolymerSimulatorTab()
        self.analysis_tab = AnalysisTab()
        self.settings_tab = SettingsTab()

        self.tabs.addTab(self.ligand_tab, "🧬 Ligand Simulator")
        self.tabs.addTab(self.solution_tab, "💧 Solution Simulator")
        self.tabs.addTab(self.membrane_tab, "🧱 Membrane Simulator (WIP)")
        self.tabs.addTab(self.polymer_tab, "🧪 Polymer Simulator (WIP)")
        self.tabs.addTab(self.analysis_tab, "📊 Analysis & Visualization")
        self.tabs.addTab(self.settings_tab, "⚙️ Settings")

        main_layout.addWidget(self.tabs)

        self.setCentralWidget(central_widget)

        # Build Menu Bar
        self._create_menu_bar()

        # Build Status Bar
        self._create_status_bar()

    def _create_menu_bar(self):
        menu_bar = self.menuBar()

        # File Menu
        file_menu = menu_bar.addMenu("📁 File")

        action_open_ws = QAction("Open Workspace Directory...", self)
        action_open_ws.triggered.connect(self.open_workspace)
        file_menu.addAction(action_open_ws)

        file_menu.addSeparator()

        action_exit = QAction("Exit App", self)
        action_exit.triggered.connect(self.close)
        file_menu.addAction(action_exit)

        # Help Menu
        help_menu = menu_bar.addMenu("❓ Help")

        action_about = QAction("About GROMACS-GUI...", self)
        action_about.triggered.connect(self.show_about_dialog)
        help_menu.addAction(action_about)

    def _create_status_bar(self):
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)

        gmx_path = self.config.get("gmx_path", "gmx")
        if is_valid_gmx(gmx_path):
            gmx_info = get_gmx_version(gmx_path)
            self.status_bar.showMessage(f"GROMACS Status: Ready ({gmx_info})")
        else:
            self.status_bar.showMessage("GROMACS Status: ⚠️ Executable not found! Please set path in Settings.")

    def open_workspace(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Workspace Directory", os.getcwd())
        if dir_path:
            self.ligand_tab.work_dir = dir_path
            self.ligand_tab.lbl_workdir.setText(dir_path)

            self.solution_tab.work_dir = dir_path
            self.solution_tab.lbl_workdir.setText(dir_path)

            self.analysis_tab.work_dir = dir_path
            self.analysis_tab.lbl_workdir.setText(dir_path)

            self.status_bar.showMessage(f"Active Workspace: {dir_path}")

    def show_about_dialog(self):
        QMessageBox.about(
            self,
            "About GROMACS-GUI",
            "<h3>GROMACS Graphical User Interface (GROMACS-GUI)</h3>"
            "<p>A modern Python & PyQt6 Graphical User Interface tool for GROMACS Molecular Dynamics simulations.</p>"
            "<p><b>Features:</b></p>"
            "<ul>"
            "<li><b>Ligand Simulator:</b> Complete Protein-Ligand complex preparation & pipeline</li>"
            "<li><b>Solution Simulator:</b> Solution system topology, box setup & MD</li>"
            "<li><b>Interactive MDP Parameter Editor:</b> Graphical parameters configuration</li>"
            "<li><b>Trajectory Analysis:</b> PBC removal (trjconv) & Matplotlib property plotting</li>"
            "<li><b>External Tools Integration:</b> VMD, PyMOL, Chimera</li>"
            "</ul>"
            "<p><i>Built with Python 3, PyQt6, and Matplotlib.</i></p>"
        )
