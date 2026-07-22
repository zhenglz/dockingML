"""
Configuration and Settings Module Tab
Manages GROMACS executable paths, external software integration, themes, and system settings.
"""

import os
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QFormLayout,
    QLineEdit, QPushButton, QComboBox, QSpinBox, QLabel,
    QFileDialog, QMessageBox
)
from PyQt6.QtGui import QFont

from src.core.config import (
    load_config, save_config, detect_gmx_path,
    is_valid_gmx, get_gmx_version
)


class SettingsTab(QWidget):
    """Application Configuration Settings Tab Interface."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.config = load_config()
        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout(self)
        main_layout.setSpacing(15)

        # GROMACS Engine Settings Box
        gmx_box = QGroupBox("⚡ GROMACS Backend Configuration")
        gmx_form = QFormLayout(gmx_box)

        gmx_path_layout = QHBoxLayout()
        self.edit_gmx_path = QLineEdit(self.config.get("gmx_path", "gmx"))
        btn_browse_gmx = QPushButton("Browse...")
        btn_browse_gmx.clicked.connect(self.browse_gmx)

        btn_detect_gmx = QPushButton("🔍 Auto-Detect")
        btn_detect_gmx.clicked.connect(self.auto_detect_gmx)

        btn_test_gmx = QPushButton("🧪 Test Executable")
        btn_test_gmx.clicked.connect(self.test_gmx)

        gmx_path_layout.addWidget(self.edit_gmx_path)
        gmx_path_layout.addWidget(btn_browse_gmx)
        gmx_path_layout.addWidget(btn_detect_gmx)
        gmx_path_layout.addWidget(btn_test_gmx)

        self.lbl_gmx_status = QLabel("Status: Click 'Test Executable' to check GROMACS response.")
        self.lbl_gmx_status.setStyleSheet("color: #888; font-style: italic;")

        gmx_form.addRow("GROMACS Executable (GMX_PATH):", gmx_path_layout)
        gmx_form.addRow("GROMACS Detection Status:", self.lbl_gmx_status)

        main_layout.addWidget(gmx_box)

        # External Software Integration
        ext_box = QGroupBox("🔬 External Visualization Tools Integration")
        ext_form = QFormLayout(ext_box)

        self.edit_vmd = QLineEdit(self.config.get("vmd_path", "vmd"))
        self.edit_pymol = QLineEdit(self.config.get("pymol_path", "pymol"))
        self.edit_chimera = QLineEdit(self.config.get("chimera_path", "chimera"))

        ext_form.addRow("VMD Executable Path:", self.edit_vmd)
        ext_form.addRow("PyMOL Executable Path:", self.edit_pymol)
        ext_form.addRow("UCSF Chimera Path:", self.edit_chimera)

        main_layout.addWidget(ext_box)

        # Application Theme & Performance Box
        app_box = QGroupBox("🎨 GUI Appearance & Performance")
        app_form = QFormLayout(app_box)

        self.combo_theme = QComboBox()
        self.combo_theme.addItems(["Modern Dark", "Modern Light", "System Default"])

        saved_theme = self.config.get("theme", "dark")
        if saved_theme == "light":
            self.combo_theme.setCurrentText("Modern Light")
        elif saved_theme == "system":
            self.combo_theme.setCurrentText("System Default")
        else:
            self.combo_theme.setCurrentText("Modern Dark")

        self.spin_threads = QSpinBox()
        self.spin_threads.setRange(0, 128)
        self.spin_threads.setValue(self.config.get("cpu_threads", 0))
        self.spin_threads.setSpecialValueText("0 (Auto / All Cores)")

        app_form.addRow("GUI Color Theme:", self.combo_theme)
        app_form.addRow("GROMACS CPU Threads (-nt):", self.spin_threads)

        main_layout.addWidget(app_box)

        # Save Configuration Button
        btn_save = QPushButton("💾 Save Configuration & Apply Settings")
        btn_save.setFont(QFont("Segoe UI", 10, QFont.Weight.Bold))
        btn_save.setStyleSheet("background-color: #1976d2; color: white; padding: 10px; border-radius: 4px;")
        btn_save.clicked.connect(self.save_settings)

        main_layout.addWidget(btn_save)
        main_layout.addStretch()

    def browse_gmx(self):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select GROMACS Executable (gmx)", "", "Executables (*.exe gmx*);;All Files (*)")
        if file_path:
            self.edit_gmx_path.setText(file_path)

    def auto_detect_gmx(self):
        detected = detect_gmx_path()
        if detected:
            self.edit_gmx_path.setText(detected)
            self.lbl_gmx_status.setText(f"✅ Auto-detected valid GROMACS at: {detected}")
            self.lbl_gmx_status.setStyleSheet("color: #4caf50; font-weight: bold;")
        else:
            self.lbl_gmx_status.setText("❌ No GROMACS executable automatically detected in system PATH.")
            self.lbl_gmx_status.setStyleSheet("color: #f44336; font-weight: bold;")

    def test_gmx(self):
        path = self.edit_gmx_path.text().strip()
        if is_valid_gmx(path):
            ver = get_gmx_version(path)
            self.lbl_gmx_status.setText(f"✅ GROMACS Verified: {ver}")
            self.lbl_gmx_status.setStyleSheet("color: #4caf50; font-weight: bold;")
            QMessageBox.information(self, "GROMACS Found", f"Successfully communicated with GROMACS!\n\n{ver}")
        else:
            self.lbl_gmx_status.setText(f"❌ Failed to run GROMACS at: {path}")
            self.lbl_gmx_status.setStyleSheet("color: #f44336; font-weight: bold;")
            QMessageBox.critical(
                self,
                "GROMACS Not Found",
                f"Could not execute GROMACS at '{path}'.\n\nPlease ensure GROMACS is installed and specify the correct path to gmx / gmx.exe in config.py or settings.",
            )

    def save_settings(self):
        path = self.edit_gmx_path.text().strip()
        theme_str = self.combo_theme.currentText()
        if "Light" in theme_str:
            t = "light"
        elif "System" in theme_str:
            t = "system"
        else:
            t = "dark"

        self.config["gmx_path"] = path
        self.config["vmd_path"] = self.edit_vmd.text().strip()
        self.config["pymol_path"] = self.edit_pymol.text().strip()
        self.config["chimera_path"] = self.edit_chimera.text().strip()
        self.config["theme"] = t
        self.config["cpu_threads"] = self.spin_threads.value()

        if save_config(self.config):
            QMessageBox.information(self, "Settings Saved", "Configuration saved successfully!")
