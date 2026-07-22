"""
Membrane Simulator Module Tab (Work In Progress)
Roadmap and placeholder interface for Membrane Protein Simulation system construction.
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QProgressBar, QListWidget, QListWidgetItem
)
from PyQt6.QtGui import QFont, QColor
from PyQt6.QtCore import Qt


class MembraneSimulatorTab(QWidget):
    """Membrane Simulator WIP Placeholder Tab."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(15)

        # Header Title
        title_box = QGroupBox("🧬 Membrane Protein Simulator — [Work In Progress]")
        t_layout = QVBoxLayout(title_box)

        lbl = QLabel("🚧 This module is currently under active development.")
        lbl.setFont(QFont("Segoe UI", 12, QFont.Weight.Bold))
        lbl.setStyleSheet("color: #ff9800;")

        desc = QLabel(
            "The Membrane Simulator will support automated lipid bilayer insertion, "
            "orientations of proteins in membranes (OPM database integration), "
            "CHARMM-GUI / Insane lipid box generation, and multi-component lipid mixtures."
        )
        desc.setWordWrap(True)

        t_layout.addWidget(lbl)
        t_layout.addWidget(desc)

        # Planned Features List
        features_box = QGroupBox("✨ Planned Features & Roadmap")
        f_layout = QVBoxLayout(features_box)

        list_widget = QListWidget()
        features = [
            "🟢 Automated Protein Orientation in Bilayer (OPM / PPM Server)",
            "🟢 Homogeneous & Heterogeneous Bilayer Creation (POPC, DPPC, DMPC, POPG, Cholesterol)",
            "🟢 Insane / Membed Insertion Pipeline Integration",
            "🟢 CHARMM36m & AMBER LIPID21 Force Field Preset Auto-configuration",
            "🟢 Area Per Lipid (APL) & Bilayer Thickness Analysis Tools",
            "🟡 Automated Asymmetric Membrane Builder (WIP)",
        ]

        for feat in features:
            item = QListWidgetItem(feat)
            list_widget.addItem(item)

        f_layout.addWidget(list_widget)

        # Progress Status
        progress_box = QGroupBox("Development Progress")
        p_layout = QVBoxLayout(progress_box)

        p_lbl = QLabel("Module Development Completion: 45%")
        p_bar = QProgressBar()
        p_bar.setValue(45)

        p_layout.addWidget(p_lbl)
        p_layout.addWidget(p_bar)

        layout.addWidget(title_box)
        layout.addWidget(features_box)
        layout.addWidget(progress_box)
        layout.addStretch()
