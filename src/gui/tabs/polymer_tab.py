"""
Polymer Simulator Module Tab (Work In Progress)
Roadmap and placeholder interface for Polymer System Construction and MD Simulation.
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QProgressBar, QListWidget, QListWidgetItem
)
from PyQt6.QtGui import QFont


class PolymerSimulatorTab(QWidget):
    """Polymer Simulator WIP Placeholder Tab."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(15)

        # Header Title
        title_box = QGroupBox("🧱 Polymer Simulator — [Work In Progress]")
        t_layout = QVBoxLayout(title_box)

        lbl = QLabel("🚧 This module is currently under active development.")
        lbl.setFont(QFont("Segoe UI", 12, QFont.Weight.Bold))
        lbl.setStyleSheet("color: #ff9800;")

        desc = QLabel(
            "The Polymer Simulator will support automated polymer chain generation, "
            "monomer topology parameterization, amorphous polymer cell packings, "
            "and glass transition temperature (Tg) / mechanical property analysis."
        )
        desc.setWordWrap(True)

        t_layout.addWidget(lbl)
        t_layout.addWidget(desc)

        # Planned Features
        features_box = QGroupBox("✨ Planned Features & Roadmap")
        f_layout = QVBoxLayout(features_box)

        list_widget = QListWidget()
        features = [
            "🟢 Automated Monomer & Oligomer Assembly (Polymer Builder)",
            "🟢 OPLS-AA & GAFF2 Parameterization for Synthetic Polymers",
            "🟢 Amorphous Cell Density Equalization & Relaxation Workflow",
            "🟢 Glass Transition Temperature (Tg) Automated Annealing Pipeline",
            "🟢 Young's Modulus & Tensile Stress Calculation Helpers",
        ]

        for feat in features:
            list_widget.addItem(QListWidgetItem(feat))

        f_layout.addWidget(list_widget)

        # Progress Status
        progress_box = QGroupBox("Development Progress")
        p_layout = QVBoxLayout(progress_box)

        p_lbl = QLabel("Module Development Completion: 35%")
        p_bar = QProgressBar()
        p_bar.setValue(35)

        p_layout.addWidget(p_lbl)
        p_layout.addWidget(p_bar)

        layout.addWidget(title_box)
        layout.addWidget(features_box)
        layout.addWidget(progress_box)
        layout.addStretch()
