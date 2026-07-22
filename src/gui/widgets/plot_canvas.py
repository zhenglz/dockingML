"""
Matplotlib Plotting Canvas Widget for GROMACS GUI
Provides interactive chart rendering, legend/axis configuration, and PNG/SVG/CSV export.
"""

import os
import pandas as pd
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPushButton,
    QLabel, QFileDialog, QMessageBox, QCheckBox, QComboBox
)
from PyQt6.QtGui import QFont

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure


class PlotCanvas(QWidget):
    """Interactive Matplotlib Canvas widget for rendering GROMACS trajectory data."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_data = None
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(5, 5, 5, 5)

        # Matplotlib Figure & Canvas
        self.figure = Figure(figsize=(7, 4), dpi=100)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)

        # Custom Controls Toolbar
        ctrl_layout = QHBoxLayout()

        self.grid_cb = QCheckBox("Show Grid")
        self.grid_cb.setChecked(True)
        self.grid_cb.stateChanged.connect(self.replot)

        self.style_combo = QComboBox()
        self.style_combo.addItems(["Default Theme", "Dark Theme", "Seaborn Light"])
        self.style_combo.currentIndexChanged.connect(self.replot)

        export_csv_btn = QPushButton("📊 Export Data (CSV)")
        export_csv_btn.clicked.connect(self.export_csv)

        export_img_btn = QPushButton("🖼️ Export Figure...")
        export_img_btn.clicked.connect(self.export_image)

        ctrl_layout.addWidget(QLabel("Plot Style:"))
        ctrl_layout.addWidget(self.style_combo)
        ctrl_layout.addWidget(self.grid_cb)
        ctrl_layout.addStretch()
        ctrl_layout.addWidget(export_csv_btn)
        ctrl_layout.addWidget(export_img_btn)

        layout.addWidget(self.toolbar)
        layout.addLayout(ctrl_layout)
        layout.addWidget(self.canvas)

    def plot_data(self, data_dict):
        """
        Plot parsed XVG data dictionary.
        data_dict = {'title': str, 'xaxis': str, 'yaxis': str, 'legends': list, 'df': DataFrame}
        """
        self.current_data = data_dict
        self.replot()

    def replot(self):
        if not self.current_data or "df" not in self.current_data:
            return

        df = self.current_data["df"]
        if df.empty or len(df.columns) < 2:
            return

        self.figure.clear()
        ax = self.figure.add_subplot(111)

        # Style customization
        style = self.style_combo.currentText()
        if style == "Dark Theme":
            self.figure.patch.set_facecolor("#2b2b2b")
            ax.set_facecolor("#1e1e1e")
            ax.tick_params(colors="white")
            ax.xaxis.label.set_color("white")
            ax.yaxis.label.set_color("white")
            ax.title.set_color("white")
            line_colors = ["#4fc3f7", "#ffb74d", "#81c784", "#e57373", "#ba68c8"]
        else:
            self.figure.patch.set_facecolor("white")
            ax.set_facecolor("#f9f9f9")
            ax.tick_params(colors="black")
            line_colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

        x_col = df.columns[0]
        y_cols = df.columns[1:]

        for i, col in enumerate(y_cols):
            c = line_colors[i % len(line_colors)]
            ax.plot(df[x_col], df[col], label=col, color=c, linewidth=1.5)

        ax.set_xlabel(self.current_data.get("xaxis", x_col), fontsize=10)
        ax.set_ylabel(self.current_data.get("yaxis", "Value"), fontsize=10)
        ax.set_title(self.current_data.get("title", "Trajectory Property"), fontsize=11, fontweight="bold")

        if self.grid_cb.isChecked():
            grid_color = "#444" if style == "Dark Theme" else "#ccc"
            ax.grid(True, linestyle="--", alpha=0.5, color=grid_color)

        if len(y_cols) > 1 or (len(y_cols) == 1 and y_cols[0] not in ["Series 1", "Y"]):
            ax.legend(loc="best", framealpha=0.8)

        self.figure.tight_layout()
        self.canvas.draw()

    def export_csv(self):
        if not self.current_data or "df" not in self.current_data or self.current_data["df"].empty:
            QMessageBox.warning(self, "No Data", "There is no data currently plotted to export.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save CSV Data", "gromacs_plot_data.csv", "CSV Files (*.csv);;All Files (*)"
        )
        if file_path:
            try:
                self.current_data["df"].to_csv(file_path, index=False)
                QMessageBox.information(self, "Exported", f"Data exported successfully to:\n{file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to export CSV:\n{e}")

    def export_image(self):
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Figure", "gromacs_plot.png", "PNG Image (*.png);;SVG Vector (*.svg);;PDF Document (*.pdf)"
        )
        if file_path:
            try:
                self.figure.savefig(file_path, dpi=300, bbox_inches="tight")
                QMessageBox.information(self, "Saved", f"Figure saved successfully to:\n{file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save figure image:\n{e}")
