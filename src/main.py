"""
GROMACS Graphical User Interface (GROMACS-GUI) Entry Point
Launches the PyQt6 application with proper headless environment setup and QSS theme stylesheet.
"""

import os
import sys

# Ensure repository root is on sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

# Ensure headless environment dynamic libraries are bootstrapped if needed
from src.core import env_setup
env_setup.prepare_headless_env()

from PyQt6.QtWidgets import QApplication
from PyQt6.QtGui import QFont
from src.core.config import load_config
from src.gui.main_window import MainWindow

DARK_STYLESHEET = """
QMainWindow, QDialog, QWidget {
    background-color: #2b2b2b;
    color: #e0e0e0;
    font-family: 'Segoe UI', Arial, sans-serif;
}

QGroupBox {
    border: 1px solid #444;
    border-radius: 6px;
    margin-top: 10px;
    padding-top: 10px;
    font-weight: bold;
    color: #1976d2;
}

QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    padding: 0 5px;
}

QPushButton {
    background-color: #3c3f41;
    color: #ffffff;
    border: 1px solid #555;
    border-radius: 4px;
    padding: 6px 12px;
    min-height: 22px;
}

QPushButton:hover {
    background-color: #484b4d;
    border-color: #1976d2;
}

QPushButton:pressed {
    background-color: #1976d2;
}

QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {
    background-color: #1e1e1e;
    color: #ffffff;
    border: 1px solid #555;
    border-radius: 4px;
    padding: 4px 8px;
    min-height: 22px;
}

QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {
    border: 1px solid #1976d2;
}

QTabWidget::pane {
    border: 1px solid #444;
    border-radius: 4px;
    background-color: #2b2b2b;
}

QTabBar::tab {
    background-color: #222;
    color: #aaa;
    padding: 8px 16px;
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
}

QTabBar::tab:selected {
    background-color: #2b2b2b;
    color: #4fc3f7;
    font-weight: bold;
    border-top: 2px solid #1976d2;
}

QProgressBar {
    border: 1px solid #444;
    border-radius: 4px;
    text-align: center;
    background-color: #1e1e1e;
}

QProgressBar::chunk {
    background-color: #1976d2;
    width: 10px;
}

QStatusBar {
    background-color: #1e1e1e;
    color: #aaa;
}
"""


def main():
    config = load_config()

    app = QApplication(sys.argv)
    app.setApplicationName("GROMACS-GUI")
    app.setStyle("Fusion")

    # Apply dark theme if requested
    theme = config.get("theme", "dark")
    if theme == "dark":
        app.setStyleSheet(DARK_STYLESHEET)

    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
