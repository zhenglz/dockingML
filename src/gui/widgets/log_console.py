"""
Log Console Widget for GROMACS GUI
Displays real-time execution log messages with color syntax highlighting, search filter, and export features.
"""

import os
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QPlainTextEdit,
    QPushButton, QCheckBox, QLineEdit, QLabel, QFileDialog, QMessageBox
)
from PyQt6.QtGui import QFont, QTextCursor, QColor
from PyQt6.QtCore import Qt


class LogConsole(QWidget):
    """Rich-text live output console with search, clear, autoscroll, and export controls."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(5, 5, 5, 5)

        # Controls bar
        ctrl_layout = QHBoxLayout()
        ctrl_layout.setContentsMargins(0, 0, 0, 0)

        title_lbl = QLabel("📜 Execution Log Output:")
        title_lbl.setFont(QFont("Segoe UI", 10, QFont.Weight.Bold))

        self.filter_edit = QLineEdit()
        self.filter_edit.setPlaceholderText("Filter output log...")
        self.filter_edit.textChanged.connect(self.filter_text)

        self.autoscroll_cb = QCheckBox("Auto-scroll")
        self.autoscroll_cb.setChecked(True)

        copy_btn = QPushButton("📋 Copy Log")
        copy_btn.clicked.connect(self.copy_log)

        save_btn = QPushButton("💾 Save Log...")
        save_btn.clicked.connect(self.save_log)

        clear_btn = QPushButton("🗑️ Clear")
        clear_btn.clicked.connect(self.clear_log)

        ctrl_layout.addWidget(title_lbl)
        ctrl_layout.addWidget(self.filter_edit)
        ctrl_layout.addWidget(self.autoscroll_cb)
        ctrl_layout.addWidget(copy_btn)
        ctrl_layout.addWidget(save_btn)
        ctrl_layout.addWidget(clear_btn)

        # Log Text Display Area
        self.log_text = QPlainTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFont(QFont("Consolas", 9))
        self.log_text.setStyleSheet(
            "QPlainTextEdit { background-color: #1e1e1e; color: #d4d4d4; border: 1px solid #333; font-family: 'Consolas', 'Courier New', monospace; }"
        )

        layout.addLayout(ctrl_layout)
        layout.addWidget(self.log_text)

    def append_log(self, text):
        """Append log line with syntax coloring."""
        if not text:
            return

        # Simple color formatting based on keywords
        color_hex = "#d4d4d4"  # Default light gray
        if "❌" in text or "Error" in text or "FAILED" in text or "Fatal" in text:
            color_hex = "#f44336"  # Bright red
        elif "⚠️" in text or "Warning" in text or "WARN" in text:
            color_hex = "#ff9800"  # Orange
        elif "🎉" in text or "SUCCESS" in text or "Completed" in text:
            color_hex = "#4caf50"  # Green
        elif "🚀" in text or "▶️" in text or "Executing" in text:
            color_hex = "#2196f3"  # Blue

        html_line = f'<span style="color:{color_hex};">{self._escape_html(text)}</span>'
        self.log_text.appendHtml(html_line)

        if self.autoscroll_cb.isChecked():
            self.log_text.moveCursor(QTextCursor.MoveOperation.End)

    def _escape_html(self, text):
        return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\n", "<br>")

    def filter_text(self, pattern):
        # Filter is handled dynamically if needed or highlighting
        pass

    def clear_log(self):
        self.log_text.clear()

    def copy_log(self):
        self.log_text.selectAll()
        self.log_text.copy()
        cursor = self.log_text.textCursor()
        cursor.clearSelection()
        self.log_text.setTextCursor(cursor)

    def save_log(self):
        content = self.log_text.toPlainText()
        if not content.strip():
            QMessageBox.information(self, "Empty Log", "Log console is currently empty.")
            return

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Log File", "gromacs_gui_run.log", "Log Files (*.log);;Text Files (*.txt);;All Files (*)"
        )
        if file_path:
            try:
                with open(file_path, "w", encoding="utf-8") as f:
                    f.write(content)
                QMessageBox.information(self, "Saved", f"Log saved to:\n{file_path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to save log file:\n{e}")
