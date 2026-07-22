"""
Analysis & Visualization Module Tab
Provides trajectory PBC processing (trjconv), property calculation (RMSD, RMSF, Gyrate, Energy),
embedded Matplotlib charting, and external visualization launcher (VMD, PyMOL, Chimera).
"""

import os
import subprocess
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QFormLayout,
    QLineEdit, QPushButton, QComboBox, QCheckBox, QLabel,
    QFileDialog, QMessageBox, QTabWidget, QSplitter
)
from PyQt6.QtGui import QFont
from PyQt6.QtCore import Qt

from src.core.config import GMX_PATH, load_config
from src.core.analysis import parse_xvg, build_trjconv_cmd, build_analysis_cmd
from src.core.gmx_runner import GMXWorker
from src.gui.widgets.plot_canvas import PlotCanvas
from src.gui.widgets.log_console import LogConsole


class AnalysisTab(QWidget):
    """Trajectory Analysis and Visualization Tab interface."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.work_dir = os.getcwd()
        self.active_worker = None
        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout(self)

        # Working directory header
        workdir_layout = QHBoxLayout()
        self.lbl_workdir = QLineEdit(self.work_dir)
        btn_browse_workdir = QPushButton("📂 Set Working Directory...")
        btn_browse_workdir.clicked.connect(self.browse_working_dir)
        workdir_layout.addWidget(QLabel("📁 Working Directory:"))
        workdir_layout.addWidget(self.lbl_workdir)
        workdir_layout.addWidget(btn_browse_workdir)
        main_layout.addLayout(workdir_layout)

        # Create Horizontal Splitter (Controls Left, Chart Right)
        splitter = QSplitter(Qt.Orientation.Horizontal)

        # Left Controls Widget
        controls_widget = QWidget()
        ctrl_layout = QVBoxLayout(controls_widget)

        # Sub-tabs for Analysis Features
        sub_tabs = QTabWidget()

        # Sub-tab 1: Trajectory Processing (trjconv)
        sub_tabs.addTab(self._create_trjconv_tab(), "🔄 PBC Removal")

        # Sub-tab 2: Calculation Analysis (RMSD/RMSF/Gyrate)
        sub_tabs.addTab(self._create_calc_tab(), "📈 Calculate Properties")

        # Sub-tab 3: XVG File Plotter
        sub_tabs.addTab(self._create_xvg_tab(), "📄 Plot XVG File")

        # Sub-tab 4: External Visualizers
        sub_tabs.addTab(self._create_external_tab(), "🚀 External Tools")

        ctrl_layout.addWidget(sub_tabs)

        # Log Console at bottom left
        self.console = LogConsole()
        ctrl_layout.addWidget(self.console)

        # Right Chart Canvas Widget
        self.canvas = PlotCanvas()

        splitter.addWidget(controls_widget)
        splitter.addWidget(self.canvas)
        splitter.setSizes([450, 650])

        main_layout.addWidget(splitter)

    def browse_working_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Working Directory", self.work_dir)
        if dir_path:
            self.work_dir = dir_path
            self.lbl_workdir.setText(dir_path)

    def browse_file(self, target_line_edit, file_filter):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", self.work_dir, file_filter)
        if file_path:
            target_line_edit.setText(file_path)

    def _create_trjconv_tab(self):
        widget = QWidget()
        form = QFormLayout(widget)

        tpr_layout = QHBoxLayout()
        self.edit_trj_tpr = QLineEdit("md.tpr")
        btn_tpr = QPushButton("Browse...")
        btn_tpr.clicked.connect(lambda: self.browse_file(self.edit_trj_tpr, "TPR Files (*.tpr);;All Files (*)"))
        tpr_layout.addWidget(self.edit_trj_tpr)
        tpr_layout.addWidget(btn_tpr)

        xtc_layout = QHBoxLayout()
        self.edit_trj_xtc = QLineEdit("md.xtc")
        btn_xtc = QPushButton("Browse...")
        btn_xtc.clicked.connect(lambda: self.browse_file(self.edit_trj_xtc, "Trajectory Files (*.xtc *.trr);;All Files (*)"))
        xtc_layout.addWidget(self.edit_trj_xtc)
        xtc_layout.addWidget(btn_xtc)

        self.edit_trj_out = QLineEdit("md_noPBC.xtc")

        self.combo_pbc = QComboBox()
        self.combo_pbc.addItems(["mol", "res", "atom", "nojump", "cluster"])

        self.cb_center = QCheckBox("Center System (-center)")
        self.cb_center.setChecked(True)

        btn_run_trjconv = QPushButton("▶️ Run PBC Removal (trjconv)")
        btn_run_trjconv.setFont(QFont("Segoe UI", 9, QFont.Weight.Bold))
        btn_run_trjconv.clicked.connect(self.run_trjconv)

        form.addRow("Structure/Topology (.tpr):", tpr_layout)
        form.addRow("Input Trajectory (.xtc):", xtc_layout)
        form.addRow("Output Trajectory Name:", self.edit_trj_out)
        form.addRow("PBC Mode (-pbc):", self.combo_pbc)
        form.addRow(self.cb_center)
        form.addRow(btn_run_trjconv)

        return widget

    def _create_calc_tab(self):
        widget = QWidget()
        form = QFormLayout(widget)

        tpr_layout = QHBoxLayout()
        self.edit_calc_tpr = QLineEdit("md.tpr")
        btn_tpr = QPushButton("Browse...")
        btn_tpr.clicked.connect(lambda: self.browse_file(self.edit_calc_tpr, "TPR Files (*.tpr);;All Files (*)"))
        tpr_layout.addWidget(self.edit_calc_tpr)
        tpr_layout.addWidget(btn_tpr)

        xtc_layout = QHBoxLayout()
        self.edit_calc_xtc = QLineEdit("md_noPBC.xtc")
        btn_xtc = QPushButton("Browse...")
        btn_xtc.clicked.connect(lambda: self.browse_file(self.edit_calc_xtc, "Trajectory Files (*.xtc *.trr);;All Files (*)"))
        xtc_layout.addWidget(self.edit_calc_xtc)
        xtc_layout.addWidget(btn_xtc)

        self.combo_analysis_type = QComboBox()
        self.combo_analysis_type.addItems([
            "RMSD (Root Mean Square Deviation)",
            "RMSF (Root Mean Square Fluctuation)",
            "Radius of Gyration (Gyrate)",
            "Hydrogen Bonds Count",
        ])

        self.combo_group = QComboBox()
        self.combo_group.addItems(["Backbone", "C-alpha", "Protein", "Ligand", "Protein-Ligand Complex"])

        btn_run_calc = QPushButton("📊 Calculate & Plot Property")
        btn_run_calc.setFont(QFont("Segoe UI", 9, QFont.Weight.Bold))
        btn_run_calc.clicked.connect(self.run_analysis)

        form.addRow("Structure/Topology (.tpr):", tpr_layout)
        form.addRow("Trajectory File (.xtc):", xtc_layout)
        form.addRow("Analysis Type:", self.combo_analysis_type)
        form.addRow("Atom Group Selection:", self.combo_group)
        form.addRow(btn_run_calc)

        return widget

    def _create_xvg_tab(self):
        widget = QWidget()
        form = QFormLayout(widget)

        xvg_layout = QHBoxLayout()
        self.edit_xvg_path = QLineEdit()
        self.edit_xvg_path.setPlaceholderText("Select any .xvg file on disk...")
        btn_xvg = QPushButton("Browse...")
        btn_xvg.clicked.connect(lambda: self.browse_file(self.edit_xvg_path, "XVG Data Files (*.xvg);;All Files (*)"))
        xvg_layout.addWidget(self.edit_xvg_path)
        xvg_layout.addWidget(btn_xvg)

        btn_plot_xvg = QPushButton("📈 Render XVG Plot")
        btn_plot_xvg.setFont(QFont("Segoe UI", 9, QFont.Weight.Bold))
        btn_plot_xvg.clicked.connect(self.plot_selected_xvg)

        form.addRow("XVG File Path:", xvg_layout)
        form.addRow(btn_plot_xvg)

        return widget

    def _create_external_tab(self):
        widget = QWidget()
        form = QFormLayout(widget)

        gro_layout = QHBoxLayout()
        self.edit_ext_gro = QLineEdit("complex.gro")
        btn_gro = QPushButton("Browse...")
        btn_gro.clicked.connect(lambda: self.browse_file(self.edit_ext_gro, "Structure Files (*.gro *.pdb);;All Files (*)"))
        gro_layout.addWidget(self.edit_ext_gro)
        gro_layout.addWidget(btn_gro)

        xtc_layout = QHBoxLayout()
        self.edit_ext_xtc = QLineEdit("md_noPBC.xtc")
        btn_xtc = QPushButton("Browse...")
        btn_xtc.clicked.connect(lambda: self.browse_file(self.edit_ext_xtc, "Trajectory Files (*.xtc *.trr);;All Files (*)"))
        xtc_layout.addWidget(self.edit_ext_xtc)
        xtc_layout.addWidget(btn_xtc)

        btn_vmd = QPushButton("👁️ Launch VMD (Visual Molecular Dynamics)")
        btn_vmd.clicked.connect(self.launch_vmd)

        btn_pymol = QPushButton("🧬 Launch PyMOL")
        btn_pymol.clicked.connect(self.launch_pymol)

        btn_chimera = QPushButton("🔬 Launch UCSF Chimera / ChimeraX")
        btn_chimera.clicked.connect(self.launch_chimera)

        form.addRow("Structure File:", gro_layout)
        form.addRow("Trajectory File:", xtc_layout)
        form.addRow(btn_vmd)
        form.addRow(btn_pymol)
        form.addRow(btn_chimera)

        return widget

    def run_trjconv(self):
        self.work_dir = self.lbl_workdir.text().strip()
        tpr = os.path.join(self.work_dir, self.edit_trj_tpr.text().strip())
        xtc = os.path.join(self.work_dir, self.edit_trj_xtc.text().strip())
        out_xtc = os.path.join(self.work_dir, self.edit_trj_out.text().strip())

        if not os.path.exists(tpr) or not os.path.exists(xtc):
            QMessageBox.warning(self, "Missing File", "Specified TPR or Trajectory file does not exist.")
            return

        pbc_mode = self.combo_pbc.currentText()
        center = self.cb_center.isChecked()

        cmd = f'echo "Protein System" | "{GMX_PATH}" trjconv -s "{tpr}" -f "{xtc}" -o "{out_xtc}" -pbc {pbc_mode}'
        if center:
            cmd += " -center"

        def on_trj_done(success, task, out):
            if success:
                QMessageBox.information(self, "PBC Removed", f"🎉 Processed trajectory saved to:\n{out_xtc}")

        self.console.clear_log()
        self.active_worker = GMXWorker([cmd], work_dir=self.work_dir, task_name="trjconv PBC Removal")
        self.active_worker.log_signal.connect(self.console.append_log)
        self.active_worker.finished_signal.connect(on_trj_done)
        self.active_worker.start()

    def run_analysis(self):
        self.work_dir = self.lbl_workdir.text().strip()
        tpr = os.path.join(self.work_dir, self.edit_calc_tpr.text().strip())
        xtc = os.path.join(self.work_dir, self.edit_calc_xtc.text().strip())

        if not os.path.exists(tpr) or not os.path.exists(xtc):
            QMessageBox.warning(self, "Missing File", "Specified TPR or Trajectory file does not exist.")
            return

        atype_text = self.combo_analysis_type.currentText()
        if "RMSD" in atype_text:
            atype = "rmsd"
            out_xvg = os.path.join(self.work_dir, "rmsd.xvg")
        elif "RMSF" in atype_text:
            atype = "rmsf"
            out_xvg = os.path.join(self.work_dir, "rmsf.xvg")
        elif "Gyrate" in atype_text:
            atype = "gyrate"
            out_xvg = os.path.join(self.work_dir, "gyrate.xvg")
        else:
            atype = "hbond"
            out_xvg = os.path.join(self.work_dir, "hbonds.xvg")

        grp = self.combo_group.currentText()
        grp_input = "4 4" if "Backbone" in grp else "1 1"

        cmd = f'echo "{grp_input}" | "{GMX_PATH}" {atype if atype != "gyrate" else "gyrate"} -s "{tpr}" -f "{xtc}" -o "{out_xvg}"'

        def on_calc_done(success, task, out):
            if success and os.path.exists(out_xvg):
                data = parse_xvg(out_xvg)
                self.canvas.plot_data(data)
                self.console.append_log(f"🎉 Plot rendered for {out_xvg}")

        self.console.clear_log()
        self.active_worker = GMXWorker([cmd], work_dir=self.work_dir, task_name=f"Calculate {atype.upper()}")
        self.active_worker.log_signal.connect(self.console.append_log)
        self.active_worker.finished_signal.connect(on_calc_done)
        self.active_worker.start()

    def plot_selected_xvg(self):
        file_path = self.edit_xvg_path.text().strip()
        if not file_path or not os.path.exists(file_path):
            QMessageBox.warning(self, "File Not Found", "Please select a valid .xvg file on disk.")
            return

        try:
            data = parse_xvg(file_path)
            self.canvas.plot_data(data)
            self.console.append_log(f"📈 Loaded and plotted XVG file: {file_path}")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to parse XVG file:\n{e}")

    def launch_vmd(self):
        cfg = load_config()
        vmd_bin = cfg.get("vmd_path", "vmd")
        gro = os.path.join(self.work_dir, self.edit_ext_gro.text().strip())
        xtc = os.path.join(self.work_dir, self.edit_ext_xtc.text().strip())

        cmd = [vmd_bin]
        if os.path.exists(gro):
            cmd.append(gro)
        if os.path.exists(xtc):
            cmd.append(xtc)

        try:
            subprocess.Popen(cmd)
            self.console.append_log(f"🚀 Launched VMD: {' '.join(cmd)}")
        except Exception as e:
            QMessageBox.critical(self, "Launch Failed", f"Could not launch VMD executable '{vmd_bin}':\n{e}")

    def launch_pymol(self):
        cfg = load_config()
        pymol_bin = cfg.get("pymol_path", "pymol")
        gro = os.path.join(self.work_dir, self.edit_ext_gro.text().strip())

        cmd = [pymol_bin]
        if os.path.exists(gro):
            cmd.append(gro)

        try:
            subprocess.Popen(cmd)
            self.console.append_log(f"🚀 Launched PyMOL: {' '.join(cmd)}")
        except Exception as e:
            QMessageBox.critical(self, "Launch Failed", f"Could not launch PyMOL executable '{pymol_bin}':\n{e}")

    def launch_chimera(self):
        cfg = load_config()
        chimera_bin = cfg.get("chimera_path", "chimera")
        gro = os.path.join(self.work_dir, self.edit_ext_gro.text().strip())

        cmd = [chimera_bin]
        if os.path.exists(gro):
            cmd.append(gro)

        try:
            subprocess.Popen(cmd)
            self.console.append_log(f"🚀 Launched Chimera: {' '.join(cmd)}")
        except Exception as e:
            QMessageBox.critical(self, "Launch Failed", f"Could not launch Chimera executable '{chimera_bin}':\n{e}")
