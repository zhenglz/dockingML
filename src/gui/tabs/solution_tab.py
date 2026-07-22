"""
Solution Simulator Module Tab
Provides solution system topology creation, solvation, box setup, and MD simulation execution.
"""

import os
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QFormLayout,
    QLineEdit, QPushButton, QComboBox, QDoubleSpinBox,
    QLabel, QFileDialog, QMessageBox, QProgressBar
)
from PyQt6.QtGui import QFont

from src.core.config import GMX_PATH, FORCEFIELDS, WATER_MODELS, BOX_SHAPES
from src.core.mdp_editor import MDPEditor
from src.core.gmx_runner import GMXWorker
from src.gui.widgets.log_console import LogConsole
from src.gui.mdp_dialog import MDPDialog


class SolutionSimulatorTab(QWidget):
    """Solution System Simulator Tab interface."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.work_dir = os.getcwd()
        self.mdp_params = {
            "em": MDPEditor.get_preset("em"),
            "nvt": MDPEditor.get_preset("nvt"),
            "npt": MDPEditor.get_preset("npt"),
            "md": MDPEditor.get_preset("md"),
        }
        self.active_worker = None
        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout(self)

        # Working Directory
        workdir_layout = QHBoxLayout()
        self.lbl_workdir = QLineEdit(self.work_dir)
        btn_browse_workdir = QPushButton("📂 Set Working Directory...")
        btn_browse_workdir.clicked.connect(self.browse_working_dir)
        workdir_layout.addWidget(QLabel("📁 Working Directory:"))
        workdir_layout.addWidget(self.lbl_workdir)
        workdir_layout.addWidget(btn_browse_workdir)
        main_layout.addLayout(workdir_layout)

        # Step 1: Solute Input & Topology Setup
        step1_box = QGroupBox("Step 1: Solute Topology & Water Box Setup")
        step1_form = QFormLayout(step1_box)

        solute_layout = QHBoxLayout()
        self.edit_solute = QLineEdit()
        self.edit_solute.setPlaceholderText("Select solute PDB/GRO file (or leave empty for pure water box)...")
        btn_solute = QPushButton("Browse...")
        btn_solute.clicked.connect(lambda: self.browse_file(self.edit_solute, "Structure Files (*.pdb *.gro);;All Files (*)"))
        solute_layout.addWidget(self.edit_solute)
        solute_layout.addWidget(btn_solute)

        self.combo_ff = QComboBox()
        for ff in FORCEFIELDS:
            self.combo_ff.addItem(ff["name"], ff["id"])

        self.combo_water = QComboBox()
        for w in WATER_MODELS:
            self.combo_water.addItem(w["name"], w["id"])

        btn_step1_run = QPushButton("⚡ Generate Solute Topology & Setup System (pdb2gmx / topol.top)")
        btn_step1_run.setFont(QFont("Segoe UI", 9, QFont.Weight.Bold))
        btn_step1_run.clicked.connect(self.run_step1_setup)

        step1_form.addRow("Solute File (PDB/GRO):", solute_layout)
        step1_form.addRow("Force Field:", self.combo_ff)
        step1_form.addRow("Water Model:", self.combo_water)
        step1_form.addRow(btn_step1_run)

        main_layout.addWidget(step1_box)

        # Step 2: Box, Solvation & Neutralization
        step2_box = QGroupBox("Step 2: Box Geometry, Solvation & Ion Neutralization")
        step2_form = QFormLayout(step2_box)

        self.combo_box_type = QComboBox()
        for b in BOX_SHAPES:
            self.combo_box_type.addItem(b["name"], b["id"])

        self.spin_box_dist = QDoubleSpinBox()
        self.spin_box_dist.setRange(0.5, 5.0)
        self.spin_box_dist.setValue(1.0)
        self.spin_box_dist.setSuffix(" nm")

        self.spin_salt_conc = QDoubleSpinBox()
        self.spin_salt_conc.setRange(0.0, 2.0)
        self.spin_salt_conc.setValue(0.15)
        self.spin_salt_conc.setSuffix(" M")

        btn_step2_run = QPushButton("💧 Solvate & Add Ions (editconf + solvate + genion)")
        btn_step2_run.setFont(QFont("Segoe UI", 9, QFont.Weight.Bold))
        btn_step2_run.clicked.connect(self.run_step2_solvation)

        step2_form.addRow("Box Shape:", self.combo_box_type)
        step2_form.addRow("Box Buffer Clearance:", self.spin_box_dist)
        step2_form.addRow("Salt Concentration:", self.spin_salt_conc)
        step2_form.addRow(btn_step2_run)

        main_layout.addWidget(step2_box)

        # Step 3: MDP Parameter Configuration & Pipeline Actions
        step3_box = QGroupBox("Step 3: Solution Simulation Pipeline Controls")
        step3_layout = QVBoxLayout(step3_box)

        cfg_layout = QHBoxLayout()
        btn_cfg_em = QPushButton("⚙️ Edit EM.mdp")
        btn_cfg_em.clicked.connect(lambda: self.open_mdp_editor("em"))

        btn_cfg_nvt = QPushButton("⚙️ Edit NVT.mdp")
        btn_cfg_nvt.clicked.connect(lambda: self.open_mdp_editor("nvt"))

        btn_cfg_npt = QPushButton("⚙️ Edit NPT.mdp")
        btn_cfg_npt.clicked.connect(lambda: self.open_mdp_editor("npt"))

        btn_cfg_md = QPushButton("⚙️ Edit MD.mdp")
        btn_cfg_md.clicked.connect(lambda: self.open_mdp_editor("md"))

        cfg_layout.addWidget(btn_cfg_em)
        cfg_layout.addWidget(btn_cfg_nvt)
        cfg_layout.addWidget(btn_cfg_npt)
        cfg_layout.addWidget(btn_cfg_md)

        run_layout = QHBoxLayout()
        btn_run_em = QPushButton("▶️ Run EM")
        btn_run_em.clicked.connect(self.run_em)

        btn_run_nvt = QPushButton("▶️ Run NVT")
        btn_run_nvt.clicked.connect(self.run_nvt)

        btn_run_npt = QPushButton("▶️ Run NPT")
        btn_run_npt.clicked.connect(self.run_npt)

        btn_run_md = QPushButton("▶️ Run Production MD")
        btn_run_md.clicked.connect(self.run_md)

        btn_run_all = QPushButton("🚀 Run Complete Solution Pipeline")
        btn_run_all.setFont(QFont("Segoe UI", 9, QFont.Weight.Bold))
        btn_run_all.setStyleSheet("background-color: #2e7d32; color: white; border-radius: 4px; padding: 6px;")
        btn_run_all.clicked.connect(self.run_all_stages)

        run_layout.addWidget(btn_run_em)
        run_layout.addWidget(btn_run_nvt)
        run_layout.addWidget(btn_run_npt)
        run_layout.addWidget(btn_run_md)
        run_layout.addWidget(btn_run_all)

        step3_layout.addLayout(cfg_layout)
        step3_layout.addLayout(run_layout)

        main_layout.addWidget(step3_box)

        # Progress Bar & Console
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        main_layout.addWidget(self.progress_bar)

        self.console = LogConsole()
        main_layout.addWidget(self.console)

    def browse_working_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Working Directory", self.work_dir)
        if dir_path:
            self.work_dir = dir_path
            self.lbl_workdir.setText(dir_path)

    def browse_file(self, target_line_edit, file_filter):
        file_path, _ = QFileDialog.getOpenFileName(self, "Select File", self.work_dir, file_filter)
        if file_path:
            target_line_edit.setText(file_path)

    def open_mdp_editor(self, stage):
        dlg = MDPDialog(stage=stage, initial_params=self.mdp_params[stage], parent=self)
        if dlg.exec() == MDPDialog.DialogCode.Accepted:
            self.mdp_params[stage] = dlg.params
            self.console.append_log(f"✅ Updated {stage.upper()} parameter config.")

    def run_step1_setup(self):
        self.work_dir = self.lbl_workdir.text().strip()
        os.makedirs(self.work_dir, exist_ok=True)

        solute_pdb = self.edit_solute.text().strip()
        ff_id = self.combo_ff.currentData()
        water_id = self.combo_water.currentData()

        if solute_pdb and os.path.exists(solute_pdb):
            cmd = [
                GMX_PATH, "pdb2gmx",
                "-f", solute_pdb,
                "-o", os.path.join(self.work_dir, "solute.gro"),
                "-p", os.path.join(self.work_dir, "topol.top"),
                "-ff", ff_id,
                "-water", water_id,
                "-ignh"
            ]
            self._execute_worker([cmd], task_name="Step 1: Solute pdb2gmx")
        else:
            self.console.append_log("ℹ️ No solute input provided; assuming pure solvent/water box pipeline setup.")

    def run_step2_solvation(self):
        self.work_dir = self.lbl_workdir.text().strip()
        solute_gro = os.path.join(self.work_dir, "solute.gro")
        topol_top = os.path.join(self.work_dir, "topol.top")

        in_struct = solute_gro if os.path.exists(solute_gro) else None
        box_type = self.combo_box_type.currentData()
        box_dist = self.spin_box_dist.value()
        salt_conc = self.spin_salt_conc.value()

        cmds = []
        box_gro = os.path.join(self.work_dir, "box.gro")
        solv_gro = os.path.join(self.work_dir, "solv.gro")
        ions_tpr = os.path.join(self.work_dir, "ions.tpr")
        ions_gro = os.path.join(self.work_dir, "ions.gro")
        em_sol_mdp = os.path.join(self.work_dir, "em_sol.mdp")

        MDPEditor.save_mdp(MDPEditor.get_preset("em"), em_sol_mdp)

        if in_struct:
            cmds.append([GMX_PATH, "editconf", "-f", in_struct, "-o", box_gro, "-bt", box_type, "-d", str(box_dist), "-c"])
            cmds.append([GMX_PATH, "solvate", "-cp", box_gro, "-cs", "spc216.gro", "-o", solv_gro, "-p", topol_top])
        else:
            # Create pure water box
            cmds.append([GMX_PATH, "solvate", "-cs", "spc216.gro", "-o", solv_gro, "-box", str(box_dist), str(box_dist), str(box_dist), "-p", topol_top])

        cmds.append([GMX_PATH, "grompp", "-f", em_sol_mdp, "-c", solv_gro, "-p", topol_top, "-o", ions_tpr, "-maxwarn", "2"])
        cmds.append(f'echo "SOL" | "{GMX_PATH}" genion -s "{ions_tpr}" -o "{ions_gro}" -p "{topol_top}" -pname NA -nname CL -neutral -conc {salt_conc}')

        self._execute_worker(cmds, task_name="Step 2: Solvation & Ions")

    def run_em(self):
        self._run_stage("em", "ions.gro", "em")

    def run_nvt(self):
        self._run_stage("nvt", "em.gro", "nvt")

    def run_npt(self):
        self._run_stage("npt", "nvt.gro", "npt")

    def run_md(self):
        self._run_stage("md", "npt.gro", "md")

    def _run_stage(self, stage, in_gro, out_prefix):
        self.work_dir = self.lbl_workdir.text().strip()
        in_path = os.path.join(self.work_dir, in_gro)
        topol = os.path.join(self.work_dir, "topol.top")

        mdp_file = os.path.join(self.work_dir, f"{stage}.mdp")
        tpr_file = os.path.join(self.work_dir, f"{out_prefix}.tpr")

        MDPEditor.save_mdp(self.mdp_params[stage], mdp_file)

        cmds = [
            [GMX_PATH, "grompp", "-f", mdp_file, "-c", in_path, "-r", in_path, "-p", topol, "-o", tpr_file, "-maxwarn", "2"],
            [GMX_PATH, "mdrun", "-deffnm", out_prefix, "-v"],
        ]
        self._execute_worker(cmds, task_name=f"Solution Stage: {stage.upper()}")

    def run_all_stages(self):
        self.work_dir = self.lbl_workdir.text().strip()
        for s in ["em", "nvt", "npt", "md"]:
            MDPEditor.save_mdp(self.mdp_params[s], os.path.join(self.work_dir, f"{s}.mdp"))

        cmds = [
            [GMX_PATH, "grompp", "-f", "em.mdp", "-c", "ions.gro", "-p", "topol.top", "-o", "em.tpr", "-maxwarn", "2"],
            [GMX_PATH, "mdrun", "-deffnm", "em", "-v"],
            [GMX_PATH, "grompp", "-f", "nvt.mdp", "-c", "em.gro", "-r", "em.gro", "-p", "topol.top", "-o", "nvt.tpr", "-maxwarn", "2"],
            [GMX_PATH, "mdrun", "-deffnm", "nvt", "-v"],
            [GMX_PATH, "grompp", "-f", "npt.mdp", "-c", "nvt.gro", "-r", "nvt.gro", "-t", "nvt.cpt", "-p", "topol.top", "-o", "npt.tpr", "-maxwarn", "2"],
            [GMX_PATH, "mdrun", "-deffnm", "npt", "-v"],
            [GMX_PATH, "grompp", "-f", "md.mdp", "-c", "npt.gro", "-t", "npt.cpt", "-p", "topol.top", "-o", "md.tpr", "-maxwarn", "2"],
            [GMX_PATH, "mdrun", "-deffnm", "md", "-v"],
        ]
        self._execute_worker(cmds, task_name="Solution Complete Workflow")

    def _execute_worker(self, commands, task_name="Task"):
        if self.active_worker and self.active_worker.isRunning():
            QMessageBox.warning(self, "Task Active", "A task is currently executing.")
            return

        self.console.clear_log()
        self.progress_bar.setValue(0)

        self.active_worker = GMXWorker(commands, work_dir=self.work_dir, task_name=task_name)
        self.active_worker.log_signal.connect(self.console.append_log)
        self.active_worker.progress_signal.connect(self.progress_bar.setValue)
        self.active_worker.start()
