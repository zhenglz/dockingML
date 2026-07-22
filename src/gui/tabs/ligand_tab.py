"""
Ligand Simulator Module Tab
Provides complete protein-ligand complex pipeline setup and workflow execution.
"""

import os
from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QFormLayout,
    QLineEdit, QPushButton, QComboBox, QDoubleSpinBox, QSpinBox,
    QLabel, QFileDialog, QMessageBox, QProgressBar
)
from PyQt6.QtGui import QFont

from src.core.config import GMX_PATH, FORCEFIELDS, WATER_MODELS, BOX_SHAPES
from src.core.mdp_editor import MDPEditor
from src.core.topology_builder import merge_gro_files, update_topology, parse_ligand_resname_from_gro
from src.core.gmx_runner import GMXWorker
from src.gui.widgets.log_console import LogConsole
from src.gui.mdp_dialog import MDPDialog


class LigandSimulatorTab(QWidget):
    """Protein-Ligand Complex Simulator Tab interface."""

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

        # Top Group: Workspace Working Directory
        workdir_layout = QHBoxLayout()
        self.lbl_workdir = QLineEdit(self.work_dir)
        btn_browse_workdir = QPushButton("📂 Set Working Directory...")
        btn_browse_workdir.clicked.connect(self.browse_working_dir)
        workdir_layout.addWidget(QLabel("📁 Working Directory:"))
        workdir_layout.addWidget(self.lbl_workdir)
        workdir_layout.addWidget(btn_browse_workdir)
        main_layout.addLayout(workdir_layout)

        # Step 1: Receptor & Ligand Topology Preparation
        step1_box = QGroupBox("Step 1: Receptor & Ligand Preparation")
        step1_form = QFormLayout(step1_box)

        # Receptor PDB File
        rec_layout = QHBoxLayout()
        self.edit_receptor = QLineEdit()
        self.edit_receptor.setPlaceholderText("Select receptor protein .pdb file...")
        btn_rec = QPushButton("Browse...")
        btn_rec.clicked.connect(lambda: self.browse_file(self.edit_receptor, "PDB Files (*.pdb);;All Files (*)"))
        rec_layout.addWidget(self.edit_receptor)
        rec_layout.addWidget(btn_rec)

        # Ligand Structure & Topology
        lig_struct_layout = QHBoxLayout()
        self.edit_ligand_struct = QLineEdit()
        self.edit_ligand_struct.setPlaceholderText("Select ligand structure .gro or .pdb file...")
        btn_lig_s = QPushButton("Browse...")
        btn_lig_s.clicked.connect(lambda: self.browse_file(self.edit_ligand_struct, "Structure Files (*.gro *.pdb);;All Files (*)"))
        lig_struct_layout.addWidget(self.edit_ligand_struct)
        lig_struct_layout.addWidget(btn_lig_s)

        lig_itp_layout = QHBoxLayout()
        self.edit_ligand_itp = QLineEdit()
        self.edit_ligand_itp.setPlaceholderText("Select ligand topology .itp file (from CGenFF / ACPYPE / OpenFF)...")
        btn_lig_i = QPushButton("Browse...")
        btn_lig_i.clicked.connect(lambda: self.browse_file(self.edit_ligand_itp, "ITP Files (*.itp);;All Files (*)"))
        lig_itp_layout.addWidget(self.edit_ligand_itp)
        lig_itp_layout.addWidget(btn_lig_i)

        # Force field & Water model
        self.combo_ff = QComboBox()
        for ff in FORCEFIELDS:
            self.combo_ff.addItem(ff["name"], ff["id"])

        self.combo_water = QComboBox()
        for w in WATER_MODELS:
            self.combo_water.addItem(w["name"], w["id"])

        btn_step1_run = QPushButton("⚡ Process Receptor & Build Complex (pdb2gmx + merge complex.gro + update topol.top)")
        btn_step1_run.setFont(QFont("Segoe UI", 9, QFont.Weight.Bold))
        btn_step1_run.clicked.connect(self.run_step1_build_complex)

        step1_form.addRow("Receptor PDB:", rec_layout)
        step1_form.addRow("Ligand Structure (.gro/.pdb):", lig_struct_layout)
        step1_form.addRow("Ligand Topology (.itp):", lig_itp_layout)
        step1_form.addRow("Protein Force Field:", self.combo_ff)
        step1_form.addRow("Water Model:", self.combo_water)
        step1_form.addRow(btn_step1_run)

        main_layout.addWidget(step1_box)

        # Step 2: Solvation & Neutralization
        step2_box = QGroupBox("Step 2: Box Definition, Solvation & Neutralization")
        step2_form = QFormLayout(step2_box)

        self.combo_box_type = QComboBox()
        for b in BOX_SHAPES:
            self.combo_box_type.addItem(b["name"], b["id"])

        self.spin_box_dist = QDoubleSpinBox()
        self.spin_box_dist.setRange(0.5, 5.0)
        self.spin_box_dist.setValue(1.0)
        self.spin_box_dist.setSingleStep(0.1)
        self.spin_box_dist.setSuffix(" nm")

        self.spin_salt_conc = QDoubleSpinBox()
        self.spin_salt_conc.setRange(0.0, 2.0)
        self.spin_salt_conc.setValue(0.15)
        self.spin_salt_conc.setSingleStep(0.05)
        self.spin_salt_conc.setSuffix(" M")

        btn_step2_run = QPushButton("💧 Solvate & Add Neutralizing Ions (editconf + solvate + genion)")
        btn_step2_run.setFont(QFont("Segoe UI", 9, QFont.Weight.Bold))
        btn_step2_run.clicked.connect(self.run_step2_solvate_and_ionize)

        step2_form.addRow("Box Geometry:", self.combo_box_type)
        step2_form.addRow("Box Buffer Distance:", self.spin_box_dist)
        step2_form.addRow("Neutralization Salt Concentration:", self.spin_salt_conc)
        step2_form.addRow(btn_step2_run)

        main_layout.addWidget(step2_box)

        # Step 3: Simulation Workflow Stage Controls
        step3_box = QGroupBox("Step 3: Molecular Dynamics Execution Pipeline")
        step3_layout = QVBoxLayout(step3_box)

        pipeline_btn_layout = QHBoxLayout()

        # Config buttons
        btn_cfg_em = QPushButton("⚙️ Edit EM.mdp")
        btn_cfg_em.clicked.connect(lambda: self.open_mdp_editor("em"))

        btn_cfg_nvt = QPushButton("⚙️ Edit NVT.mdp")
        btn_cfg_nvt.clicked.connect(lambda: self.open_mdp_editor("nvt"))

        btn_cfg_npt = QPushButton("⚙️ Edit NPT.mdp")
        btn_cfg_npt.clicked.connect(lambda: self.open_mdp_editor("npt"))

        btn_cfg_md = QPushButton("⚙️ Edit MD.mdp")
        btn_cfg_md.clicked.connect(lambda: self.open_mdp_editor("md"))

        pipeline_btn_layout.addWidget(btn_cfg_em)
        pipeline_btn_layout.addWidget(btn_cfg_nvt)
        pipeline_btn_layout.addWidget(btn_cfg_npt)
        pipeline_btn_layout.addWidget(btn_cfg_md)

        run_actions_layout = QHBoxLayout()

        btn_run_em = QPushButton("▶️ Run EM")
        btn_run_em.clicked.connect(self.run_em_stage)

        btn_run_nvt = QPushButton("▶️ Run NVT")
        btn_run_nvt.clicked.connect(self.run_nvt_stage)

        btn_run_npt = QPushButton("▶️ Run NPT")
        btn_run_npt.clicked.connect(self.run_npt_stage)

        btn_run_md = QPushButton("▶️ Run Production MD")
        btn_run_md.clicked.connect(self.run_md_stage)

        btn_run_all = QPushButton("🚀 Run Complete Workflow (EM ➔ NVT ➔ NPT ➔ MD)")
        btn_run_all.setFont(QFont("Segoe UI", 9, QFont.Weight.Bold))
        btn_run_all.setStyleSheet("background-color: #2e7d32; color: white; border-radius: 4px; padding: 6px;")
        btn_run_all.clicked.connect(self.run_complete_pipeline)

        run_actions_layout.addWidget(btn_run_em)
        run_actions_layout.addWidget(btn_run_nvt)
        run_actions_layout.addWidget(btn_run_npt)
        run_actions_layout.addWidget(btn_run_md)
        run_actions_layout.addWidget(btn_run_all)

        step3_layout.addLayout(pipeline_btn_layout)
        step3_layout.addLayout(run_actions_layout)

        main_layout.addWidget(step3_box)

        # Progress Bar
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        self.progress_bar.setValue(0)
        main_layout.addWidget(self.progress_bar)

        # Console Log Window
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
            self.console.append_log(f"✅ Updated {stage.upper()} parameter configuration.")

    def run_step1_build_complex(self):
        receptor_pdb = self.edit_receptor.text().strip()
        ligand_struct = self.edit_ligand_struct.text().strip()
        ligand_itp = self.edit_ligand_itp.text().strip()

        if not os.path.exists(receptor_pdb):
            QMessageBox.warning(self, "Missing Input", "Please select a valid receptor PDB file.")
            return
        if not os.path.exists(ligand_struct):
            QMessageBox.warning(self, "Missing Input", "Please select a valid ligand structure (.gro or .pdb) file.")
            return
        if not os.path.exists(ligand_itp):
            QMessageBox.warning(self, "Missing Input", "Please select a valid ligand topology (.itp) file.")
            return

        self.work_dir = self.lbl_workdir.text().strip()
        os.makedirs(self.work_dir, exist_ok=True)

        ff_id = self.combo_ff.currentData()
        water_id = self.combo_water.currentData()

        self.console.append_log("🚀 Step 1: Processing Receptor and Merging Complex Topology...")

        # Step 1 Commands
        # 1. Run pdb2gmx on receptor
        rec_gro = os.path.join(self.work_dir, "receptor.gro")
        topol_top = os.path.join(self.work_dir, "topol.top")

        cmd_pdb2gmx = [
            GMX_PATH, "pdb2gmx",
            "-f", receptor_pdb,
            "-o", rec_gro,
            "-p", topol_top,
            "-ff", ff_id,
            "-water", water_id,
            "-ignh"
        ]

        def post_process_complex(success, task_name, output):
            if not success:
                return

            try:
                # Merge GRO files
                complex_gro = os.path.join(self.work_dir, "complex.gro")
                merge_gro_files(rec_gro, ligand_struct, complex_gro)
                self.console.append_log(f"🎉 Created merged complex structure: {complex_gro}")

                # Copy ligand itp into working directory if needed
                local_itp = os.path.join(self.work_dir, os.path.basename(ligand_itp))
                if os.path.abspath(ligand_itp) != os.path.abspath(local_itp):
                    import shutil
                    shutil.copy(ligand_itp, local_itp)

                # Update topology topol.top
                resname = parse_ligand_resname_from_gro(ligand_struct)
                update_topology(topol_top, local_itp, ligand_resname=resname, count=1)
                self.console.append_log(f"🎉 Updated topology file '{topol_top}' with ligand '{resname}' include!")
                QMessageBox.information(self, "Step 1 Complete", "Complex structure (complex.gro) and topology (topol.top) built successfully!")

            except Exception as e:
                self.console.append_log(f"❌ Error during complex building: {e}")
                QMessageBox.critical(self, "Error", f"Failed to build complex:\n{e}")

        self._execute_worker([cmd_pdb2gmx], task_name="Step 1: Build Receptor & Topology", on_finish=post_process_complex)

    def run_step2_solvate_and_ionize(self):
        self.work_dir = self.lbl_workdir.text().strip()
        complex_gro = os.path.join(self.work_dir, "complex.gro")
        topol_top = os.path.join(self.work_dir, "topol.top")

        if not os.path.exists(complex_gro) or not os.path.exists(topol_top):
            QMessageBox.warning(self, "Missing Files", "complex.gro or topol.top not found in working directory. Please run Step 1 first.")
            return

        box_type = self.combo_box_type.currentData()
        box_dist = self.spin_box_dist.value()
        salt_conc = self.spin_salt_conc.value()

        # Step 2 commands
        box_gro = os.path.join(self.work_dir, "box.gro")
        solv_gro = os.path.join(self.work_dir, "solv.gro")
        ions_tpr = os.path.join(self.work_dir, "ions.tpr")
        ions_gro = os.path.join(self.work_dir, "ions.gro")
        em_sol_mdp = os.path.join(self.work_dir, "em_sol.mdp")

        # Save temporary em_sol.mdp for grompp
        MDPEditor.save_mdp(MDPEditor.get_preset("em"), em_sol_mdp, "EM MDP for Genion")

        cmd_editconf = [GMX_PATH, "editconf", "-f", complex_gro, "-o", box_gro, "-bt", box_type, "-d", str(box_dist), "-c"]
        cmd_solvate = [GMX_PATH, "solvate", "-cp", box_gro, "-cs", "spc216.gro", "-o", solv_gro, "-p", topol_top]
        cmd_grompp = [GMX_PATH, "grompp", "-f", em_sol_mdp, "-c", solv_gro, "-p", topol_top, "-o", ions_tpr, "-maxwarn", "2"]
        cmd_genion = f'echo "SOL" | "{GMX_PATH}" genion -s "{ions_tpr}" -o "{ions_gro}" -p "{topol_top}" -pname NA -nname CL -neutral -conc {salt_conc}'

        cmds = [cmd_editconf, cmd_solvate, cmd_grompp, cmd_genion]

        def post_process_solvation(success, task_name, output):
            if success:
                QMessageBox.information(self, "Step 2 Complete", "System box created, solvated, and ions neutralised (ions.gro) successfully!")

        self._execute_worker(cmds, task_name="Step 2: Solvation & Ions", on_finish=post_process_solvation)

    def run_em_stage(self):
        self._run_single_stage("em", in_gro="ions.gro", out_prefix="em")

    def run_nvt_stage(self):
        self._run_single_stage("nvt", in_gro="em.gro", out_prefix="nvt")

    def run_npt_stage(self):
        self._run_single_stage("npt", in_gro="nvt.gro", out_prefix="npt")

    def run_md_stage(self):
        self._run_single_stage("md", in_gro="npt.gro", out_prefix="md")

    def _run_single_stage(self, stage, in_gro, out_prefix):
        self.work_dir = self.lbl_workdir.text().strip()
        topol = os.path.join(self.work_dir, "topol.top")
        in_gro_path = os.path.join(self.work_dir, in_gro)

        if not os.path.exists(in_gro_path):
            QMessageBox.warning(self, "Missing Structure", f"Input file '{in_gro}' not found. Please run previous step first.")
            return

        mdp_file = os.path.join(self.work_dir, f"{stage}.mdp")
        tpr_file = os.path.join(self.work_dir, f"{out_prefix}.tpr")
        out_gro = os.path.join(self.work_dir, f"{out_prefix}.gro")

        MDPEditor.save_mdp(self.mdp_params[stage], mdp_file, title=f"{stage.upper()} Simulation Stage")

        cmd_grompp = [GMX_PATH, "grompp", "-f", mdp_file, "-c", in_gro_path, "-r", in_gro_path, "-p", topol, "-o", tpr_file, "-maxwarn", "2"]
        cmd_mdrun = [GMX_PATH, "mdrun", "-deffnm", out_prefix, "-v"]

        def on_stage_done(success, task, out):
            if success:
                self.console.append_log(f"🎉 Stage '{stage.upper()}' finished! Output gro: {out_gro}")

        self._execute_worker([cmd_grompp, cmd_mdrun], task_name=f"Stage: {stage.upper()}", on_finish=on_stage_done)

    def run_complete_pipeline(self):
        """Execute complete pipeline sequence from EM to NVT to NPT to Production MD."""
        self.work_dir = self.lbl_workdir.text().strip()
        ions_gro = os.path.join(self.work_dir, "ions.gro")
        topol = os.path.join(self.work_dir, "topol.top")

        if not os.path.exists(ions_gro):
            QMessageBox.warning(self, "Missing Ions GRO", "ions.gro not found in working directory. Please complete Step 1 & Step 2 first.")
            return

        # Prepare all MDP files
        for stage in ["em", "nvt", "npt", "md"]:
            mdp_path = os.path.join(self.work_dir, f"{stage}.mdp")
            MDPEditor.save_mdp(self.mdp_params[stage], mdp_path, title=f"{stage.upper()} Stage Config")

        cmds = [
            # 1. EM
            [GMX_PATH, "grompp", "-f", "em.mdp", "-c", "ions.gro", "-p", "topol.top", "-o", "em.tpr", "-maxwarn", "2"],
            [GMX_PATH, "mdrun", "-deffnm", "em", "-v"],
            # 2. NVT
            [GMX_PATH, "grompp", "-f", "nvt.mdp", "-c", "em.gro", "-r", "em.gro", "-p", "topol.top", "-o", "nvt.tpr", "-maxwarn", "2"],
            [GMX_PATH, "mdrun", "-deffnm", "nvt", "-v"],
            # 3. NPT
            [GMX_PATH, "grompp", "-f", "npt.mdp", "-c", "nvt.gro", "-r", "nvt.gro", "-t", "nvt.cpt", "-p", "topol.top", "-o", "npt.tpr", "-maxwarn", "2"],
            [GMX_PATH, "mdrun", "-deffnm", "npt", "-v"],
            # 4. Production MD
            [GMX_PATH, "grompp", "-f", "md.mdp", "-c", "npt.gro", "-t", "npt.cpt", "-p", "topol.top", "-o", "md.tpr", "-maxwarn", "2"],
            [GMX_PATH, "mdrun", "-deffnm", "md", "-v"],
        ]

        def on_all_done(success, task, summary):
            if success:
                QMessageBox.information(self, "Pipeline Complete", "🎉 Complete Protein-Ligand MD Workflow finished successfully!")

        self._execute_worker(cmds, task_name="Complete Pipeline Execution", on_finish=on_all_done)

    def _execute_worker(self, commands, task_name="Task", on_finish=None):
        if self.active_worker and self.active_worker.isRunning():
            QMessageBox.warning(self, "Task Running", "A task is currently executing. Please wait or cancel it before starting another.")
            return

        self.console.clear_log()
        self.progress_bar.setValue(0)

        self.active_worker = GMXWorker(commands, work_dir=self.work_dir, task_name=task_name)
        self.active_worker.log_signal.connect(self.console.append_log)
        self.active_worker.progress_signal.connect(self.progress_bar.setValue)

        if on_finish:
            self.active_worker.finished_signal.connect(on_finish)

        self.active_worker.start()
