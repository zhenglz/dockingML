"""
Visual MDP Parameter Editor Dialog for GROMACS GUI
Enables graphical configuration of integrator, nsteps, emtol, coupling algorithms, and restraints.
"""

from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QHBoxLayout, QTabWidget, QWidget,
    QFormLayout, QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox,
    QPushButton, QLabel, QPlainTextEdit, QFileDialog, QMessageBox,
    QGroupBox
)
from PyQt6.QtGui import QFont, QFontDatabase
from PyQt6.QtCore import Qt

from src.core.mdp_editor import MDPEditor


class MDPDialog(QDialog):
    """Graphical MDP Configuration Editor Dialog."""

    def __init__(self, stage="em", initial_params=None, parent=None):
        super().__init__(parent)
        self.stage = stage
        self.params = initial_params if initial_params else MDPEditor.get_preset(stage)
        self.init_ui()

    def init_ui(self):
        self.setWindowTitle(f"⚙️ MDP Parameter Editor — [{self.stage.upper()} Stage]")
        self.resize(750, 550)

        layout = QVBoxLayout(self)

        # Tabs for parameter categories
        self.tabs = QTabWidget()

        # Tab 1: Run Control
        self.tabs.addTab(self._create_run_control_tab(), "🏃 Run Control")

        # Tab 2: Output Control
        self.tabs.addTab(self._create_output_control_tab(), "📡 Output Control")

        # Tab 3: Nonbonded & Electrostatics
        self.tabs.addTab(self._create_nonbonded_tab(), "⚡ Nonbonded")

        # Tab 4: Temperature & Pressure Coupling
        self.tabs.addTab(self._create_coupling_tab(), "🌡️ Thermostat & Barostat")

        # Tab 5: Constraints & Defines
        self.tabs.addTab(self._create_constraints_tab(), "🔒 Constraints & Restraints")

        # Tab 6: Raw Text Preview
        self.tabs.addTab(self._create_raw_preview_tab(), "📝 Raw MDP Preview")

        self.tabs.currentChanged.connect(self.update_raw_preview)

        layout.addWidget(self.tabs)

        # Bottom Buttons
        btn_layout = QHBoxLayout()

        load_btn = QPushButton("📂 Load File...")
        load_btn.clicked.connect(self.load_from_file)

        preset_btn = QPushButton("🔄 Reset Stage Preset")
        preset_btn.clicked.connect(self.reset_preset)

        save_as_btn = QPushButton("💾 Save To File...")
        save_as_btn.clicked.connect(self.save_to_file)

        apply_btn = QPushButton("✅ Apply & Close")
        apply_btn.setDefault(True)
        apply_btn.clicked.connect(self.apply_and_accept)

        cancel_btn = QPushButton("Cancel")
        cancel_btn.clicked.connect(self.reject)

        btn_layout.addWidget(load_btn)
        btn_layout.addWidget(preset_btn)
        btn_layout.addWidget(save_as_btn)
        btn_layout.addStretch()
        btn_layout.addWidget(cancel_btn)
        btn_layout.addWidget(apply_btn)

        layout.addLayout(btn_layout)

        # Populate form fields with current parameter dictionary
        self.load_fields_from_dict(self.params)

    def _create_run_control_tab(self):
        widget = QWidget()
        form = QFormLayout(widget)

        self.edit_title = QLineEdit()
        self.combo_integrator = QComboBox()
        self.combo_integrator.addItems(["steep", "cg", "md", "md-vv", "sd", "bd"])

        self.edit_nsteps = QLineEdit()
        self.edit_dt = QLineEdit()
        self.edit_emtol = QLineEdit()
        self.edit_emstep = QLineEdit()

        form.addRow("Title Comment:", self.edit_title)
        form.addRow("Integrator (integrator):", self.combo_integrator)
        form.addRow("Number of Steps (nsteps):", self.edit_nsteps)
        form.addRow("Time Step dt (ps):", self.edit_dt)
        form.addRow("EM Tolerance emtol (kJ/mol/nm):", self.edit_emtol)
        form.addRow("EM Initial Step Size emstep (nm):", self.edit_emstep)

        return widget

    def _create_output_control_tab(self):
        widget = QWidget()
        form = QFormLayout(widget)

        self.edit_nstxout_comp = QLineEdit()
        self.edit_nstenergy = QLineEdit()
        self.edit_nstlog = QLineEdit()

        form.addRow("Compressed Trajectory Saving Interval (nstxout-compressed):", self.edit_nstxout_comp)
        form.addRow("Energy Writing Interval (nstenergy):", self.edit_nstenergy)
        form.addRow("Log Output Writing Interval (nstlog):", self.edit_nstlog)

        return widget

    def _create_nonbonded_tab(self):
        widget = QWidget()
        form = QFormLayout(widget)

        self.combo_cutoff_scheme = QComboBox()
        self.combo_cutoff_scheme.addItems(["Verlet", "Group"])

        self.edit_nstlist = QLineEdit()

        self.combo_coulombtype = QComboBox()
        self.combo_coulombtype.addItems(["PME", "Cut-off", "Ewald", "Reaction-Field"])

        self.edit_rcoulomb = QLineEdit()
        self.edit_rvdw = QLineEdit()
        self.edit_pme_order = QLineEdit()
        self.edit_fourierspacing = QLineEdit()

        form.addRow("Cutoff Scheme (cutoff-scheme):", self.combo_cutoff_scheme)
        form.addRow("Neighbor List Update Frequency (nstlist):", self.edit_nstlist)
        form.addRow("Electrostatics Type (coulombtype):", self.combo_coulombtype)
        form.addRow("Coulomb Cutoff rcoulomb (nm):", self.edit_rcoulomb)
        form.addRow("VdW Cutoff rvdw (nm):", self.edit_rvdw)
        form.addRow("PME Order (pme_order):", self.edit_pme_order)
        form.addRow("Fourier Grid Spacing (fourierspacing nm):", self.edit_fourierspacing)

        return widget

    def _create_coupling_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        # Thermostat Group
        t_group = QGroupBox("🌡️ Temperature Coupling (Thermostat)")
        t_form = QFormLayout(t_group)

        self.combo_tcoupl = QComboBox()
        self.combo_tcoupl.addItems(["no", "V-rescale", "Nose-Hoover", "Berendsen", "andersen"])

        self.edit_tc_grps = QLineEdit()
        self.edit_tau_t = QLineEdit()
        self.edit_ref_t = QLineEdit()

        t_form.addRow("Temperature Coupling Type (Tcoupl):", self.combo_tcoupl)
        t_form.addRow("Coupling Groups (tc-grps):", self.edit_tc_grps)
        t_form.addRow("Time Constant tau_t (ps):", self.edit_tau_t)
        t_form.addRow("Reference Temp ref_t (K):", self.edit_ref_t)

        # Barostat Group
        p_group = QGroupBox("⏲️ Pressure Coupling (Barostat)")
        p_form = QFormLayout(p_group)

        self.combo_pcoupl = QComboBox()
        self.combo_pcoupl.addItems(["no", "C-rescale", "Parrinello-Rahman", "Berendsen", "MTTK"])

        self.combo_pcoupltype = QComboBox()
        self.combo_pcoupltype.addItems(["isotropic", "semiisotropic", "anisotropic", "surface-tension"])

        self.edit_tau_p = QLineEdit()
        self.edit_ref_p = QLineEdit()
        self.edit_compressibility = QLineEdit()

        p_form.addRow("Pressure Coupling Type (Pcoupl):", self.combo_pcoupl)
        p_form.addRow("Coupling Geometry (pcoupltype):", self.combo_pcoupltype)
        p_form.addRow("Time Constant tau_p (ps):", self.edit_tau_p)
        p_form.addRow("Reference Pressure ref_p (bar):", self.edit_ref_p)
        p_form.addRow("Isothermal Compressibility (bar^-1):", self.edit_compressibility)

        layout.addWidget(t_group)
        layout.addWidget(p_group)

        return widget

    def _create_constraints_tab(self):
        widget = QWidget()
        form = QFormLayout(widget)

        self.combo_constraints = QComboBox()
        self.combo_constraints.addItems(["none", "h-bonds", "all-bonds", "h-angles", "all-angles"])

        self.combo_algorithm = QComboBox()
        self.combo_algorithm.addItems(["lincs", "shake"])

        self.edit_define = QLineEdit()

        form.addRow("Bond Constraints (constraints):", self.combo_constraints)
        form.addRow("Constraint Algorithm:", self.combo_algorithm)
        form.addRow("Pre-processor Defines (define):", self.edit_define)

        return widget

    def _create_raw_preview_tab(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)

        self.preview_text = QPlainTextEdit()
        self.preview_text.setFont(QFont("Consolas", 9))
        self.preview_text.setStyleSheet(
            "QPlainTextEdit { background-color: #1e1e1e; color: #a9b7c6; font-family: 'Consolas', monospace; }"
        )

        layout.addWidget(QLabel("Raw MDP File Contents Preview:"))
        layout.addWidget(self.preview_text)

        return widget

    def load_fields_from_dict(self, d):
        """Populate GUI form fields from parameter dictionary."""
        self.edit_title.setText(d.get("title", f"GROMACS {self.stage.upper()} Stage"))

        # Combos helper
        def set_combo(combo, val):
            idx = combo.findText(str(val), Qt.MatchFlag.MatchExactly)
            if idx >= 0:
                combo.setCurrentIndex(idx)

        set_combo(self.combo_integrator, d.get("integrator", "steep"))
        self.edit_nsteps.setText(str(d.get("nsteps", "50000")))
        self.edit_dt.setText(str(d.get("dt", "0.002")))
        self.edit_emtol.setText(str(d.get("emtol", "500.0")))
        self.edit_emstep.setText(str(d.get("emstep", "0.01")))

        self.edit_nstxout_comp.setText(str(d.get("nstxout-compressed", "5000")))
        self.edit_nstenergy.setText(str(d.get("nstenergy", "5000")))
        self.edit_nstlog.setText(str(d.get("nstlog", "5000")))

        set_combo(self.combo_cutoff_scheme, d.get("cutoff-scheme", "Verlet"))
        self.edit_nstlist.setText(str(d.get("nstlist", "20")))
        set_combo(self.combo_coulombtype, d.get("coulombtype", "PME"))
        self.edit_rcoulomb.setText(str(d.get("rcoulomb", "1.0")))
        self.edit_rvdw.setText(str(d.get("rvdw", "1.0")))
        self.edit_pme_order.setText(str(d.get("pme_order", "4")))
        self.edit_fourierspacing.setText(str(d.get("fourierspacing", "0.12")))

        set_combo(self.combo_tcoupl, d.get("Tcoupl", "V-rescale"))
        self.edit_tc_grps.setText(str(d.get("tc-grps", "Protein_Ligand Water_and_ions")))
        self.edit_tau_t.setText(str(d.get("tau_t", "0.1 0.1")))
        self.edit_ref_t.setText(str(d.get("ref_t", "300 300")))

        set_combo(self.combo_pcoupl, d.get("Pcoupl", "no"))
        set_combo(self.combo_pcoupltype, d.get("pcoupltype", "isotropic"))
        self.edit_tau_p.setText(str(d.get("tau_p", "2.0")))
        self.edit_ref_p.setText(str(d.get("ref_p", "1.0")))
        self.edit_compressibility.setText(str(d.get("compressibility", "4.5e-5")))

        set_combo(self.combo_constraints, d.get("constraints", "h-bonds"))
        set_combo(self.combo_algorithm, d.get("constraint_algorithm", "lincs"))
        self.edit_define.setText(str(d.get("define", "")))

    def build_dict_from_fields(self):
        """Gather values from GUI form fields into parameter dictionary."""
        d = self.params.copy()

        if self.edit_title.text().strip():
            d["title"] = self.edit_title.text().strip()

        d["integrator"] = self.combo_integrator.currentText()
        d["nsteps"] = self.edit_nsteps.text().strip()
        if "dt" in d or self.stage != "em":
            d["dt"] = self.edit_dt.text().strip()
        if self.stage == "em" or self.combo_integrator.currentText() in ["steep", "cg"]:
            d["emtol"] = self.edit_emtol.text().strip()
            d["emstep"] = self.edit_emstep.text().strip()

        d["nstxout-compressed"] = self.edit_nstxout_comp.text().strip()
        d["nstenergy"] = self.edit_nstenergy.text().strip()
        d["nstlog"] = self.edit_nstlog.text().strip()

        d["cutoff-scheme"] = self.combo_cutoff_scheme.currentText()
        d["nstlist"] = self.edit_nstlist.text().strip()
        d["coulombtype"] = self.combo_coulombtype.currentText()
        d["rcoulomb"] = self.edit_rcoulomb.text().strip()
        d["rvdw"] = self.edit_rvdw.text().strip()
        d["pme_order"] = self.edit_pme_order.text().strip()
        d["fourierspacing"] = self.edit_fourierspacing.text().strip()

        d["Tcoupl"] = self.combo_tcoupl.currentText()
        if d["Tcoupl"] != "no":
            d["tc-grps"] = self.edit_tc_grps.text().strip()
            d["tau_t"] = self.edit_tau_t.text().strip()
            d["ref_t"] = self.edit_ref_t.text().strip()

        d["Pcoupl"] = self.combo_pcoupl.currentText()
        if d["Pcoupl"] != "no":
            d["pcoupltype"] = self.combo_pcoupltype.currentText()
            d["tau_p"] = self.edit_tau_p.text().strip()
            d["ref_p"] = self.edit_ref_p.text().strip()
            d["compressibility"] = self.edit_compressibility.text().strip()

        d["constraints"] = self.combo_constraints.currentText()
        d["constraint_algorithm"] = self.combo_algorithm.currentText()

        if self.edit_define.text().strip():
            d["define"] = self.edit_define.text().strip()
        elif "define" in d:
            del d["define"]

        self.params = d
        return d

    def update_raw_preview(self):
        d = self.build_dict_from_fields()
        raw_text = MDPEditor.generate_mdp_text(d, title=f"GROMACS MDP [{self.stage.upper()}]")
        self.preview_text.setPlainText(raw_text)

    def reset_preset(self):
        self.params = MDPEditor.get_preset(self.stage)
        self.load_fields_from_dict(self.params)
        QMessageBox.information(self, "Reset", f"Form reset to default {self.stage.upper()} stage preset.")

    def load_from_file(self):
        path, _ = QFileDialog.getOpenFileName(self, "Select MDP File", "", "MDP Files (*.mdp);;All Files (*)")
        if path:
            parsed = MDPEditor.parse_mdp(path)
            if parsed:
                self.params = parsed
                self.load_fields_from_dict(parsed)
                QMessageBox.information(self, "Loaded", f"Loaded parameters from:\n{path}")

    def save_to_file(self):
        d = self.build_dict_from_fields()
        path, _ = QFileDialog.getSaveFileName(self, "Save MDP File", f"{self.stage}.mdp", "MDP Files (*.mdp);;All Files (*)")
        if path:
            MDPEditor.save_mdp(d, path)
            QMessageBox.information(self, "Saved", f"MDP saved to:\n{path}")

    def apply_and_accept(self):
        self.build_dict_from_fields()
        self.accept()
