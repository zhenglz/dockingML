"""
Async Process Runner for GROMACS Commands
Uses PyQt6 QThread to execute shell commands non-blockingly with real-time log output streaming.
"""

import os
import subprocess
from PyQt6.QtCore import QThread, pyqtSignal


class GMXWorker(QThread):
    """Worker thread for running single commands or pipelined GROMACS workflow steps."""

    log_signal = pyqtSignal(str)          # Emits live stdout/stderr log lines
    error_signal = pyqtSignal(str)        # Emits error messages
    progress_signal = pyqtSignal(int)     # Emits integer percentage or stage index
    finished_signal = pyqtSignal(bool, str, str)  # Emits (success, task_name, output_summary)

    def __init__(self, command_list, work_dir=None, task_name="GROMACS Task", parent=None):
        super().__init__(parent)
        self.commands = command_list if isinstance(command_list, list) else [command_list]
        self.work_dir = work_dir or os.getcwd()
        self.task_name = task_name
        self._is_cancelled = False
        self.process = None

    def run(self):
        """Execute commands in order."""
        total_steps = len(self.commands)
        self.log_signal.emit(f"🚀 Starting Task: {self.task_name}\n📁 Directory: {self.work_dir}\n" + "─" * 60)

        for idx, cmd in enumerate(self.commands):
            if self._is_cancelled:
                self.log_signal.emit("⚠️ Task was cancelled by user.")
                self.finished_signal.emit(False, self.task_name, "Task cancelled by user.")
                return

            if isinstance(cmd, str):
                cmd_args = cmd
                shell = True
            else:
                cmd_args = [str(x) for x in cmd]
                shell = False

            cmd_str = cmd if isinstance(cmd, str) else " ".join(cmd_args)
            self.log_signal.emit(f"▶️ [{idx + 1}/{total_steps}] Executing: {cmd_str}")

            try:
                self.process = subprocess.Popen(
                    cmd_args,
                    cwd=self.work_dir,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                    universal_newlines=True,
                    shell=shell,
                )

                if self.process.stdout:
                    for line in iter(self.process.stdout.readline, ""):
                        if self._is_cancelled:
                            self.process.terminate()
                            self.log_signal.emit("⚠️ Terminating process...")
                            break
                        self.log_signal.emit(line.rstrip())

                self.process.wait()
                return_code = self.process.returncode

                if return_code != 0:
                    err_msg = f"❌ Error in step {idx + 1} (Exit Code: {return_code}): {cmd_str}"
                    self.log_signal.emit(err_msg)
                    self.error_signal.emit(err_msg)
                    self.finished_signal.emit(False, self.task_name, err_msg)
                    return

                pct = int(((idx + 1) / total_steps) * 100)
                self.progress_signal.emit(pct)

            except Exception as e:
                err_msg = f"❌ Exception running command '{cmd_str}': {e}"
                self.log_signal.emit(err_msg)
                self.error_signal.emit(err_msg)
                self.finished_signal.emit(False, self.task_name, err_msg)
                return

        self.log_signal.emit("🎉 Task Completed Successfully!\n" + "─" * 60)
        self.finished_signal.emit(True, self.task_name, "Completed successfully.")

    def cancel(self):
        """Cancel process execution."""
        self._is_cancelled = True
        if self.process and self.process.poll() is None:
            try:
                self.process.terminate()
            except Exception:
                pass
