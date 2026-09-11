"""Dialog shown when the `spliceai` package isn't importable -- explains why
it isn't bundled (see LICENSE-THIRD-PARTY.md: source is PolyForm Strict, no
redistribution; models are CC BY-NC) and helps the user install it.

Installed app (frozen): Setup normally installs SpliceAI already, so this only
appears if that didn't happen (e.g. the PC was offline during Setup). It
offers a one-click download of the pinned, hash-checked wheel from PyPI, or
installing that same file from disk -- no Python, pip or command line needed
(see spliceai_setup.install_spliceai_wheel). From source: the original pip
flow, with a copyable manual command.

Not a hard gate -- SpliceAI is only required for LIVE scoring; a user working
entirely from precomputed data, or just downloading reference files, can
still use the rest of the app. So this dialog always offers a way to
continue without installing right now ("Continue without SpliceAI"), same as
"Done" once it's actually detected -- both close the dialog and let startup
proceed; only the persistent status indicator in the main window (see
main_window.py) and the specific error a user hits if they actually try to
live-score without it (worker.py) differ based on which one they picked.
"""
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QApplication, QDialog, QFileDialog, QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit, QProgressBar,
    QPushButton, QVBoxLayout,
)

from . import spliceai_setup as setup

EXPLANATION_HTML = (
    "<p><b>SpliceAI isn't installed.</b></p><p>This program is an interface to SpliceAI, "
    "but does not include it -- SpliceAI's source is licensed under PolyForm Strict (no "
    "redistribution) and its trained models under CC BY-NC 4.0, so each user obtains it "
    "directly from Illumina under their own acceptance of those terms, rather than "
    "receiving it bundled here. See LICENSE-THIRD-PARTY.md for the full "
    "picture.</p><p><b>Live scoring</b> (any variant not already covered by precomputed "
    "data) needs it installed. The rest of the app -- loading/checking VCFs, downloading "
    "the reference FASTA/SnpEff/MANE data -- works without it.</p>"
)

FROZEN_EXPLANATION_HTML = (
    "<p><b>SpliceAI isn't installed yet.</b></p>"
    "<p>Setup normally installs it; this happens when the computer wasn't connected to the "
    "internet during Setup.</p>"
    f"<p>SpliceAI is made by Illumina and isn't included in this program. <b>Download and "
    f"install</b> fetches the official SpliceAI {setup.SPLICEAI_VERSION} package (about 16 MB) "
    "once from pypi.org and checks that it's the published file. SpliceAI is licensed for "
    "academic and non-commercial use only (source: PolyForm Strict 1.0.0, models: "
    "CC BY-NC 4.0).</p>"
    "<p><b>Live scoring</b> needs it. The rest of the app works without it.</p>"
)


# Runs the install off the GUI thread, streaming each output line back
# through `line` so the dialog can show live progress.
class InstallWorker(QThread):
    line = Signal(str)
    finished_ok = Signal(bool)
    failed = Signal(str)

    def __init__(self, task, parent=None):
        super().__init__(parent)
        self.task = task

    def run(self):
        try:
            self.task(on_line=self.line.emit)
        except Exception as exc:
            self.failed.emit(f"{type(exc).__name__}: {exc}")
            return
        self.finished_ok.emit(setup.refresh_after_install())


class SpliceAISetupDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        frozen = setup.is_frozen()
        self.setWindowTitle("SpliceAI Setup Required")
        self.resize(620, 460)
        self._worker = None
        self.installed = setup.is_spliceai_installed()

        layout = QVBoxLayout(self)

        explanation_label = QLabel(FROZEN_EXPLANATION_HTML if frozen else EXPLANATION_HTML)
        explanation_label.setTextFormat(Qt.RichText)
        explanation_label.setWordWrap(True)
        layout.addWidget(explanation_label)

        link = QLabel(f'<a href="{setup.SPLICEAI_REPO_URL}">{setup.SPLICEAI_REPO_URL}</a>')
        link.setOpenExternalLinks(True)
        layout.addWidget(link)

        self.command_edit = None
        if not frozen:
            layout.addWidget(QLabel('Manual install -- run this command once, then click "Re-check" below:'))
            cmd_row = QHBoxLayout()
            self.command_edit = QLineEdit(setup.display_install_command())
            self.command_edit.setReadOnly(True)
            cmd_row.addWidget(self.command_edit, stretch=1)
            self.copy_button = QPushButton("Copy")
            self.copy_button.clicked.connect(self._on_copy)
            cmd_row.addWidget(self.copy_button)
            layout.addLayout(cmd_row)

        self._auto_task = setup.auto_install_task()
        self.auto_button = None
        self.file_button = None
        self.progress_bar = None
        self.log_box = None
        if self._auto_task:
            install_row = QHBoxLayout()
            self.auto_button = QPushButton("Download and install" if frozen else "Install automatically")
            self.auto_button.clicked.connect(self._on_auto_install)
            install_row.addWidget(self.auto_button)
            if frozen:
                self.file_button = QPushButton("Install from file...")
                self.file_button.setToolTip(f"Use a copy of {setup.SPLICEAI_WHEEL_NAME}, e.g. from a USB drive.")
                self.file_button.clicked.connect(self._on_install_from_file)
                install_row.addWidget(self.file_button)
            install_row.addStretch(1)
            layout.addLayout(install_row)

            self.progress_bar = QProgressBar()
            self.progress_bar.setRange(0, 0)
            self.progress_bar.setVisible(False)
            layout.addWidget(self.progress_bar)

            self.log_box = QPlainTextEdit()
            self.log_box.setReadOnly(True)
            self.log_box.setMaximumHeight(140)
            self.log_box.setVisible(False)
            layout.addWidget(self.log_box)
        else:
            no_python_note = QLabel(
                "No usable Python installation was found on this system for automatic "
                "install -- use the manual command above."
            )
            no_python_note.setWordWrap(True)
            layout.addWidget(no_python_note)

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)
        self._refresh_status_label()

        button_row = QHBoxLayout()
        self.recheck_button = QPushButton("Re-check")
        self.recheck_button.clicked.connect(self._on_recheck)
        button_row.addWidget(self.recheck_button)
        button_row.addStretch(1)
        self.decline_button = QPushButton("Continue without SpliceAI")
        self.decline_button.clicked.connect(self.reject)
        button_row.addWidget(self.decline_button)
        self.done_button = QPushButton("Done")
        self.done_button.clicked.connect(self.accept)
        self.done_button.setEnabled(self.installed)
        button_row.addWidget(self.done_button)
        layout.addLayout(button_row)

    def _refresh_status_label(self):
        if self.installed:
            self.status_label.setText("✓ SpliceAI is installed and importable.")
            self.status_label.setStyleSheet("color: #1a7f37;")
        else:
            self.status_label.setText("✗ SpliceAI is not currently importable.")
            self.status_label.setStyleSheet("color: #b00020;")

    def _on_copy(self):
        QApplication.clipboard().setText(self.command_edit.text())

    def _on_recheck(self):
        self.installed = setup.refresh_after_install()
        self._refresh_status_label()
        self.done_button.setEnabled(self.installed)

    def _set_busy(self, busy):
        self.auto_button.setEnabled(not busy)
        if self.file_button:
            self.file_button.setEnabled(not busy)
        self.recheck_button.setEnabled(not busy)
        self.progress_bar.setVisible(busy)

    def _start(self, task):
        self._set_busy(True)
        self.log_box.setVisible(True)
        self.log_box.clear()
        self.status_label.setText("Installing...")
        self.status_label.setStyleSheet("")

        self._worker = InstallWorker(task)
        self._worker.line.connect(self.log_box.appendPlainText)
        self._worker.finished_ok.connect(self._on_install_finished)
        self._worker.failed.connect(self._on_install_failed)
        self._worker.start()

    def _on_auto_install(self):
        self._start(self._auto_task)

    def _on_install_from_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Select the SpliceAI package", "", f"SpliceAI package ({setup.SPLICEAI_WHEEL_NAME});;All files (*)"
        )
        if path:
            self._start(lambda on_line=None: setup.install_spliceai_wheel(on_line=on_line, wheel_path=path))

    def _on_install_finished(self, installed):
        self._worker = None
        self._set_busy(False)
        self.installed = installed
        self._refresh_status_label()
        self.done_button.setEnabled(self.installed)
        if not installed:
            self.status_label.setText(
                "The install reported success, but SpliceAI still isn't importable -- see the log above."
            )
            self.status_label.setStyleSheet("color: #b00020;")

    def _on_install_failed(self, message):
        self._worker = None
        self._set_busy(False)
        self.status_label.setText(f"Install failed: {message}")
        self.status_label.setStyleSheet("color: #b00020;")
