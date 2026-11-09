"""Dialog shown when the `spliceai` package isn't importable -- explains why
it isn't bundled (see LICENSE-THIRD-PARTY.md: source is PolyForm Strict, no
redistribution; models are CC BY-NC) and helps the user install it.

Setup normally installs SpliceAI already, so this only appears if that didn't
happen (e.g. the PC was offline during Setup), or when running from source. It
offers a one-click download of the pinned, hash-checked wheel from PyPI, or
installing that same file from disk -- no Python, pip or command line needed
(see spliceai_setup.install_spliceai_wheel).

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
    QDialog, QFileDialog, QHBoxLayout, QLabel, QPlainTextEdit, QProgressBar,
    QPushButton, QVBoxLayout,
)

from . import spliceai_setup as setup
from .i18n import t

EXPLANATION_HTML = t(
    "<p><b>SpliceAI is not installed yet.</b></p><p>Setup normally installs it; this occurs when the computer "
    "was not connected to the internet during Setup, or when the program is run from source.</p><p>SpliceAI is "
    "developed by Illumina and is not included in this program. <b>Download and install</b> downloads the "
    "official SpliceAI {version} package (about 16 MB) once from pypi.org and verifies that it is the published "
    "file. SpliceAI is licensed for academic and non-commercial use only (source code: PolyForm Strict 1.0.0, "
    "models: CC BY-NC 4.0).</p><p><b>Live scoring</b> requires SpliceAI; the rest of the program works without "
    "it.</p>",
    version=setup.SPLICEAI_VERSION)


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
        self.setWindowTitle(t("SpliceAI Setup Required"))
        self.resize(620, 460)
        self._worker = None
        self.installed = setup.is_spliceai_installed()

        layout = QVBoxLayout(self)

        explanation_label = QLabel(EXPLANATION_HTML)
        explanation_label.setTextFormat(Qt.RichText)
        explanation_label.setWordWrap(True)
        layout.addWidget(explanation_label)

        link = QLabel(f'<a href="{setup.SPLICEAI_REPO_URL}">{setup.SPLICEAI_REPO_URL}</a>')
        link.setOpenExternalLinks(True)
        layout.addWidget(link)

        install_row = QHBoxLayout()
        self.auto_button = QPushButton(t("Download and install"))
        self.auto_button.clicked.connect(self._on_auto_install)
        install_row.addWidget(self.auto_button)
        self.file_button = QPushButton(t("Install from file..."))
        self.file_button.setToolTip(t("Use a copy of {name}, e.g. from a USB drive.", name=setup.SPLICEAI_WHEEL_NAME))
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

        self.status_label = QLabel("")
        self.status_label.setWordWrap(True)
        layout.addWidget(self.status_label)
        self._refresh_status_label()

        button_row = QHBoxLayout()
        self.recheck_button = QPushButton(t("Re-check"))
        self.recheck_button.clicked.connect(self._on_recheck)
        button_row.addWidget(self.recheck_button)
        button_row.addStretch(1)
        self.decline_button = QPushButton(t("Continue without SpliceAI"))
        self.decline_button.clicked.connect(self.reject)
        button_row.addWidget(self.decline_button)
        self.done_button = QPushButton(t("Done"))
        self.done_button.clicked.connect(self.accept)
        self.done_button.setEnabled(self.installed)
        button_row.addWidget(self.done_button)
        layout.addLayout(button_row)

    def _refresh_status_label(self):
        if self.installed:
            self.status_label.setText(t("✓ SpliceAI is installed and importable."))
            self.status_label.setStyleSheet("color: #1a7f37;")
        else:
            self.status_label.setText(t("✗ SpliceAI is not currently importable."))
            self.status_label.setStyleSheet("color: #b00020;")

    def _on_recheck(self):
        self.installed = setup.refresh_after_install()
        self._refresh_status_label()
        self.done_button.setEnabled(self.installed)

    def _set_busy(self, busy):
        self.auto_button.setEnabled(not busy)
        self.file_button.setEnabled(not busy)
        self.recheck_button.setEnabled(not busy)
        self.progress_bar.setVisible(busy)

    def _start(self, task):
        self._set_busy(True)
        self.log_box.setVisible(True)
        self.log_box.clear()
        self.status_label.setText(t("Installing..."))
        self.status_label.setStyleSheet("")

        self._worker = InstallWorker(task)
        self._worker.line.connect(self.log_box.appendPlainText)
        self._worker.finished_ok.connect(self._on_install_finished)
        self._worker.failed.connect(self._on_install_failed)
        self._worker.start()

    def _on_auto_install(self):
        self._start(setup.install_spliceai_wheel)

    def _on_install_from_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, t("Select the SpliceAI package"), "", t("SpliceAI package ({name});;All files (*)", name=setup.SPLICEAI_WHEEL_NAME)
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
                t("The install reported success, but SpliceAI still isn't importable -- see the log above.")
            )
            self.status_label.setStyleSheet("color: #b00020;")

    def _on_install_failed(self, message):
        self._worker = None
        self._set_busy(False)
        self.status_label.setText(t("Install failed: {message}", message=message))
        self.status_label.setStyleSheet("color: #b00020;")
