"""Dialog for downloading SnpEff (core distribution + RefSeq database for the
selected build) -- runs in a background thread so the GUI doesn't freeze,
same pattern as ReferenceDownloadDialog for the FASTA.
"""
from PySide6.QtCore import QThread, Signal
from PySide6.QtWidgets import (
    QDialog, QFileDialog, QHBoxLayout, QLabel, QLineEdit, QMessageBox, QProgressBar,
    QPushButton, QVBoxLayout,
)

from spliceai_pipeline.snpeff import DEFAULT_SNPEFF_DIR

from .reference_download import DownloadCancelled
from .snpeff_download import fetch_snpeff

PHASE_LABELS = {
    "java": "Downloading Java (SnpEff needs it)",
    "java_unzip": "Unpacking Java",
    "download": "Downloading SnpEff",
    "unzip": "Unzipping",
    "database": "Downloading database",
}


class SnpEffDownloadWorker(QThread):
    progress = Signal(str, object, object)
    finished_ok = Signal(str)
    failed = Signal(str)
    cancelled = Signal()

    def __init__(self, build, dest_dir, parent=None):
        super().__init__(parent)
        self.build = build
        self.dest_dir = dest_dir
        self._cancel_requested = False

    def cancel(self):
        self._cancel_requested = True

    def run(self):
        try:
            def on_progress(phase, done, total):
                self.progress.emit(phase, done, total)

            jar_path = fetch_snpeff(
                self.dest_dir, self.build, on_progress=on_progress,
                should_cancel=lambda: self._cancel_requested,
            )
            self.finished_ok.emit(jar_path)
        except DownloadCancelled:
            self.cancelled.emit()
        except Exception as exc:
            self.failed.emit(str(exc))


class SnpEffDownloadDialog(QDialog):
    def __init__(self, parent=None, build="hg19", dest_dir=None):
        super().__init__(parent)
        self.setWindowTitle("Download SnpEff")
        self.resize(560, 220)
        self._worker = None
        self._result_jar_path = None
        self.build = build

        layout = QVBoxLayout(self)

        dir_row = QHBoxLayout()
        dir_row.addWidget(QLabel("Install to folder:"))
        self.dir_edit = QLineEdit()
        self.dir_edit.setText(dest_dir or DEFAULT_SNPEFF_DIR)
        dir_row.addWidget(self.dir_edit, stretch=1)
        self.browse_button = QPushButton("Browse...")
        self.browse_button.clicked.connect(self._on_browse)
        dir_row.addWidget(self.browse_button)
        layout.addLayout(dir_row)

        note = QLabel(
            f"Downloads SnpEff (~65 MB) plus the RefSeq-transcript database for {build}"
            " (a few hundred MB), and installs both -- this can take a few minutes "
            "depending on your connection. SnpEff runs on Java 21+: if this PC doesn't have it, "
            "a private copy of Java (~50 MB) is downloaded automatically too -- nothing to "
            "install by hand."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 1)
        layout.addWidget(self.progress_bar)
        self.status_label = QLabel("")
        layout.addWidget(self.status_label)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self._on_cancel_clicked)
        button_row.addWidget(self.cancel_button)
        self.download_button = QPushButton("Download")
        self.download_button.clicked.connect(self._on_download_clicked)
        button_row.addWidget(self.download_button)
        layout.addLayout(button_row)

    def _on_browse(self):
        path = QFileDialog.getExistingDirectory(self, "Select destination folder", self.dir_edit.text())
        if path:
            self.dir_edit.setText(path)

    def _set_inputs_enabled(self, enabled):
        self.dir_edit.setEnabled(enabled)
        self.browse_button.setEnabled(enabled)
        self.download_button.setEnabled(enabled)

    def _on_download_clicked(self):
        dest_dir = self.dir_edit.text().strip()
        if not dest_dir:
            QMessageBox.warning(self, "No folder selected", "Choose a destination folder first.")
            return

        self._set_inputs_enabled(False)
        self.status_label.setText("Starting...")
        self.progress_bar.setRange(0, 0)

        self._worker = SnpEffDownloadWorker(self.build, dest_dir)
        self._worker.progress.connect(self._on_progress)
        self._worker.finished_ok.connect(self._on_finished)
        self._worker.failed.connect(self._on_failed)
        self._worker.cancelled.connect(self._on_cancelled)
        self._worker.start()

    def _on_cancel_clicked(self):
        if self._worker is not None:
            self._worker.cancel()
            self.status_label.setText("Cancelling...")
            self.cancel_button.setEnabled(False)
        else:
            self.reject()

    def _on_progress(self, phase, done, total):
        label = PHASE_LABELS.get(phase, phase)
        if total:
            self.progress_bar.setRange(0, total)
            self.progress_bar.setValue(done or 0)
            pct = (done or 0) / total * 100
            self.status_label.setText(f"{label}... {pct:.0f}% ({(done or 0) / 1e6:.0f} MB / {total / 1e6:.0f} MB)")
        else:
            self.progress_bar.setRange(0, 0)
            self.status_label.setText(f"{label}...")

    def _on_finished(self, jar_path):
        self._result_jar_path = jar_path
        self._worker = None
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(1)
        self.status_label.setText(f"Done: {jar_path}")
        QMessageBox.information(self, "SnpEff ready", f"SnpEff and the {self.build} database are installed:\n{jar_path}")
        self.accept()

    def _on_failed(self, message):
        self._worker = None
        self._set_inputs_enabled(True)
        self.cancel_button.setEnabled(True)
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(0)
        self.status_label.setText("Failed.")
        QMessageBox.critical(self, "Download failed", message)

    def _on_cancelled(self):
        self._worker = None
        self._set_inputs_enabled(True)
        self.cancel_button.setEnabled(True)
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(0)
        self.status_label.setText("Cancelled.")

    def result_jar_path(self):
        return self._result_jar_path
