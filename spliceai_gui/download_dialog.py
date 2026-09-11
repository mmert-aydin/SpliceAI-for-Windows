"""Dialog for downloading a reference FASTA (UCSC goldenPath), decompressing
it, and building its .fai index -- runs in a background thread so the GUI
doesn't freeze during the ~1GB download.
"""
import os

from PySide6.QtCore import QThread, Signal
from PySide6.QtWidgets import (
    QComboBox, QDialog, QFileDialog, QHBoxLayout, QLabel, QLineEdit, QMessageBox,
    QProgressBar, QPushButton, QVBoxLayout,
)

from . import config
from .reference_download import REFERENCE_URLS, DownloadCancelled, fetch_reference

PHASE_LABELS = {"download": "Downloading", "decompress": "Decompressing", "index": "Building index"}


class ReferenceDownloadWorker(QThread):
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

            fasta_path = fetch_reference(
                self.build, self.dest_dir, on_progress=on_progress,
                should_cancel=lambda: self._cancel_requested,
            )
            self.finished_ok.emit(fasta_path)
        except DownloadCancelled:
            self.cancelled.emit()
        except Exception as exc:
            self.failed.emit(str(exc))


class ReferenceDownloadDialog(QDialog):
    def __init__(self, parent=None, default_build="hg19"):
        super().__init__(parent)
        self.setWindowTitle("Download reference FASTA")
        self.resize(560, 220)
        self._worker = None
        self._result_path = None

        layout = QVBoxLayout(self)

        build_row = QHBoxLayout()
        build_row.addWidget(QLabel("Build:"))
        self.build_combo = QComboBox()
        self.build_combo.addItems(list(REFERENCE_URLS.keys()))
        if default_build in REFERENCE_URLS:
            self.build_combo.setCurrentText(default_build)
        build_row.addWidget(self.build_combo)
        build_row.addStretch(1)
        layout.addLayout(build_row)

        dir_row = QHBoxLayout()
        dir_row.addWidget(QLabel("Save to folder:"))
        self.dir_edit = QLineEdit()
        self.dir_edit.setText(self._default_dest_dir())
        dir_row.addWidget(self.dir_edit, stretch=1)
        self.browse_button = QPushButton("Browse...")
        self.browse_button.clicked.connect(self._on_browse)
        dir_row.addWidget(self.browse_button)
        layout.addLayout(dir_row)

        note = QLabel(
            "Downloads ~1 GB (compressed) from UCSC's goldenPath, decompresses to ~3 GB, "
            "then builds the .fai index. This can take a while depending on your "
            "connection -- the progress bar below will track it."
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

    def _default_dest_dir(self):
        saved = config.load().get("last_reference_download_dir")
        if saved:
            return saved
        return os.path.join(os.path.expanduser("~"), "SpliceAI_reference_data")

    def _on_browse(self):
        path = QFileDialog.getExistingDirectory(self, "Select destination folder", self.dir_edit.text())
        if path:
            self.dir_edit.setText(path)

    def _set_inputs_enabled(self, enabled):
        self.build_combo.setEnabled(enabled)
        self.dir_edit.setEnabled(enabled)
        self.browse_button.setEnabled(enabled)
        self.download_button.setEnabled(enabled)

    def _on_download_clicked(self):
        dest_dir = self.dir_edit.text().strip()
        if not dest_dir:
            QMessageBox.warning(self, "No folder selected", "Choose a destination folder first.")
            return

        settings = config.load()
        settings["last_reference_download_dir"] = dest_dir
        config.save(settings)

        self._set_inputs_enabled(False)
        self.status_label.setText("Starting...")
        self.progress_bar.setRange(0, 0)

        build = self.build_combo.currentText()
        self._worker = ReferenceDownloadWorker(build, dest_dir)
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

    def _on_finished(self, fasta_path):
        self._result_path = fasta_path
        self._worker = None
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(1)
        self.status_label.setText(f"Done: {fasta_path}")
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

    def result_path(self):
        return self._result_path

    def result_build(self):
        """The build actually downloaded -- the user can change it in this
        dialog, so it isn't necessarily the one it was opened for."""
        return self.build_combo.currentText()
