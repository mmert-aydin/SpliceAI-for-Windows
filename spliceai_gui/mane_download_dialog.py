"""Dialog for downloading the current NCBI MANE Select summary file --
mirrors ReferenceDownloadDialog/SnpEffDownloadDialog's background-thread +
progress pattern, but for a single small (~1MB) file with no byte-level
progress available: mane.fetch_mane_summary() does one plain response.read()
rather than a chunked loop (see its own docstring for why that's fine at
this size), so there's nothing to report a byte count for -- progress here
is a single indeterminate phase, like the "index"/"unzip" phases of the
FASTA/SnpEff dialogs.

Unlike the FASTA and SnpEff downloads, the destination folder is NOT
user-browsable here: mane_dir has no corresponding GUI field anywhere else
in this app (the pipeline always reads spliceai_pipeline.mane.DEFAULT_MANE_DIR
-- see cli.py/worker.py, neither ever passes a mane_dir override), so letting
someone pick a different folder would silently produce a file the pipeline
never looks at. Always downloads to DEFAULT_MANE_DIR, which is itself
frozen-aware (see app_paths.py) -- the packaged exe and a source checkout
each get their own correct default automatically.
"""
from PySide6.QtCore import QThread, Signal
from PySide6.QtWidgets import (
    QDialog, QHBoxLayout, QLabel, QMessageBox, QProgressBar, QPushButton,
    QVBoxLayout,
)
from spliceai_pipeline.mane import DEFAULT_MANE_DIR, ManeDownloadError, fetch_mane_summary


class ManeDownloadWorker(QThread):
    finished_ok = Signal(str)
    failed = Signal(str)

    def run(self):
        try:
            path = fetch_mane_summary(DEFAULT_MANE_DIR)
            self.finished_ok.emit(path)
        except ManeDownloadError as exc:
            self.failed.emit(str(exc))
        except Exception as exc:
            self.failed.emit(f"{type(exc).__name__}: {exc}")


class ManeDownloadDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Download MANE Select data")
        self.resize(520, 180)
        self._worker = None
        self._result_path = None

        layout = QVBoxLayout(self)

        note = QLabel(
            "Downloads the current NCBI MANE Select summary file (~1 MB) -- the exact "
            "filename is discovered from NCBI's own directory listing rather than "
            "hardcoded, so this always gets the latest version. Used to prefer the "
            "canonical transcript when a variant overlaps several during SnpEff "
            "annotation; optional otherwise."
        )
        note.setWordWrap(True)
        layout.addWidget(note)

        dest_row = QHBoxLayout()
        dest_row.addWidget(QLabel("Save to:"))
        dest_row.addWidget(QLabel(DEFAULT_MANE_DIR))
        dest_row.addStretch(1)
        layout.addLayout(dest_row)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 1)
        layout.addWidget(self.progress_bar)
        self.status_label = QLabel("")
        layout.addWidget(self.status_label)

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        self.cancel_button = QPushButton("Cancel")
        self.cancel_button.clicked.connect(self.reject)
        button_row.addWidget(self.cancel_button)
        self.download_button = QPushButton("Download")
        self.download_button.clicked.connect(self._on_download_clicked)
        button_row.addWidget(self.download_button)
        layout.addLayout(button_row)

    def _on_download_clicked(self):
        self.download_button.setEnabled(False)
        self.status_label.setText("Downloading...")
        self.progress_bar.setRange(0, 0)

        self._worker = ManeDownloadWorker(self)
        self._worker.finished_ok.connect(self._on_finished)
        self._worker.failed.connect(self._on_failed)
        self._worker.start()

    def _on_finished(self, path):
        self._result_path = path
        self._worker = None
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(1)
        self.status_label.setText(f"Done: {path}")
        QMessageBox.information(self, "MANE Select data ready", f"Downloaded:\n{path}")
        self.accept()

    def _on_failed(self, message):
        self._worker = None
        self.download_button.setEnabled(True)
        self.progress_bar.setRange(0, 1)
        self.progress_bar.setValue(0)
        self.status_label.setText("Failed.")
        QMessageBox.critical(self, "Download failed", message)

    def result_path(self):
        return self._result_path
