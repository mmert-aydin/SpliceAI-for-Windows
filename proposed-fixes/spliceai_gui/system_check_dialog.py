"""First-run "can this machine actually run an analysis?" dialog.

Shown by main() before the main window on first launch, and again on any
later launch where a REQUIRED check fails (see config's "setup_complete").
Renders preflight.run_checks() as a list of rows -- a status glyph, a title,
a one-line detail, and (where we can do something about it) a Fix button
that opens the relevant existing downloader/setup dialog. "Continue" stays
disabled while any REQUIRED check is still failing; "Re-check" re-runs
everything after a fix.

No pipeline or check logic lives here -- this file is purely Qt plumbing
around preflight.py, the same split first_launch_dialog.py / license_dialog.py
use.
"""

from __future__ import annotations

import webbrowser

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget,
)

from spliceai_gui import config
from spliceai_gui.preflight import (
    CheckId,
    Severity,
    Status,
    has_blocking_failure,
    run_checks,
    visible_results,
)

_GLYPH = {
    Status.OK: ("✓", "#2e7d32"),
    Status.MISSING: ("!", "#ed6c02"),
    Status.ERROR: ("✕", "#c62828"),
}
# A blocking REQUIRED failure gets the red glyph regardless of MISSING/ERROR.
_BLOCKING_COLOR = "#c62828"

# Adoptium Temurin is the no-registration, GPL+CE JDK/JRE SnpEff is happy with.
_JAVA_HELP_URL = "https://adoptium.net/temurin/releases/?version=21"
_JAVA_WINGET = "winget install --id EclipseAdoptium.Temurin.21.JRE -e"


class _CheckRow(QFrame):
    def __init__(self, result, on_fix):
        super().__init__()
        self.setFrameShape(QFrame.StyledPanel)
        row = QHBoxLayout(self)

        glyph, color = _GLYPH[result.status]
        if result.blocking:
            color = _BLOCKING_COLOR
        badge = QLabel(glyph)
        badge.setStyleSheet(f"color: {color}; font-weight: bold; font-size: 15px;")
        badge.setFixedWidth(18)
        badge.setAlignment(Qt.AlignTop)
        row.addWidget(badge)

        text_col = QVBoxLayout()
        title = QLabel(result.title)
        title.setStyleSheet("font-weight: bold;")
        text_col.addWidget(title)
        if result.detail:
            detail = QLabel(result.detail)
            detail.setWordWrap(True)
            detail.setStyleSheet("color: palette(mid);")
            text_col.addWidget(detail)
        if result.severity is Severity.OPTIONAL and result.status is not Status.OK:
            tag = QLabel("optional")
            tag.setStyleSheet("color: palette(mid); font-style: italic;")
            text_col.addWidget(tag)
        row.addLayout(text_col, 1)

        if result.status is not Status.OK and result.fixable:
            fix_button = QPushButton(result.fix_label)
            fix_button.clicked.connect(lambda: on_fix(result))
            fix_button.setFixedHeight(28)
            row.addWidget(fix_button, 0, Qt.AlignTop)


class SystemCheckDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("System check")
        self.setModal(True)
        self.resize(620, 520)

        self._results = []

        outer = QVBoxLayout(self)

        heading = QLabel(
            "<b>Checking this computer for everything an analysis needs.</b><br>"
            "Anything marked in orange or red can be set up below before you "
            "start. Optional items can be left for later."
        )
        heading.setWordWrap(True)
        heading.setTextFormat(Qt.RichText)
        outer.addWidget(heading)

        self._scroll = QScrollArea()
        self._scroll.setWidgetResizable(True)
        outer.addWidget(self._scroll, 1)

        button_row = QHBoxLayout()
        self.recheck_button = QPushButton("Re-check")
        self.recheck_button.clicked.connect(self.refresh)
        button_row.addWidget(self.recheck_button)
        button_row.addStretch(1)

        self.quit_button = QPushButton("Quit")
        self.quit_button.clicked.connect(self.reject)
        button_row.addWidget(self.quit_button)

        self.continue_button = QPushButton("Continue")
        self.continue_button.setDefault(True)
        self.continue_button.clicked.connect(self._on_continue)
        button_row.addWidget(self.continue_button)
        outer.addLayout(button_row)

        self.refresh()

    # -- rendering ---------------------------------------------------------- #

    def refresh(self):
        self._results = run_checks(config.load())
        container = QWidget()
        col = QVBoxLayout(container)
        for result in visible_results(self._results):
            col.addWidget(_CheckRow(result, self._on_fix))
        col.addStretch(1)
        self._scroll.setWidget(container)

        blocking = has_blocking_failure(self._results)
        self.continue_button.setEnabled(not blocking)
        self.continue_button.setToolTip(
            "Resolve the items marked in red first." if blocking else ""
        )

    # -- actions ----------------------------------------------------------- #

    def _on_continue(self):
        settings = config.load()
        settings["setup_complete"] = not has_blocking_failure(self._results)
        config.save(settings)
        self.accept()

    def _on_fix(self, result):
        handler = {
            CheckId.SPLICEAI: self._fix_spliceai,
            CheckId.REFERENCE_FASTA: self._fix_reference,
            CheckId.JAVA: self._fix_java,
            CheckId.SNPEFF: self._fix_snpeff,
            CheckId.MANE: self._fix_mane,
        }.get(result.id)
        if handler is None:
            return
        try:
            handler(result)
        except Exception as exc:  # noqa: BLE001 -- a failed fix must not kill the dialog
            QMessageBox.critical(self, "Setup failed", str(exc))
        self.refresh()

    def _fix_spliceai(self, _result):
        from spliceai_gui.spliceai_setup_dialog import SpliceAISetupDialog

        SpliceAISetupDialog(self).exec()

    def _fix_reference(self, result):
        from spliceai_gui.download_dialog import ReferenceDownloadDialog

        build = result.context.get("build", "hg19")
        dlg = ReferenceDownloadDialog(self, build)
        if dlg.exec() == QDialog.Accepted and dlg.result_path():
            settings = config.load()
            settings["fasta_path"] = dlg.result_path()
            config.save(settings)

    def _fix_java(self, _result):
        box = QMessageBox(self)
        box.setIcon(QMessageBox.Information)
        box.setWindowTitle("Install Java")
        box.setText(
            "SnpEff needs a Java 21+ runtime, which cannot be bundled with this "
            "app. The quickest way to get one:\n\n"
            f"    {_JAVA_WINGET}\n\n"
            "Run that in a terminal, or download an installer from the Adoptium "
            "page. Then click Re-check."
        )
        open_page = box.addButton("Open Adoptium page", QMessageBox.AcceptRole)
        box.addButton(QMessageBox.Close)
        box.exec()
        if box.clickedButton() is open_page:
            webbrowser.open(_JAVA_HELP_URL)

    def _fix_snpeff(self, result):
        from spliceai_gui.snpeff_download_dialog import SnpEffDownloadDialog

        build = result.context.get("build", "hg19")
        dest_dir = result.context.get("snpeff_dir") or ""
        dlg = SnpEffDownloadDialog(self, build, dest_dir)
        if dlg.exec() == QDialog.Accepted and dest_dir:
            settings = config.load()
            settings["snpeff_dir"] = dest_dir
            config.save(settings)

    def _fix_mane(self, _result):
        from spliceai_gui.mane_download_dialog import ManeDownloadDialog

        ManeDownloadDialog(self).exec()


def run_system_check(app) -> bool:
    """Show the check modally (keeping the app alive while no main window
    exists). Returns True to proceed to the main window, False to quit."""
    previous = app.quitOnLastWindowClosed()
    app.setQuitOnLastWindowClosed(False)
    try:
        return SystemCheckDialog().exec() == QDialog.Accepted
    finally:
        app.setQuitOnLastWindowClosed(previous)


def should_show_on_startup() -> bool:
    """First launch (no config yet, or setup never completed) -> always show.
    Later launches -> only when a REQUIRED check is currently failing."""
    settings = config.load()
    if not settings.get("setup_complete"):
        return True
    return has_blocking_failure(run_checks(settings))
