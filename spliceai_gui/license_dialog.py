"""Full-text third-party license viewer -- opened on demand from
FirstLaunchDialog's "Read the full license..." button, not shown
automatically and not any kind of gate. Purely informational: a normal,
freely scrollable read-only text view of LICENSE-THIRD-PARTY.md (covering
SpliceAI, SnpEff, NCBI MANE, and the UCSC reference genome) with a single
Close button. Closing it (Close button, X, Alt+F4 -- all equivalent) just
dismisses this window and returns to whatever's underneath (the consent
dialog, or the main window if already past it); nothing else in the app
reads or cares how/whether this was ever opened.
"""
from PySide6.QtWidgets import QDialog, QHBoxLayout, QLabel, QPushButton, QTextEdit, QVBoxLayout

from .assets_paths import LICENSE_THIRD_PARTY_PATH

LICENSE_DIALOG_TITLE = "Third-Party Licenses -- SpliceAI Variant Scoring"

HEADER_HTML = (
    "<b>Third-Party Licenses</b><br>Covers the licensing terms for SpliceAI, SnpEff, NCBI "
    "MANE, and the UCSC reference genome -- none of which are bundled with this "
    "application."
)


def _load_license_text():
    """Returns LICENSE-THIRD-PARTY.md's raw Markdown text, or a clear
    placeholder explaining it couldn't be found -- opening this viewer must
    never itself crash the app, just degrade to an honest error message."""
    try:
        return LICENSE_THIRD_PARTY_PATH.read_text(encoding="utf-8")
    except OSError as exc:
        return (
            "# Third-Party Licenses\n\n**Could not load LICENSE-THIRD-PARTY.md** (expected "
            "at `"
            f"{LICENSE_THIRD_PARTY_PATH}`): {exc}"
            "\n\nThis file documents the licensing terms for SpliceAI, SnpEff, NCBI MANE, "
            "and the UCSC reference genome, none of which are bundled with this "
            "application -- see the project's repository for the current text."
        )


class LicenseDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(LICENSE_DIALOG_TITLE)
        self.setModal(True)
        self.resize(760, 680)

        layout = QVBoxLayout(self)

        header = QLabel(HEADER_HTML)
        header.setWordWrap(True)
        layout.addWidget(header)

        self.text_view = QTextEdit()
        self.text_view.setReadOnly(True)
        self.text_view.setMarkdown(_load_license_text())
        layout.addWidget(self.text_view, stretch=1)
        self.text_view.setFocus()

        button_row = QHBoxLayout()
        button_row.addStretch(1)
        self.close_button = QPushButton("Close")
        self.close_button.setDefault(True)
        self.close_button.clicked.connect(self.accept)
        button_row.addWidget(self.close_button)
        layout.addLayout(button_row)
