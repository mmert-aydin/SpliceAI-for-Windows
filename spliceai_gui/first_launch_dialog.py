"""Mandatory notice shown once per app launch, before the main window opens.

Not a "don't show again" dialog -- it's meant to run every time the process
starts (a startup notice, not a one-time onboarding step), so there is
deliberately no persisted "seen" flag here. It only needs to never reappear
within the same running session, which falls out naturally from main.py only
ever constructing it once, before the main window exists.

Two distinct outcomes, both exercised by tests -- not just the accept path:
- Accept ("I have read and accept the terms"): dialog closes, main() proceeds
  to build and show MainWindow as normal.
- Decline: dialog closes, but main() returns immediately without ever
  constructing MainWindow -- the process exits with no main window shown.

Also has a "Read the full license..." button that opens LicenseDialog (the
complete LICENSE-THIRD-PARTY.md, freely scrollable, no gate) on demand --
purely optional, never shown automatically, and doesn't affect Accept/
Decline above in any way; see license_dialog.py.
"""
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QFontDatabase
from PySide6.QtWidgets import QDialog, QHBoxLayout, QLabel, QPushButton, QVBoxLayout

from .license_dialog import LicenseDialog

FIRST_LAUNCH_TITLE = "Welcome to SpliceAI Variant Scoring."

FIRST_LAUNCH_SECTIONS = [
    (
        "Purpose of Development",
        "This tool was built to make SpliceAI usable as a standard local Windows application, without requiring "
        "Linux, WSL, or manual command-line operation. It grew out of practical need during variant analysis work, "
        "where the setup overhead between a VCF file and a scored result was the actual bottleneck, not the science "
        "itself.",
    ),
    (
        "About SpliceAI",
        "SpliceAI is a neural network model, developed by Illumina, that predicts the likelihood a genetic variant "
        "disrupts normal RNA splicing. This program is an interface to that model; it does not alter or replace it.",
    ),
    (
        "About SnpEff",
        "This tool also uses SnpEff (open source, MIT license) to annotate each variant with its affected transcript "
        "and coding-sequence position. Credit to its original author, Pablo Cingolani.",
    ),
    (
        "Terms of Use",
        "This software is provided for academic and research use only, without warranty of any kind, express or "
        "implied. It is not a diagnostic tool, and its output must not be used for clinical decision-making without "
        "independent verification by qualified personnel. SpliceAI itself is licensed separately by Illumina, Inc.; "
        "use of this program does not grant any rights under that license, and the user is solely responsible for "
        "compliance with it. The developer assumes no liability for any use of this software or reliance on its "
        "output. By using this program, the user accepts all associated risk.",
    ),
]


# Where to send bugs and suggestions -- shown just above the signature.
CONTACT_EMAIL = "aydinn.mmert@gmail.com"
CONTACT_LEAD = "Bugs, fixes and suggestions"
CONTACT_TEXT = f"{CONTACT_LEAD}: {CONTACT_EMAIL}"

# Shown right-aligned under the notice, in a handwriting-style font when one
# is installed (see FirstLaunchDialog.__init__).
SIGNATURE_TEXT = "— MMA"

# Handwriting-style fonts, tried in order for the signature line; see
# _pick_signature_font_family().


SIGNATURE_FONT_CANDIDATES = (
    "Segoe Script", "Segoe Print", "Lucida Handwriting", "Brush Script MT", "Monotype Corsiva",
)


def _build_first_launch_html():
    parts = [f"<p><b>{FIRST_LAUNCH_TITLE}</b></p>"]
    for header, body in FIRST_LAUNCH_SECTIONS:
        parts.append(f"<p><b>{header}</b><br>{body}</p>")
    return "".join(parts)


FIRST_LAUNCH_HTML = _build_first_launch_html()

# The same notice as plain text (title and sections separated by blank lines).


FIRST_LAUNCH_TEXT = FIRST_LAUNCH_TITLE + "\n\n" + "\n\n".join(
    f"{header}\n{body}" for header, body in FIRST_LAUNCH_SECTIONS
)


def _pick_signature_font_family():
    """Returns the first available candidate script/cursive font actually
    installed on this machine, or None if none of them are -- callers must
    handle None by falling back to the default font (still italicized)."""
    available = set(QFontDatabase.families())
    for candidate in SIGNATURE_FONT_CANDIDATES:
        if candidate in available:
            return candidate
    return None


class FirstLaunchDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Welcome to SpliceAI Variant Scoring")
        self.setModal(True)
        self.resize(560, 460)

        layout = QVBoxLayout(self)

        body_label = QLabel(FIRST_LAUNCH_HTML)
        body_label.setTextFormat(Qt.RichText)
        body_label.setWordWrap(True)
        layout.addWidget(body_label)

        self.contact_label = QLabel(
            f'{CONTACT_LEAD}: <a href="mailto:{CONTACT_EMAIL}">{CONTACT_EMAIL}</a>'
        )
        self.contact_label.setTextFormat(Qt.RichText)
        self.contact_label.setOpenExternalLinks(True)
        self.contact_label.setToolTip("Opens your e-mail program")
        layout.addWidget(self.contact_label, alignment=Qt.AlignRight)

        self.signature_label = QLabel(SIGNATURE_TEXT)
        signature_family = _pick_signature_font_family()
        font = QFont(signature_family) if signature_family else QFont(self.signature_label.font())
        font.setItalic(True)
        font.setPointSize(13)
        self.signature_label.setFont(font)
        layout.addWidget(self.signature_label, alignment=Qt.AlignRight)

        full_license_row = QHBoxLayout()
        self.full_license_button = QPushButton("Read the full license...")
        self.full_license_button.setFlat(True)
        self.full_license_button.setCursor(Qt.PointingHandCursor)
        self.full_license_button.setStyleSheet("QPushButton { color: palette(link); text-align: left; }")
        self.full_license_button.clicked.connect(self._on_read_full_license)
        full_license_row.addWidget(self.full_license_button)
        full_license_row.addStretch(1)
        layout.addLayout(full_license_row)

        button_row = QHBoxLayout()
        self.decline_button = QPushButton("Decline")
        self.decline_button.clicked.connect(self.reject)
        button_row.addWidget(self.decline_button)

        button_row.addStretch(1)

        self.accept_button = QPushButton("I have read and accept the terms")
        self.accept_button.setDefault(True)
        self.accept_button.clicked.connect(self.accept)
        button_row.addWidget(self.accept_button)

        layout.addLayout(button_row)

    def _on_read_full_license(self):
        LicenseDialog(self).exec()
