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
from PySide6.QtWidgets import (
    QApplication, QComboBox, QDialog, QFrame, QHBoxLayout, QLabel, QPushButton, QScrollArea, QVBoxLayout,
)

from . import config, i18n
from .i18n import N_, t
from .license_dialog import LicenseDialog

FIRST_LAUNCH_TITLE = N_("Welcome to SpliceAI Variant Scoring.")

FIRST_LAUNCH_SECTIONS = [
    (
        None,   # the purpose, right under the welcome line, needs no heading
        N_("This program was developed so that users who wish to work with SpliceAI can examine their variants without "
           "Linux, WSL or the command line, and so that SpliceAI's raw output can be brought into a form suitable for "
           "analysis."),
    ),
    (
        N_("About SpliceAI"),
        N_("SpliceAI is a deep-learning model developed by Illumina that predicts whether a genetic variant alters "
           "normal RNA splicing. This program provides an interface to the model; it does not modify or replace it."),
    ),
    (
        N_("About Sequence to cDNA"),
        N_("The Sequence to cDNA tab takes a sample sequence of one gene, or variants entered as coordinates, and "
           "uses SpliceAI to predict how it will be spliced. This makes it possible to compare the reference exons, "
           "the normal splicing predicted for the reference sequence, and the predicted effect of the variant."),
    ),
    (
        N_("Terms of Use"),
        N_("This software is provided for academic and research use only, without warranty of any kind, express or "
        "implied. It is not a diagnostic tool, and its output must not be used for clinical decision-making without "
        "independent verification by qualified personnel. SpliceAI itself is licensed separately by Illumina, Inc.; "
        "use of this program does not grant any rights under that license, and the user is solely responsible for "
        "compliance with it. The developer assumes no liability for any use of this software or reliance on its "
        "output. By using this program, the user accepts all associated risk."),
    ),
]


# Where to send bugs and suggestions -- shown just above the signature.
CONTACT_EMAIL = "aydinn.mmert@gmail.com"
CONTACT_LEAD = N_("Bug reports and suggestions")

# Shown right-aligned under the notice, in a handwriting-style font when one
# is installed (see FirstLaunchDialog.__init__).
SIGNATURE_TEXT = "— MMA"

# Handwriting-style fonts, tried in order for the signature line; see
# _pick_signature_font_family().


SIGNATURE_FONT_CANDIDATES = (
    "Segoe Script", "Segoe Print", "Lucida Handwriting", "Brush Script MT", "Monotype Corsiva",
)


def first_launch_html():
    """The notice in the window's language (i18n.t), built when shown: the language can be changed right here."""
    parts = [f"<p><b>{t(FIRST_LAUNCH_TITLE)}</b></p>"]
    for header, body in FIRST_LAUNCH_SECTIONS:
        parts.append(f"<p><b>{t(header)}</b><br>{t(body)}</p>" if header else f"<p>{t(body)}</p>")
    return "".join(parts)


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
    """The notice, with the language choice (English / Türkçe) at its top right: the first thing a new user
    sees is in a language they read, and the window after it opens in the same language."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setModal(True)
        self.resize(560, 480)

        layout = QVBoxLayout(self)
        top = QHBoxLayout()
        top.addStretch(1)
        self.language_combo = QComboBox()
        for code, name in i18n.LANGUAGES.items():
            self.language_combo.addItem(name, code)
        self.language_combo.setCurrentIndex(list(i18n.LANGUAGES).index(i18n.language()))
        self.language_combo.setToolTip("Language / Dil")
        self.language_combo.currentIndexChanged.connect(self._on_language_changed)
        top.addWidget(self.language_combo)
        layout.addLayout(top)

        self.body_label = QLabel()
        self.body_label.setTextFormat(Qt.RichText)
        self.body_label.setWordWrap(True)
        # In a scroll area: on a short screen (a laptop at 150 %) the text scrolls and the buttons stay in sight.
        self.body_scroll = QScrollArea()
        self.body_scroll.setWidgetResizable(True)
        self.body_scroll.setFrameShape(QFrame.NoFrame)
        self.body_scroll.setWidget(self.body_label)
        layout.addWidget(self.body_scroll, stretch=1)

        self.contact_label = QLabel()
        self.contact_label.setTextFormat(Qt.RichText)
        self.contact_label.setOpenExternalLinks(True)
        layout.addWidget(self.contact_label, alignment=Qt.AlignRight)

        self.signature_label = QLabel(SIGNATURE_TEXT)
        signature_family = _pick_signature_font_family()
        font = QFont(signature_family) if signature_family else QFont(self.signature_label.font())
        font.setItalic(True)
        font.setPointSize(13)
        self.signature_label.setFont(font)
        layout.addWidget(self.signature_label, alignment=Qt.AlignRight)

        full_license_row = QHBoxLayout()
        self.full_license_button = QPushButton()
        self.full_license_button.setFlat(True)
        self.full_license_button.setCursor(Qt.PointingHandCursor)
        self.full_license_button.setStyleSheet("QPushButton { color: palette(link); text-align: left; }")
        self.full_license_button.clicked.connect(self._on_read_full_license)
        full_license_row.addWidget(self.full_license_button)
        full_license_row.addStretch(1)
        layout.addLayout(full_license_row)

        button_row = QHBoxLayout()
        self.decline_button = QPushButton()
        self.decline_button.clicked.connect(self.reject)
        button_row.addWidget(self.decline_button)
        button_row.addStretch(1)
        self.accept_button = QPushButton()
        self.accept_button.setDefault(True)
        self.accept_button.clicked.connect(self.accept)
        button_row.addWidget(self.accept_button)
        layout.addLayout(button_row)
        self._show_texts()

    def _show_texts(self):
        self.setWindowTitle(t("Welcome to SpliceAI Variant Scoring"))
        self.body_label.setText(first_launch_html())
        self.contact_label.setText(f'{t(CONTACT_LEAD)}: <a href="mailto:{CONTACT_EMAIL}">{CONTACT_EMAIL}</a>')
        self.contact_label.setToolTip(t("Opens your e-mail program"))
        self.full_license_button.setText(t("Read the full license..."))
        self.decline_button.setText(t("Decline"))
        self.accept_button.setText(t("I have read and accept the terms"))
        if self.isVisible():
            self._fit_text()

    def showEvent(self, event):
        super().showEvent(event)
        self._fit_text()

    def _fit_text(self):
        """Wrapped text doesn't make the window taller by itself (the Turkish one is longer, and the language can
        change while it is open): as tall as the whole text needs, so it never scrolls on a screen it fits -- on a
        shorter one the window is the screen's height and the text scrolls (theme.fit_to_screen). Measured once
        shown, at the text's real width."""
        from .theme import fit_to_screen

        self.layout().activate()
        view = self.body_scroll.viewport()
        extra = self.body_label.heightForWidth(view.width()) - view.height()
        if extra > 0:
            self.resize(self.width(), self.height() + extra)
            fit_to_screen(self)

    def _on_language_changed(self):
        """Takes effect at once: nothing else of the program is built yet (main.py imports it after this)."""
        settings = config.load()
        settings["language"] = i18n._language = self.language_combo.currentData()
        config.save(settings)
        i18n.install_qt_translation(QApplication.instance())
        self._show_texts()

    def _on_read_full_license(self):
        LicenseDialog(self).exec()
