"""The small "(?)" icon that opens an explanation, shared by every tab and
dialog (it used to be MainWindow._make_help_button).

Same pattern everywhere, on purpose: a hover-only tooltip was tried first for
the Mode icon and turned out to be unreliable, so every help icon has a real
click handler, and the tooltip only previews the text.
"""
from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtGui import Qt as GuiQt     # its mightBeRichText (not in QtCore's Qt)
from PySide6.QtWidgets import QDialog, QDialogButtonBox, QHBoxLayout, QLabel, QTextBrowser, QToolButton, QVBoxLayout

from .assets_paths import LOGO_PNG

ABOUT_DIALOG_LOGO_SIZE = 72


def make_help_button(parent, dialog_title, help_text, font_size="10pt", logo=False, label=None):
    """A "(?)" button that pops up help_text (plain or rich text). logo=True
    shows the app logo at the top of the popup; label puts a word before the
    "(?)", for several helps side by side."""
    button = QToolButton(parent)
    button.setText(f"{label} (?)" if label else "(?)")
    button.setAutoRaise(True)
    button.setCursor(Qt.PointingHandCursor)
    button.setStyleSheet(
        f"QToolButton {{ color: palette(link); font-weight: bold; font-size: {font_size}; border: none; }}"
    )
    button.setToolTip(help_text)
    button.clicked.connect(lambda: help_dialog(parent, dialog_title, help_text, logo=logo).exec())
    return button


class _HelpDialog(QDialog):
    def showEvent(self, event):
        """Once shown, at the text's real width: as tall as the text, so it never scrolls on a screen it fits."""
        from .theme import fit_to_screen

        super().showEvent(event)
        self.layout().activate()
        view = self.body.viewport()
        extra = int(self.body.document().size().height()) + 2 - view.height()
        if extra > 0:
            self.resize(self.width(), self.height() + extra)
            fit_to_screen(self)


def help_dialog(parent, title, text, logo=False, level=None):
    """A help (or notice) text in a window: message-box wide, as tall as the text, and scrolling when the text is
    taller than the screen (theme.fit_to_screen trims the window to it) -- a QMessageBox can't scroll, so a long
    help ran off small screens."""
    from . import theme

    d = _HelpDialog(parent)
    d.setWindowTitle(title)
    v = QVBoxLayout(d)
    row = QHBoxLayout()
    if logo and LOGO_PNG.exists():
        icon = QLabel()
        icon.setPixmap(QPixmap(str(LOGO_PNG)).scaled(ABOUT_DIALOG_LOGO_SIZE, ABOUT_DIALOG_LOGO_SIZE,
                                                      Qt.KeepAspectRatio, Qt.SmoothTransformation))
        row.addWidget(icon, alignment=Qt.AlignTop)
    body = QTextBrowser()
    body.setOpenExternalLinks(True)
    body.setFrameShape(QTextBrowser.NoFrame)
    body.setStyleSheet(f"background: transparent; color: {theme.CRIT if level == 'critical' else theme.INK};")
    if GuiQt.mightBeRichText(text):     # as QMessageBox decided it
        body.setHtml(text)
    else:
        body.setPlainText(text)
    row.addWidget(body, stretch=1)
    d.body = body
    v.addLayout(row)
    buttons = QDialogButtonBox(QDialogButtonBox.Ok)
    buttons.accepted.connect(d.accept)
    v.addWidget(buttons)
    width = 560
    body.document().setTextWidth(width - (ABOUT_DIALOG_LOGO_SIZE + 12 if logo else 0) - 24)
    d.resize(width, int(body.document().size().height()) + 80)
    return d


def make_notice_button(parent, level, title):
    """Spec 15 item 5: a "?" (cautions) or red "!" (critical problems) icon; its
    tooltip is the notices' text, and a click shows them in full."""
    from . import theme

    button = QToolButton(parent)
    button.setText("!" if level == "critical" else "?")
    button.setAutoRaise(True)
    button.setCursor(Qt.PointingHandCursor)
    color = theme.CRIT if level == "critical" else theme.TEAL
    button.setStyleSheet(f"QToolButton {{ color: white; background: {color}; font-weight: bold; font-size: 10pt; "
                         "border: none; border-radius: 9px; min-width: 18px; min-height: 18px; }")
    button.clicked.connect(lambda: help_dialog(parent, title, button.toolTip(), level=level).exec())
    button.setVisible(False)
    return button


def notice_icon(level):
    """A small round "?" or "!" for a table cell (the results' notices)."""
    from PySide6.QtGui import QColor, QFont, QIcon, QPainter

    from . import theme

    pixmap = QPixmap(16, 16)
    pixmap.fill(Qt.transparent)
    p = QPainter(pixmap)
    p.setRenderHint(QPainter.Antialiasing)
    p.setBrush(QColor(theme.CRIT if level == "critical" else theme.TEAL))
    p.setPen(Qt.NoPen)
    p.drawEllipse(1, 1, 14, 14)
    font = QFont()
    font.setBold(True)
    font.setPointSize(8)
    p.setFont(font)
    p.setPen(QColor("white"))
    p.drawText(pixmap.rect(), Qt.AlignCenter, "!" if level == "critical" else "?")
    p.end()
    return QIcon(pixmap)
