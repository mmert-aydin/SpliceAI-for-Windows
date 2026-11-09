"""The small "(?)" icon that opens an explanation, shared by every tab and
dialog (it used to be MainWindow._make_help_button).

Same pattern everywhere, on purpose: a hover-only tooltip was tried first for
the Mode icon and turned out to be unreliable, so every help icon has a real
click handler, and the tooltip only previews the text.
"""
from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import QMessageBox, QToolButton

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
    if logo:
        button.clicked.connect(lambda: show_help_with_logo(parent, dialog_title, help_text))
    else:
        button.clicked.connect(lambda: QMessageBox.information(parent, dialog_title, help_text))
    return button


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
    button.clicked.connect(lambda: (QMessageBox.warning if level == "critical" else QMessageBox.information)(
        parent, title, button.toolTip()))
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


def show_help_with_logo(parent, dialog_title, help_text):
    box = QMessageBox(parent)
    box.setWindowTitle(dialog_title)
    box.setText(help_text)
    box.setStandardButtons(QMessageBox.Ok)
    if LOGO_PNG.exists():
        pixmap = QPixmap(str(LOGO_PNG)).scaled(
            ABOUT_DIALOG_LOGO_SIZE, ABOUT_DIALOG_LOGO_SIZE, Qt.KeepAspectRatio, Qt.SmoothTransformation,
        )
        box.setIconPixmap(pixmap)
    else:
        box.setIcon(QMessageBox.Information)
    box.exec()
