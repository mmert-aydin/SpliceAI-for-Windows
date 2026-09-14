"""The application's look: one palette, one stylesheet, applied once at startup.

Before this, colours were written inline at each widget that needed one --
about a dozen different blues, greens and reds, set with setStyleSheet() calls
scattered through main_window. They are collected here instead, so there is one
place to change how the program looks and one set of colours to keep coherent.

The rules it follows:

- One accent (TEAL) for everything interactive. Green, amber and red are kept
  for meaning only: ready / warning / problem, and the score severity bands in
  score_delegate.py. When something is red here, it is red because it matters.
- Neutrals are cooled very slightly toward the accent rather than left pure
  grey, so the greys look chosen next to the teal instead of muddy.
- Fonts come from Windows: Segoe UI for text, Consolas for numbers and
  identifiers, which line up in columns. Nothing is bundled and nothing is
  downloaded -- see MONO_FAMILIES.

Qt stylesheets are not CSS. Checkbox and radio indicators are deliberately left
unstyled: giving them a border without also supplying an image makes them
disappear on some Windows themes, and a tick you can't see is worse than a
plain one.
"""
from PySide6.QtGui import QFont

# --- palette ----------------------------------------------------------------
GROUND = "#eef2f3"        # window behind everything
SURFACE = "#ffffff"       # input fields, table, raised things
SURFACE_2 = "#f6f8f9"     # buttons at rest, table header, alternating rows
SURFACE_3 = "#e6edef"     # pressed / track
LINE = "#d3dde0"          # borders
LINE_SOFT = "#e4ebed"     # grid lines, hairlines inside a panel

INK = "#0f1a1e"           # primary text
INK_2 = "#48595f"         # secondary text
INK_3 = "#74878e"         # labels, hints, disabled

TEAL = "#0e7c86"          # the one accent
TEAL_DARK = "#0a626a"     # pressed
TEAL_SOFT = "#d8edef"     # selection, soft fills
TEAL_LINE = "#9fcdd1"

OK = "#1a7f37"
OK_SOFT = "#dcf2e3"
WARN = "#a8621a"
WARN_SOFT = "#fbeada"
CRIT = "#a3181f"
CRIT_SOFT = "#fadfe0"

# Consolas and Cascadia Mono ship with Windows; Courier New is the last resort.
MONO_FAMILIES = ("Consolas", "Cascadia Mono", "Courier New")

_mono_font = None


def mono_font(point_size=None):
    """A monospaced font for numbers, positions, transcripts and paths, so
    digits line up down a column. Built once and reused -- the table asks for
    it per cell."""
    global _mono_font
    if _mono_font is None:
        font = QFont()
        font.setFamilies(list(MONO_FAMILIES))
        font.setStyleHint(QFont.Monospace)
        _mono_font = font
    if point_size is None:
        return _mono_font
    sized = QFont(_mono_font)
    sized.setPointSize(point_size)
    return sized


# --- ready-made label styles, for the few places that set one by hand --------
# Kept here so "what colour is a warning" has a single answer.
STATUS_OK = f"color: {OK}; font-weight: 600;"
STATUS_WARN = f"color: {WARN};"
STATUS_CRIT = f"color: {CRIT}; font-weight: 600;"
STATUS_MUTED = f"color: {INK_3};"


def pill(kind="muted"):
    """A small rounded status chip, as used in the top bar."""
    background, color, border = {
        "ok": (OK_SOFT, OK, OK),
        "warn": (WARN_SOFT, WARN, WARN),
        "crit": (CRIT_SOFT, CRIT, CRIT),
        "muted": (SURFACE_2, INK_2, LINE),
    }[kind]
    return (
        f"QPushButton, QLabel {{ background: {background}; color: {color};"
        f" border: 1px solid {border}; border-radius: 9px; padding: 2px 10px;"
        f" font-weight: 600; text-align: left; }}"
    )


STYLESHEET = f"""
QWidget {{
    background: {GROUND};
    color: {INK};
}}

QToolTip {{
    background: {SURFACE};
    color: {INK};
    border: 1px solid {LINE};
    padding: 5px 7px;
}}

/* Sections: a small label and a hairline instead of a boxed frame, so there is
   less ink around each group and more attention left for what is inside it. */
QGroupBox {{
    border: none;
    border-top: 1px solid {LINE};
    margin-top: 16px;
    padding-top: 10px;
    font-weight: 600;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 0px;
    padding: 0 0 4px 0;
    color: {INK_3};
    font-size: 8pt;
    font-weight: 600;
}}

QLabel {{ background: transparent; }}

QPushButton {{
    background: {SURFACE_2};
    border: 1px solid {LINE};
    border-radius: 5px;
    padding: 5px 13px;
    color: {INK};
}}
QPushButton:hover {{ background: {SURFACE}; border-color: {TEAL_LINE}; }}
QPushButton:pressed {{ background: {SURFACE_3}; }}
QPushButton:disabled {{ color: {INK_3}; border-color: {LINE_SOFT}; background: {SURFACE_2}; }}
QPushButton:focus {{ border-color: {TEAL}; }}

/* The one button that starts the work. */
QPushButton#runButton {{
    background: {TEAL};
    border-color: {TEAL};
    color: #ffffff;
    font-weight: 600;
    padding: 5px 20px;
}}
QPushButton#runButton:hover {{ background: {TEAL_DARK}; border-color: {TEAL_DARK}; }}
QPushButton#runButton:disabled {{
    background: {SURFACE_2}; border-color: {LINE}; color: {INK_3};
}}

QLineEdit, QPlainTextEdit, QTextEdit, QComboBox, QSpinBox {{
    background: {SURFACE};
    border: 1px solid {LINE};
    border-radius: 5px;
    padding: 4px 8px;
    selection-background-color: {TEAL_SOFT};
    selection-color: {INK};
}}
QLineEdit:focus, QPlainTextEdit:focus, QTextEdit:focus, QComboBox:focus {{
    border-color: {TEAL};
}}
QLineEdit:disabled, QComboBox:disabled {{ background: {SURFACE_2}; color: {INK_3}; }}
QComboBox::drop-down {{ border: none; width: 18px; }}

QTableView {{
    background: {SURFACE};
    alternate-background-color: {SURFACE_2};
    border: 1px solid {LINE};
    border-radius: 5px;
    gridline-color: {LINE_SOFT};
    selection-background-color: {TEAL_SOFT};
    selection-color: {INK};
}}
QTableView::item {{ padding: 3px 6px; }}
QHeaderView::section {{
    background: {SURFACE_2};
    color: {INK_3};
    border: none;
    border-bottom: 1px solid {LINE};
    border-right: 1px solid {LINE_SOFT};
    padding: 5px 6px;
    font-weight: 600;
}}
QHeaderView::section:hover {{ color: {TEAL}; }}
QTableCornerButton::section {{ background: {SURFACE_2}; border: none; }}

QProgressBar {{
    background: {SURFACE_3};
    border: none;
    border-radius: 3px;
    height: 6px;
    text-align: center;
    color: transparent;
}}
QProgressBar::chunk {{ background: {TEAL}; border-radius: 3px; }}

/* The divider between panes: two short lines, so it reads as something to
   drag rather than a gap that happens to be there. */
QSplitter::handle {{
    background: {GROUND};
    border-top: 1px solid {LINE_SOFT};
    border-bottom: 1px solid {LINE_SOFT};
}}
QSplitter::handle:vertical {{ height: 7px; }}
QSplitter::handle:hover {{ background: {TEAL_SOFT}; border-color: {TEAL_LINE}; }}
QSplitter::handle:pressed {{ background: {TEAL_LINE}; }}

QScrollArea {{ background: {GROUND}; border: none; }}

/* Scrollbars. Windows' own are a different grey from everything else here and
   make the window look older than it is. These are a plain track with a solid
   handle in the ink colour, going full accent while dragged -- no arrow
   buttons, which nobody clicks and which are what date the native ones. */
QScrollBar:vertical {{
    background: {SURFACE_2};
    width: 12px;
    margin: 0;
    border: none;
    border-left: 1px solid {LINE_SOFT};
}}
QScrollBar:horizontal {{
    background: {SURFACE_2};
    height: 12px;
    margin: 0;
    border: none;
    border-top: 1px solid {LINE_SOFT};
}}
QScrollBar::handle:vertical {{
    background: {INK_3};
    border-radius: 4px;
    min-height: 28px;
    margin: 2px;
}}
QScrollBar::handle:horizontal {{
    background: {INK_3};
    border-radius: 4px;
    min-width: 28px;
    margin: 2px;
}}
QScrollBar::handle:hover {{ background: {INK_2}; }}
QScrollBar::handle:pressed {{ background: {TEAL}; }}
/* The arrow buttons and the empty track above/below the handle: Qt draws its
   own unless each is given a zero size and no background. */
QScrollBar::add-line, QScrollBar::sub-line {{
    width: 0; height: 0; background: none; border: none;
}}
QScrollBar::add-page, QScrollBar::sub-page {{ background: none; }}

QMenu {{ background: {SURFACE}; border: 1px solid {LINE}; padding: 4px; }}
QMenu::item {{ padding: 5px 24px 5px 12px; border-radius: 4px; }}
QMenu::item:selected {{ background: {TEAL_SOFT}; color: {INK}; }}

QDialog {{ background: {GROUND}; }}
"""


def apply_theme(app):
    """Applies the stylesheet to the whole application. Called once from
    main(), before any window is built."""
    app.setStyleSheet(STYLESHEET)
