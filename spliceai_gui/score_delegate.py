"""Draws the Max Score column as a number with a small bar behind it, and a
severity stripe down the first column of each row.

The point is that a table of 1,950 rows can be scanned rather than read. The
colour bands are SpliceAI's own published cutoffs, not invented ones:

    < 0.20   nothing to look at        (grey)
    >= 0.20  high recall               (amber)
    >= 0.50  recommended threshold     (orange-red)
    >= 0.80  high precision            (red)

These are the same four values the threshold buttons above the table use, so a
row's colour and the filter agree by construction (THRESHOLDS below is what
both read).

Custom painting, not a stylesheet: Qt stylesheets can colour a cell but can't
draw inside one. A QStyledItemDelegate is the supported way, and it changes
nothing about the model -- sorting, filtering, copying and export all still see
the plain number, so this is presentation only.
"""
from PySide6.QtCore import QRectF, Qt
from PySide6.QtGui import QColor, QPainter
from PySide6.QtWidgets import QStyle, QStyledItemDelegate

# (lower bound, colour) darkest last; see module docstring for what they mean.
THRESHOLDS = (
    (0.80, "#a3181f"),
    (0.50, "#c2521c"),
    (0.20, "#a8621a"),
)
NEUTRAL = "#8b9aa1"

BAR_HEIGHT = 4
BAR_WIDTH = 46
BAR_GAP = 8
STRIPE_WIDTH = 3
STRIPE_INSET = 4


def severity_color(score):
    """The colour for one max score, or None when there is no score at all."""
    if score is None:
        return None
    for lower, color in THRESHOLDS:
        if score >= lower:
            return QColor(color)
    return QColor(NEUTRAL)


def _score_of(index):
    """The row's max score as a float, or None -- read through the model's own
    sort value, so it works the same behind the filter proxy."""
    from .table_model import SORT_ROLE

    value = index.data(SORT_ROLE)
    if value is None or value == float("-inf"):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


class ScoreBarDelegate(QStyledItemDelegate):
    """Max Score column: the number, right-aligned as before, with a bar to its
    left filled in proportion to the score."""

    def paint(self, painter, option, index):
        score = _score_of(index)
        if score is None:
            super().paint(painter, option, index)
            return

        # Let Qt draw the background (selection, alternating rows, hover) so
        # the cell keeps the palette and theme the rest of the table has.
        style_option = option.__class__(option)
        self.initStyleOption(style_option, index)
        style_option.text = ""
        widget = option.widget
        style = widget.style() if widget else QStyle.__class__
        if widget:
            style.drawControl(QStyle.CE_ItemViewItem, style_option, painter, widget)

        rect = option.rect.adjusted(4, 0, -6, 0)
        color = severity_color(score)

        painter.save()
        painter.setRenderHint(QPainter.Antialiasing, True)

        bar_rect = QRectF(rect.left(), rect.center().y() - BAR_HEIGHT / 2 + 1, BAR_WIDTH, BAR_HEIGHT)
        track = QColor(color)
        track.setAlpha(48)
        painter.setPen(Qt.NoPen)
        painter.setBrush(track)
        painter.drawRoundedRect(bar_rect, BAR_HEIGHT / 2, BAR_HEIGHT / 2)
        filled = QRectF(bar_rect)
        filled.setWidth(max(BAR_HEIGHT, bar_rect.width() * min(max(score, 0.0), 1.0)))
        painter.setBrush(color)
        painter.drawRoundedRect(filled, BAR_HEIGHT / 2, BAR_HEIGHT / 2)

        text_rect = rect.adjusted(int(BAR_WIDTH + BAR_GAP), 0, 0, 0)
        # A selected row paints its own high-contrast text colour; outside a
        # selection the severity colour carries the meaning.
        selected = bool(style_option.state & QStyle.State_Selected)
        painter.setPen(style_option.palette.highlightedText().color() if selected else color)
        font = painter.font()
        font.setBold(score >= THRESHOLDS[-1][0])
        painter.setFont(font)
        painter.drawText(text_rect, Qt.AlignRight | Qt.AlignVCenter, f"{score:.2f}")
        painter.restore()

    def sizeHint(self, option, index):
        size = super().sizeHint(option, index)
        size.setWidth(max(size.width(), BAR_WIDTH + BAR_GAP + 44))
        return size


class SeverityStripeDelegate(QStyledItemDelegate):
    """First column: a short coloured bar at the very left of the row, so the
    eye can find the rows worth reading without reading any of them."""

    def __init__(self, max_score_column, parent=None):
        super().__init__(parent)
        self._max_score_column = max_score_column

    def paint(self, painter, option, index):
        super().paint(painter, option, index)
        score_index = index.siblingAtColumn(self._max_score_column)
        color = severity_color(_score_of(score_index)) if score_index.isValid() else None
        if color is None:
            return
        rect = option.rect
        painter.save()
        painter.setRenderHint(QPainter.Antialiasing, True)
        painter.setPen(Qt.NoPen)
        painter.setBrush(color)
        painter.drawRoundedRect(
            QRectF(rect.left() + 1, rect.top() + STRIPE_INSET, STRIPE_WIDTH,
                   max(2, rect.height() - 2 * STRIPE_INSET)),
            STRIPE_WIDTH / 2, STRIPE_WIDTH / 2,
        )
        painter.restore()
