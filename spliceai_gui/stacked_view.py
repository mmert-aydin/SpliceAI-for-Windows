"""The viewer's stacked base-level view (spec 8.2 A): one custom widget, no
item per base, drawing only the columns in view.

Rows, top to bottom: legend, the events' labels, ruler, reference DNA, reference scores, the
canonical transcript (and, on request, the curated alternatives), the
wild type's predicted transcript (and on request its alternatives), the
sample's sequence, the sample's scores and the delta, and the sample's
predicted transcripts (candidate 1, the others on request). Score rows are
heat strips (intensity = probability, acceptor and donor in two fixed hues)
when zoomed out -- from per-pixel maxima of a pyramid built once --, bars
with the site threshold as a dashed line in between, and letters and
numbers at base level.

Spec 15: the screen is drawn from tiles 256 pixels wide, rendered in a
background thread and cached, so panning only moves pictures already made
and the view opens at once with what is ready; while zooming, the tiles of
the last zoom stand in, scaled, until the new ones arrive. Labels, legend
and cursor are drawn on top. TrackPainter draws the same for a tile and for
a PNG / SVG export.

Spec 17 item 5: tiles are rendered at the screen's device pixel ratio and placed
on whole device pixels. The ruler's numbers and the exon labels aren't in the
tiles: they are drawn on top at their final size in every frame, so they stay
sharp at any scaling and while scaled stand-ins are shown. The inserted columns
of an insertion are a labelled gap in the reference rows ("193 bp inserted"),
the transcripts running on across it as a dashed line. The gene's strand is in the viewer's header and as
chevrons on intron lines; in transcript orientation a minus-strand gene is
drawn reverse-complemented (ViewData.flipped), exon 1 on the left. Hovering
shows a short tooltip; a click (or Shift + hover) sends the full details to
the viewer's side panel.

Spec 16: every applied event is a translucent band through all rows in its
type's colour (an insertion: a line in the reference rows, the inserted
block in the sample's), labelled above the ruler with its type, size,
coordinates and place in the canonical transcript; "Show origin" outlines
where inserted bases come from. Focus shows only the reference, WT predicted
1 and sample 1 rows. The delta row is drawn on top of the tiles, scaled to
the largest |delta| in view, with positions of |delta| >= 0.1 marked.
"""
import math
import time
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QFont, QFontMetrics, QImage, QPainter, QPen
from PySide6.QtWidgets import QMenu, QToolTip, QWidget

from spliceai_pipeline import cdna_analysis, haplotype

from . import theme
from .i18n import N_, t

BASE_COLORS = {"A": "#2e9d46", "C": "#2167c4", "G": "#d38a00", "T": "#d1262b", "N": "#8a8a8a", "-": "#ffffff"}
# Exon identity colours (by the reference exon an exon matches); novel exons are grey and hatched.
EXON_COLORS = ["#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377", "#BBBBBB"]
NOVEL = "#555555"
FOOTER = 34             # px under an exported image, for the research-use note (two lines)
CHANGED = "#f4b6b6"     # substituted bases, and the gap of deleted ones, in the sample's row
LEGEND, RULER, DNA, SCORES, TRANSCRIPT, DELTA, EVENTS = range(7)
# Spec 16 item 1: one colour per event type, in the bands, the labels, the legend and the schematic.
EVENT_COLORS = {"deletion": "#D62728", "duplication": "#2CA02C", "inversion": "#8E44AD", "insertion": "#C79100",
                "complex": "#0E9AA7", "SNV": "#D0459A"}
EVENT_NAMES = {"deletion": "deletion", "duplication": "duplication", "inversion": "inversion",
               "insertion": "insertion", "complex": "complex / delins", "SNV": "SNV"}
DELTA_MARK = 0.1        # |delta| marked on the delta row
LANE = 14               # px per line of event labels
GAP_MIN_PX = 6          # an inserted gap narrower than this isn't labelled or dashed
MAX_PX_PER_COLUMN = 40
TILE = 256              # px per tile
MAX_TILES = 160         # cached; about 100 MB at the default height
CHEVRON_EVERY = 48      # px between chevrons on an intron line
TICK_EVERY = 110        # px between ruler ticks, about
ROW_GAP = 9             # px between rows
# A background per group of rows, so the reference, the WT prediction and the sample read apart; every sample
# candidate shares one.
GROUP_BACKGROUNDS = {"ref": "#eef3f9", "wt": "#eef7ef", "pat": "#fdf4e8"}


def exon_color(identity):
    return QColor(EXON_COLORS[(identity - 1) % (len(EXON_COLORS) - 1)]) if identity else QColor(NOVEL)


def rightwards(data):
    """Transcription runs left to right on screen (plus strand, or a minus-strand gene flipped)."""
    return (data.strand == "+") != bool(getattr(data, "reversed", False))


def orientation_text(data, tr=None):
    """The orientation, in English (exported images) or through tr (i18n.t, the window's label)."""
    tr = tr or (lambda text, **fields: text.format(**fields))
    strand = tr(N_("+ strand")) if data.strand == "+" else tr(N_("− strand"))
    if getattr(data, "reversed", False):
        return tr(N_("transcript orientation 5'→3' ({strand}, reverse-complemented)"), strand=strand)
    if rightwards(data):
        return tr(N_("genomic orientation ({strand}, 5'→3' points right)"), strand=strand)
    return tr(N_("genomic orientation ({strand}, 5'→3' points left)"), strand=strand)


def event_color(kind, alpha=255):
    color = QColor(EVENT_COLORS.get(kind, "#555555"))
    color.setAlpha(alpha)
    return color


def exon_label(fm, identity, px, duplicate=False):
    """Spec 16: what an exon block px wide has room for -- "Exon 14", else "14", else nothing."""
    names = (f"Exon {identity}", str(identity)) if identity else ("novel", "n")
    if duplicate:
        names = (f"duplicate of exon {identity}", f"dup {identity}", str(identity))
    return next((t for t in names if fm.horizontalAdvance(t) + 4 <= px), "")


def _nice(x):
    """1, 2 or 5 times a power of ten, at least x."""
    step = 10 ** max(0, int(math.floor(math.log10(max(x, 1)))))
    for mult in (1, 2, 5, 10):
        if step * mult >= x:
            return step * mult
    return step * 10


class TrackPainter:
    """Draws the rows for one stretch of columns: a tile (labels=False) or a whole export."""

    def __init__(self, data, rows, first, cols_per_px, labels=True, delta_scale=None):
        """delta_scale: what the delta row's full height stands for; None (tiles) leaves the row out -- the
        view draws it on top, scaled to what is in view."""
        self.data, self._rows, self.first, self.cols_per_px, self.labels = data, rows, first, cols_per_px, labels
        self.delta_scale = delta_scale
        self.background = "#ffffff"         # what a heat strip fades into: the row group's background

    def x_of(self, c):
        return (c - self.first) / self.cols_per_px

    def columns(self, width):
        return max(0, int(self.first)), min(self.data.n - 1, int(math.ceil(self.first + width * self.cols_per_px)))

    def paint(self, p, width):
        y = 2
        c0, c1 = self.columns(width)
        p.setFont(view_font())
        for kind, h, payload in self._rows:
            rect = QRectF(0, y, width, h)
            group = row_group(kind, payload)
            self.background = GROUP_BACKGROUNDS[group] if group else "#ffffff"
            if group:
                p.fillRect(QRectF(0, y - ROW_GAP / 2 + 1, width, h + ROW_GAP - 2), QColor(GROUP_BACKGROUNDS[group]))
            if kind == LEGEND:
                if self.labels:
                    legend(p, rect, self.data)
            elif kind == RULER:
                self._ruler(p, rect, c0, c1)
            elif kind == DNA:
                self._dna(p, rect, c0, c1, payload)
            elif kind == SCORES:
                self._scores(p, QRectF(0, y + 12, width, (h - 12) / 2), "ref_acc" if payload == "ref" else "pat_acc",
                             theme.ACCEPTOR, c0, c1)
                self._scores(p, QRectF(0, y + 12 + (h - 12) / 2, width, (h - 12) / 2),
                             "ref_don" if payload == "ref" else "pat_don", theme.DONOR, c0, c1)
                self._hatch_unscored(p, QRectF(0, y + 12, width, h - 12), c0, c1)
            elif kind == DELTA:
                if self.delta_scale is not None:
                    self.delta(p, rect, c0, c1, self.delta_scale)
            elif kind == TRANSCRIPT:
                self._transcript(p, rect, payload, c0, c1)
            elif kind == EVENTS and self.labels:
                event_labels(p, rect, self.data, self.x_of)
            if self.labels:
                row_label(p, rect, kind, payload, self.data)
            y += h + ROW_GAP
        if self.labels:
            event_bands(p, self.data, self._rows, self.x_of, width)
            insertion_gaps(p, self.data, self._rows, self.x_of, width)
            self.text(p, width)

    def text(self, p, width):
        """What is written on the rows and never put in a tile: the ruler's numbers and the exon labels."""
        c0, c1 = self.columns(width)
        y = 2
        for kind, h, payload in self._rows:
            if kind == RULER:
                self._ruler_text(p, QRectF(0, y, width, h), c0, c1)
            elif kind == TRANSCRIPT:
                self._exon_text(p, QRectF(0, y, width, h), payload, c0, c1)
            y += h + ROW_GAP

    def delta(self, p, rect, c0, c1, scale):
        """The delta row: acceptor above, donor below, the full half-height standing for `scale`; positions
        with |delta| >= DELTA_MARK get a tick along the top (their x pixels are returned)."""
        half = (rect.height() - 12) / 2
        marks = self._scores(p, QRectF(0, rect.y() + 12, rect.width(), half), "delta_acc", theme.ACCEPTOR, c0, c1,
                             signed=True, scale=scale)
        marks = np.union1d(marks, self._scores(p, QRectF(0, rect.y() + 12 + half, rect.width(), half), "delta_don",
                                               theme.DONOR, c0, c1, signed=True, scale=scale))
        if len(marks):
            breaks = np.nonzero(np.diff(marks) > 1)[0]
            for a, b in zip(np.r_[marks[0], marks[breaks + 1]], np.r_[marks[breaks], marks[-1]]):
                p.fillRect(QRectF(a - 1, rect.y() + 7, b - a + 3, 5), QColor(theme.CRIT))
        return marks

    def _ticks(self, c0, c1):
        """[(position, x)] of the ruler: round steps that depend on the zoom only, so tiles join up."""
        d = self.data
        lo_c, hi_c = max(0, c0 - int(80 * self.cols_per_px) - 1), min(d.n - 1, c1 + int(80 * self.cols_per_px) + 1)
        refs = d.ref[lo_c:hi_c + 1]
        placed = refs[refs > 0]
        if not len(placed):
            return []
        step = _nice(TICK_EVERY * self.cols_per_px)
        lo, hi = int(placed.min()), int(placed.max())
        return [(tick, self.x_of(d.column_of_ref(tick) + 0.5)) for tick in range(-(-lo // step) * step, hi + 1, step)]

    def _ruler(self, p, rect, c0, c1):
        p.setPen(QColor(theme.INK_2))
        p.drawLine(QPointF(0, rect.bottom() - 6), QPointF(rect.width(), rect.bottom() - 6))
        for _tick, x in self._ticks(c0, c1):
            p.drawLine(QPointF(x, rect.bottom() - 10), QPointF(x, rect.bottom() - 6))

    def _ruler_text(self, p, rect, c0, c1):
        """Genomic positions; inside inserted sequence its haplotype labels (INS1+n)."""
        d = self.data
        p.setPen(QColor(theme.INK_2))
        for tick, x in self._ticks(c0, c1):
            p.drawText(QRectF(round(x) - 60, rect.y(), 120, 12), Qt.AlignCenter, f"{tick:,}")
        if 1 / self.cols_per_px >= 8:
            for c in range(c0, c1 + 1):
                if d.kind[c] == haplotype.INSERTED and (c == 0 or d.kind[c - 1] != haplotype.INSERTED):
                    x = round(self.x_of(c))
                    p.setPen(QColor(theme.WARN))
                    p.drawText(QRectF(x, rect.y() + 11, 140, 10), Qt.AlignLeft, d.label(c))
                    p.setPen(QColor(theme.INK_2))

    def _dna(self, p, rect, c0, c1, which):
        d = self.data
        bases = d.ref_base if which == "ref" else d.pat_base
        px = 1 / self.cols_per_px
        if px < 1:
            p.fillRect(rect, QColor("#e9edef"))
            if which == "pat":
                self._spans(p, rect, c0, c1, d.kind == haplotype.INSERTED, "#fff2a8")
                self._spans(p, rect, c0, c1, d.kind == haplotype.INVERTED, "#e8d5ff")
                self._spans(p, rect, c0, c1, (d.kind == haplotype.SUBSTITUTED) | (d.kind == haplotype.DELETED),
                            "#d1262b")
            return
        for c in range(c0, c1 + 1):
            x = self.x_of(c)
            b = bases[c]
            color = QColor(BASE_COLORS.get(b, "#8a8a8a"))
            if which == "pat":
                if d.kind[c] == haplotype.INSERTED:
                    p.fillRect(QRectF(x, rect.y(), px, rect.height()), QColor("#fff2a8"))
                elif d.kind[c] == haplotype.INVERTED:
                    p.fillRect(QRectF(x, rect.y(), px, rect.height()), QColor("#e8d5ff"))
                elif d.kind[c] in (haplotype.SUBSTITUTED, haplotype.DELETED):
                    p.fillRect(QRectF(x, rect.y(), px, rect.height()), QColor(CHANGED))
                if not d.provided[c]:
                    color = QColor("#cccccc")
            if px >= 8:
                p.setPen(color)
                p.drawText(QRectF(x, rect.y(), px, rect.height()), Qt.AlignCenter,
                           ("◄" if d.kind[c] == haplotype.INVERTED and px < 14 else b) if b != "-" else "-")
            elif b != "-":
                p.fillRect(QRectF(x, rect.y() + 3, max(px, 1), rect.height() - 6), color)

    def _spans(self, p, rect, c0, c1, mask, color):
        idx = np.nonzero(mask[c0:c1 + 1])[0]
        if len(idx):
            x = (idx + c0 - self.first) / self.cols_per_px
            for xi in np.unique(x.astype(int)):
                p.fillRect(QRectF(xi, rect.y(), 1, rect.height()), QColor(color))

    def _scores(self, p, rect, track, hue, c0, c1, signed=False, scale=1.0):
        """One track: a heat strip when zoomed out, bars when zoomed in; values are drawn as value / scale.
        Returns the x pixels (ints) where |value| >= DELTA_MARK, for signed tracks."""
        d = self.data
        width = int(math.ceil(rect.width()))
        color = QColor(hue)
        px = 1 / self.cols_per_px
        if px < 2:
            levels = d.pyramid(track)
            level = max(0, min(len(levels) - 1, int(math.floor(math.log2(max(self.cols_per_px, 1))))))
            arr = levels[level]
            starts = (self.first + np.arange(width) * self.cols_per_px).astype(np.int64)
            ends = np.minimum((self.first + (np.arange(width) + 1) * self.cols_per_px).astype(np.int64), d.n)
            lo = np.clip(starts >> level, 0, len(arr) - 1)
            hi = np.clip(np.maximum((ends - 1) >> level, lo), 0, len(arr) - 1)
            # Largest value (by size) per pixel: three pyramid entries cover a pixel at most.
            a, m, b = arr[lo], arr[(lo + hi) // 2], arr[hi]
            v = np.where(np.abs(a) >= np.abs(b), a, b)
            v = np.where(np.abs(v) >= np.abs(m), v, m)
            v[starts >= d.n] = 0
            intensity = np.clip(np.abs(v) / scale, 0, 1)
            rgba = np.empty((1, width, 4), np.uint8)
            base = np.array([color.red(), color.green(), color.blue()], np.float32)
            bg = QColor(self.background)
            white = np.array([bg.red(), bg.green(), bg.blue()], np.float32)
            mix = white[None, :] * (1 - intensity[:, None]) + base[None, :] * intensity[:, None]
            rgba[0, :, 0:3] = mix[:, ::-1].astype(np.uint8)       # QImage.Format_RGB32 is BGRA in memory
            rgba[0, :, 3] = 255
            image = QImage(rgba.data, width, 1, width * 4, QImage.Format_RGB32)
            p.drawImage(QRectF(rect.x(), rect.y(), width, rect.height()), image)
            return np.nonzero(np.abs(v) >= DELTA_MARK)[0] if signed else np.array([], int)
        values = d.tracks[track]
        mid = rect.center().y() if signed else rect.bottom()
        unit = rect.height() / (2 if signed else 1)
        marks = []
        p.setPen(Qt.NoPen)
        p.setBrush(QBrush(color))
        for c in range(c0, c1 + 1):
            v = values[c]
            if np.isnan(v) or v == 0:
                continue
            x = self.x_of(c)
            h = float(np.clip(v / scale, -1, 1)) * unit
            p.drawRect(QRectF(x + px * 0.1, mid - max(h, 0), max(px * 0.8, 1), abs(h)))
            if signed and abs(v) >= DELTA_MARK:
                marks.append(int(x + px / 2))
            if px >= 34 and abs(v) >= 0.001:
                p.setPen(QColor(theme.INK))
                p.drawText(QRectF(x, rect.y(), px, 10), Qt.AlignCenter, f"{v:.3f}")
                p.setPen(Qt.NoPen)
        if not signed:
            p.setPen(QPen(QColor(theme.INK_3), 1, Qt.DashLine))
            ty = rect.bottom() - d.threshold * rect.height()
            p.drawLine(QPointF(rect.x(), ty), QPointF(rect.right(), ty))
        p.setBrush(Qt.NoBrush)
        return np.array(marks, int)

    def _hatch_unscored(self, p, rect, c0, c1):
        d = self.data
        with d._lock:
            runs = d.__dict__.get("_unscored")
            if runs is None:
                runs = d.__dict__["_unscored"] = d.unscored()
        for a, b, _why in runs:
            if b < c0 or a > c1:
                continue
            x0, x1 = self.x_of(max(a, c0)), self.x_of(min(b, c1) + 1)
            p.fillRect(QRectF(x0, rect.y(), max(x1 - x0, 1), rect.height()), QBrush(QColor(theme.INK_3), Qt.BDiagPattern))

    def _transcript(self, p, rect, row, c0, c1):
        body = QRectF(rect.x(), rect.y() + 12, rect.width(), 14)
        member = row.membership_range(c0, c1)
        span = np.nonzero(member != 0)[0]
        if len(span):
            x0, x1 = self.x_of(c0 + span[0]), self.x_of(c0 + span[-1] + 1)
            cy = body.center().y()
            p.setPen(QPen(QColor(theme.INK_2), 1))
            p.drawLine(QPointF(x0, cy), QPointF(x1, cy))
            self._chevrons(p, x0, x1, cy, c0, member)
        cds = getattr(row, "cds", None)
        for blk in row.blocks:
            if blk.last < c0 or blk.first > c1:
                continue
            x0, x1 = self.x_of(max(blk.first, c0)), self.x_of(min(blk.last, c1) + 1)
            color = exon_color(blk.identity)
            if cds is not None:     # UTR thin, coding part full height
                p.fillRect(QRectF(x0, body.y() + 4, max(x1 - x0, 1), body.height() - 8), color)
                a, b = max(blk.first, cds[0], c0), min(blk.last, cds[1], c1)
                if a <= b:
                    p.fillRect(QRectF(self.x_of(a), body.y(), max(self.x_of(b + 1) - self.x_of(a), 1),
                                      body.height()), color)
                continue
            p.fillRect(QRectF(x0, body.y(), max(x1 - x0, 1), body.height()), color)
            if not blk.identity:
                p.fillRect(QRectF(x0, body.y(), max(x1 - x0, 1), body.height()), QBrush(QColor("white"), Qt.FDiagPattern))
        # Beneath: this transcript's scores at its exonic positions.
        strip = QRectF(rect.x(), rect.y() + 27, rect.width(), 6)
        if 1 / self.cols_per_px < 2:
            return
        acc = self.data.tracks["pat_acc" if row.on_haplotype else "ref_acc"]
        don = self.data.tracks["pat_don" if row.on_haplotype else "ref_don"]
        for c in range(c0, c1 + 1):
            if member[c - c0] > 0:
                v = max(np.nan_to_num(acc[c]), np.nan_to_num(don[c]))
                if v > 0.01:
                    color = QColor(theme.ACCEPTOR if np.nan_to_num(acc[c]) >= np.nan_to_num(don[c]) else theme.DONOR)
                    color.setAlphaF(min(1.0, float(v)))
                    p.fillRect(QRectF(self.x_of(c), strip.y(), 1 / self.cols_per_px, strip.height()), color)

    def _exon_text(self, p, rect, row, c0, c1):
        """Spec 16: a number on every exon there is room for -- "Exon 14", else "14". Named after the reference
        exon it matches, so the same exon reads the same everywhere."""
        fm = p.fontMetrics()
        for blk in row.blocks:
            if blk.last < c0 or blk.first > c1:
                continue
            x0, x1 = self.x_of(blk.first), self.x_of(blk.last + 1)
            text = exon_label(fm, blk.identity, x1 - x0, blk.duplicate)
            if text:
                # The part in view, so a long exon's label stays on screen when it is zoomed into.
                a, b = max(x0, 0), min(x1, rect.width())
                centre = (a + b) / 2 if b - a > fm.horizontalAdvance(text) + 4 else (x0 + x1) / 2
                p.setPen(QColor(theme.INK) if exon_color(blk.identity).lightness() > 130 else QColor("white"))
                p.drawText(QRectF(round(centre) - 80, rect.y() + 12, 160, 14), Qt.AlignCenter, text)

    def _chevrons(self, p, x0, x1, cy, c0, member):
        """UCSC-style arrows on the intron line, pointing the way the gene is transcribed, at fixed screen
        steps (so tiles join up), only where the line is intron."""
        right = rightwards(self.data)
        origin = self.first / self.cols_per_px                    # this tile's left edge, in whole-view pixels
        k = math.ceil((origin + x0) / CHEVRON_EVERY)
        p.setPen(QPen(QColor(theme.INK_2), 1))
        while k * CHEVRON_EVERY - origin <= x1:
            x = k * CHEVRON_EVERY - origin
            c = int(self.first + x * self.cols_per_px)
            if 0 <= c - c0 < len(member) and member[c - c0] < 0:
                d = 3 if right else -3
                p.drawLine(QPointF(x - d, cy - 3), QPointF(x, cy))
                p.drawLine(QPointF(x - d, cy + 3), QPointF(x, cy))
            k += 1


def row_group(kind, payload):
    """"ref", "wt" or "pat": which group's background a row is drawn on (None: legend, labels, ruler)."""
    if kind in (DNA, SCORES):
        return "ref" if payload == "ref" else "pat"
    if kind == DELTA:
        return "pat"
    if kind == TRANSCRIPT:
        return {"reference": "ref", "alternative": "ref", "wt": "wt"}.get(payload.group, "pat")
    return None


def view_font():
    """The views' small text."""
    font = theme.mono_font(8)
    font.setHintingPreference(QFont.PreferFullHinting)
    return font


def _legend_layout(fm, data, width):
    """([(x, line, colour or None, text)], lines): the legend's entries flowed into lines `width` px wide; the
    event types (spec 16) on a line of their own."""
    entries = [(QColor(c), t) for c, t in ((theme.ACCEPTOR, "acceptor"), (theme.DONOR, "donor"), (NOVEL, "novel exon"),
                                           ("#fff2a8", "inserted"), ("#e8d5ff", "inverted"),
                                           (CHANGED, "changed / deleted"), ("#cccccc", "from the reference"))]
    entries += [(None, f"site threshold {data.threshold:g} (dashed)"), (None, "hatched: not scored")]
    kinds = list(dict.fromkeys(info["kind"] for info in getattr(data, "event_info", [])))
    if kinds:
        entries += [("new line", "Events:")] + [(event_color(k, 110), EVENT_NAMES.get(k, k)) for k in kinds]
    out, x, line = [], 4, 0
    for color, text in entries:
        if isinstance(color, str):
            x, line, color = 4, line + 1, None
        w = (14 if color is not None else 0) + fm.horizontalAdvance(text) + 12
        if x > 4 and x + w > width:
            x, line = 4, line + 1
        out.append((x, line, color, text))
        x += w
    return out, line + 1


def legend_height(data, width):
    return _legend_layout(QFontMetrics(view_font()), data, width)[1] * 16 + 2


def legend(p, rect, data):
    """The bases' and sites' colours; under them (spec 16) the types of the events shown."""
    p.setPen(QColor(theme.INK_2))
    for x, line, color, text in _legend_layout(p.fontMetrics(), data, rect.width())[0]:
        y = rect.y() + line * 16
        if color is not None:
            p.fillRect(QRectF(rect.x() + x, y + 4, 10, 10), color)
        p.drawText(QRectF(rect.x() + x + (13 if color is not None else 0), y, 400, 18),
                   Qt.AlignLeft | Qt.AlignVCenter, text)


def _event_spans(data, x_of, width):
    """[(x0, x1, info, insertion)] of the events on screen, at least 2 px wide."""
    out = []
    for (first, last, _text), info in zip(data.events, getattr(data, "event_info", [])):
        x0, x1 = x_of(first), x_of(last + 1)
        if x1 - x0 < 2:
            x0, x1 = (x0 + x1) / 2 - 1, (x0 + x1) / 2 + 1
        if x1 >= 0 and x0 <= width:
            out.append((x0, x1, info, info["event"].end < info["event"].start))
    return out


def _origin_spans(data, x_of, width):
    """[(x0, x1, info)] of the regions on screen that events' inserted bases come from."""
    out = []
    for info in getattr(data, "event_info", []):
        for r0, r1, _o in info["sources"]:
            a, b = data.columns_of(r0, r1)
            x0, x1 = x_of(a), x_of(b + 1)
            if x1 >= 0 and x0 <= width:
                out.append((x0, max(x1, x0 + 2), info))
    return out


def event_labels(p, rect, data, x_of, placed=None, origins=False):
    """Spec 16 item 1: above the ruler, each event's label ("DEL 133 bp · chr19:... · intron 13 (c....)"), over
    its band, in two lines when labels would overlap; with origins, "origin of the DUP 121 bp" over each region
    inserted bases come from (item 2). placed: a list that gets (rect, info) of each event's label."""
    fm = p.fontMetrics()
    ends = [-1e9] * max(1, int(rect.height() // LANE))
    items = [(x0, x1, info, info["label"], False) for x0, x1, info, _ins in _event_spans(data, x_of, rect.width())]
    if origins:
        items += [(x0, x1, info, f"origin of the {info['label'].split(' · ')[0]}", True)
                  for x0, x1, info in _origin_spans(data, x_of, rect.width())]
    for x0, x1, info, label, origin in sorted(items, key=lambda s: s[0]):
        text = fm.elidedText(label, Qt.ElideRight, int(rect.width() - 24))
        w = fm.horizontalAdvance(text) + 18
        x = min(max((x0 + x1) / 2 - w / 2, 2), rect.width() - w - 2)
        lane = next((i for i, end in enumerate(ends) if end < x), None)
        if lane is None:            # every line is taken here: the short name only
            lane = min(range(len(ends)), key=lambda i: ends[i])
            text = label.split(" · ")[0]
            w = fm.horizontalAdvance(text) + 18
            x = max(x, ends[lane] + 4)
        ends[lane] = x + w
        box = QRectF(x, rect.y() + lane * LANE, w, LANE - 1)
        p.fillRect(box, QColor(255, 255, 255, 230))
        if origin:
            p.setPen(QPen(event_color(info["kind"]), 1, Qt.DashLine))
            p.setBrush(Qt.NoBrush)
            p.drawRect(QRectF(x + 3, box.y() + 3, 8, 8))
        else:
            p.fillRect(QRectF(x + 3, box.y() + 3, 8, 8), event_color(info["kind"]))
        p.setPen(event_color(info["kind"]).darker(140))
        p.drawText(box.adjusted(14, 0, 0, 0), Qt.AlignLeft | Qt.AlignVCenter, text)
        if placed is not None and not origin:
            placed.append((box, info))


def event_bands(p, data, rows, x_of, width, origins=False):
    """Spec 16 items 1 and 2: a translucent band through every row over each event, in its type's colour; for
    an insertion (a duplication too), a line where it goes in the reference rows and the band over the inserted
    bases in the sample's. origins: also a dashed outline over where its inserted bases come from."""
    y, ruler, split, bottom = 2, None, None, 2
    for kind, h, payload in rows:
        if kind == RULER:
            ruler = y
        if kind == DNA and payload == "pat":
            split = y - 2
        y += h + ROW_GAP
        bottom = y - ROW_GAP
    if ruler is None:
        return
    split = split or bottom
    for x0, x1, info, insertion in _event_spans(data, x_of, width):
        fill, edge = event_color(info["kind"], 46), event_color(info["kind"], 170)
        if insertion:
            p.fillRect(QRectF(x0 - 1, ruler, 2, split - ruler), edge)
            top = split
        else:
            top = ruler
        p.fillRect(QRectF(x0, top, x1 - x0, bottom - top), fill)
        p.setPen(QPen(edge, 1))
        p.drawLine(QPointF(x0, top), QPointF(x0, bottom))
        p.drawLine(QPointF(x1, top), QPointF(x1, bottom))
    if origins:
        for x0, x1, info in _origin_spans(data, x_of, width):
            p.setPen(QPen(event_color(info["kind"]), 2, Qt.DashLine))
            p.setBrush(Qt.NoBrush)
            p.drawRect(QRectF(x0, ruler, x1 - x0, bottom - ruler))


def insertion_gaps(p, data, rows, x_of, width):
    """Spec 17 item 5: an insertion's columns in the reference rows are a gap. It is labelled in the reference
    DNA row ("193 bp inserted"), and a transcript that runs across it goes on as a dashed line."""
    gaps = []
    for (first, last, _text), info in zip(data.events, getattr(data, "event_info", [])):
        e = info["event"]
        x0, x1 = x_of(first), x_of(last + 1)
        if e.end < e.start and x1 >= 0 and x0 <= width and x1 - x0 >= GAP_MIN_PX:
            gaps.append((first, last, x0, x1, info, f"{len(e.inserted):,} bp inserted"))
    y = 2
    for kind, h, payload in rows if gaps else ():
        for first, last, x0, x1, info, text in gaps:
            if kind == DNA and payload == "ref":
                w = p.fontMetrics().horizontalAdvance(text) + 8
                box = QRectF(round((x0 + x1) / 2 - w / 2), y + 1, w, h - 2)
                p.fillRect(box, QColor(255, 255, 255, 230))
                p.setPen(event_color(info["kind"]).darker(140))
                p.drawText(box, Qt.AlignCenter, text)
            elif kind == TRANSCRIPT and not payload.on_haplotype and first > 0 and last + 1 < data.n \
                    and payload.membership_range(first - 1, first - 1)[0] and payload.membership_range(last + 1, last + 1)[0]:
                p.fillRect(QRectF(x0, y + 12, x1 - x0, 14), QColor(theme.SURFACE))
                p.setPen(QPen(QColor(theme.INK_2), 1, Qt.DashLine))
                p.drawLine(QPointF(x0, y + 19), QPointF(x1, y + 19))
        y += h + ROW_GAP


def row_label(p, rect, kind, payload, data, menu=False):
    """The row's name at its left, on a light backing. menu (on screen, not in an exported image): a sample
    candidate's name ends in a "▾" when a click on it lists other candidates."""
    if kind in (RULER, LEGEND, EVENTS):
        return
    text = {DNA: "Reference DNA" if payload == "ref" else "Sample sequence",
            SCORES: "Reference scores" if payload == "ref" else "Sample scores",
            DELTA: "Delta"}.get(kind)
    if kind == TRANSCRIPT:
        text = payload.name + (f"  relative score {payload.relative:.3f}" if payload.relative is not None else "")
        if menu and payload.group == "patient" and sum(r.group == "patient" for r in data.rows) > 1:
            text += "  ▾"       # a click on the name picks another sample candidate
    if not text:
        return
    w = p.fontMetrics().horizontalAdvance(text) + 8
    p.fillRect(QRectF(rect.x() + 2, rect.y(), w, 12), QColor(255, 255, 255, 200))
    p.setPen(QColor(theme.INK_2))
    p.drawText(QRectF(rect.x() + 4, rect.y(), rect.width(), 12), Qt.AlignLeft | Qt.AlignVCenter, text)


def event_pairs(info):
    """An event's (key, value) rows for the side panel."""
    parts = info["label"].split(" · ")
    pairs = [(N_("Type, size"), parts[0]), (N_("Where"), t(info["where"]) + (f" ({info['c']})" if info["c"] else "")),
             (N_("Coordinates"), parts[1]), (N_("Nearest sites"), info["sites"])]
    if info["origin"]:
        pairs.append((N_("Origin"), info["origin"]))
    if info["effect"]:
        pairs.append((N_("Sample 1 there"), info["effect"].split(": ", 1)[1]))
    return pairs


class StackedView(QWidget):
    hovered = Signal(int)       # the genomic column under the pointer
    moved = Signal()            # a frame was drawn at another place or zoom
    details = Signal(object)    # everything about the base, exon or event clicked (or Shift-hovered): sections
    sample_changed = Signal(int)    # the sample candidate shown
    _tile_ready = Signal(object, object)

    def __init__(self, data, parent=None):
        super().__init__(parent)
        self.genomic = data             # columns in genomic order (what the rest of the viewer speaks)
        self.data = data                # what is drawn: genomic, or flipped (transcript orientation)
        self._flipped = None
        self.show_alternatives = False      # curated alternative NM_ (annotation)
        self.show_spliceai_alternatives = False
        self.focus = False                  # spec 16: only the reference, WT predicted 1 and sample 1 rows
        self.sample_rank = 1                # the sample candidate shown (a click on its name picks another)
        self._before_whole = None           # (first, cols_per_px) the whole-gene toggle goes back to
        self.show_origin = False            # spec 16 item 2: outline where inserted bases come from
        self.start_view = None              # genomic columns to open on (else the whole gene)
        self.label_boxes = []               # (rect, event info) of the event labels drawn
        self._legend = (None, None, 18)     # (width, data, height) of the legend, flowed to the width
        self.first, self.cols_per_px = 0.0, 1.0
        self.cursor_col = None
        self._drag = None
        self._moved = False
        self.setMouseTracking(True)
        self.setMinimumHeight(300)
        self._fitted = False
        self._tiles = OrderedDict()     # (layout, zoom, index) -> QImage
        self._pending = set()
        self._recent_zooms = []         # most recent last
        self._pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix="viewer-tiles")
        self._tile_ready.connect(self._on_tile)
        self.frame_ms = []              # paint times, for the speed report

    # --- orientation -----------------------------------------------------------------------
    def set_transcript_orientation(self, on):
        """Spec 15 item 1: a minus-strand gene drawn 5'->3' (reverse-complemented), or genomic. Keeps the
        same stretch in view."""
        want = bool(on) and self.genomic.strand == "-"
        if want == bool(self.data.reversed):
            return
        a, b = self.visible_genomic()
        if want and self._flipped is None:
            self._flipped = self.genomic.flipped()
        self.data = self._flipped if want else self.genomic
        self._clear_tiles()
        self.show_genomic(a, b, margin=0)

    def to_display(self, c):
        return self.data.n - 1 - c if self.data.reversed else c

    def to_genomic(self, c):
        return self.to_display(c)

    def show_genomic(self, first, last, margin=0.1):
        """zoom_to() for genomic columns."""
        a, b = sorted((self.to_display(first), self.to_display(last)))
        self.zoom_to(a, b, margin)

    def visible_genomic(self):
        a, b = self.visible()
        return tuple(sorted((self.to_genomic(a), self.to_genomic(b))))

    # --- geometry ---------------------------------------------------------------------------
    def rows(self):
        """[(kind, height, payload)] of what is drawn, top to bottom."""
        d = self.data
        alternatives = self.show_alternatives and not self.focus
        others = self.show_spliceai_alternatives and not self.focus
        if self._legend[:2] != (self.width(), id(d)):
            self._legend = (self.width(), id(d), legend_height(d, self.width()))
        out = [(LEGEND, self._legend[2], None)]
        if d.events:
            sources = any(info["sources"] for info in getattr(d, "event_info", []))
            out.append((EVENTS, LANE * (2 if len(d.events) > 1 or sources else 1), None))
        out += [(RULER, 24, None), (DNA, 16, "ref"), (SCORES, 44, "ref")]
        for r in d.rows:
            if r.group == "reference" or (r.group == "alternative" and alternatives) \
                    or (r.group == "wt" and (r.rank == 1 or others)):
                out.append((TRANSCRIPT, 34, r))
        out += [(DNA, 16, "pat"), (SCORES, 44, "pat"), (DELTA, 40, None)]
        for r in d.rows:
            if r.group == "patient" and (r.rank == self.sample_rank or others):
                out.append((TRANSCRIPT, 34, r))
        return out

    def total_height(self):
        return sum(h + ROW_GAP for _, h, _ in self.rows())

    def fit(self):
        self.first, self.cols_per_px = 0.0, max(self.data.n / max(self.width() - 2, 1), 1 / MAX_PX_PER_COLUMN)
        self.update()

    def zoom_to(self, first, last, margin=0.1):
        span = max(last - first + 1, 1)
        pad = span * margin
        self.cols_per_px = max((span + 2 * pad) / max(self.width(), 1), 1 / MAX_PX_PER_COLUMN)
        self.first = max(0.0, first - pad)
        self._clamp()
        self.update()

    def visible(self):
        """(first, last) column in view (of what is drawn)."""
        a = max(0, int(self.first))
        return a, min(self.data.n - 1, int(math.ceil(self.first + self.width() * self.cols_per_px)))

    def col_at(self, x):
        return int(self.first + x * self.cols_per_px)

    def x_of(self, c):
        return (c - self.first) / self.cols_per_px

    def _clamp(self):
        self.cols_per_px = min(max(self.cols_per_px, 1 / MAX_PX_PER_COLUMN),
                               max(self.data.n / max(self.width(), 1), 1 / MAX_PX_PER_COLUMN))
        # Zoom levels are kept on a fine grid, so tiles of one zoom are found again.
        self.cols_per_px = float(f"{self.cols_per_px:.6g}")
        self.first = min(max(0.0, self.first), max(0.0, self.data.n - self.width() * self.cols_per_px))

    # --- tiles ------------------------------------------------------------------------------
    def _layout(self):
        return (id(self.data), self.show_alternatives, self.show_spliceai_alternatives, self.focus, self.sample_rank,
                self._legend[2],
                self.devicePixelRatioF())

    def _clear_tiles(self):
        self._tiles.clear()
        self._pending.clear()
        self._recent_zooms.clear()

    def _render(self, key, data, rows, first, cpp, height):
        """In the drawing thread: one tile."""
        try:
            if key[1] != self.cols_per_px or key[0] != self._layout():
                self._tile_ready.emit(key, None)        # the view moved on: dropped, asked for again if needed
                return
            ratio = key[0][-1]
            image = QImage(math.ceil(TILE * ratio), math.ceil(height * ratio), QImage.Format_ARGB32_Premultiplied)
            image.setDevicePixelRatio(ratio)
            image.fill(QColor(theme.SURFACE))
            p = QPainter(image)
            TrackPainter(data, rows, first, cpp, labels=False).paint(p, TILE)
            p.end()
            self._tile_ready.emit(key, image)
        except RuntimeError:            # the widget was closed meanwhile
            pass

    def _on_tile(self, key, image):
        self._pending.discard(key)
        if image is None or key[0] != self._layout():
            if image is None:
                self.update()           # if that zoom is shown again, its tiles are asked for again
            return
        self._tiles[key] = image
        while len(self._tiles) > MAX_TILES / key[0][-1] ** 2:
            self._tiles.popitem(last=False)
        self.update()

    def _request(self, index):
        key = (self._layout(), self.cols_per_px, index)
        if key in self._tiles:
            self._tiles.move_to_end(key)
            return self._tiles[key]
        if key not in self._pending:
            self._pending.add(key)
            self._pool.submit(self._render, key, self.data, self.rows(), index * TILE * self.cols_per_px,
                              self.cols_per_px, self.total_height())
        return None

    def _stand_in(self, p, index, height):
        """A missing tile: the last zoom's tiles covering it, scaled; else a blank."""
        a = index * TILE * self.cols_per_px
        b = a + TILE * self.cols_per_px
        drew = False
        for zoom in reversed(self._recent_zooms):
            if zoom == self.cols_per_px:
                continue
            span = TILE * zoom
            for j in range(int(a // span), int(b // span) + 1):
                old = self._tiles.get((self._layout(), zoom, j))
                if old is not None:
                    x0, x1 = self.x_of(j * span), self.x_of((j + 1) * span)
                    p.drawImage(QRectF(x0, 0, x1 - x0, height), old)
                    drew = True
            if drew:
                return

    def paintEvent(self, _event):
        t0 = time.perf_counter()
        if not self._fitted and self.width() > 10:
            self._fitted = True
            if self.start_view:
                self.show_genomic(*self.start_view, margin=0.02)
            else:
                self.fit()
        p = QPainter(self)
        p.fillRect(self.rect(), QColor(theme.SURFACE))
        height = self.total_height()
        ratio = self.devicePixelRatioF()
        if self.cols_per_px not in self._recent_zooms:
            self._recent_zooms.append(self.cols_per_px)
            del self._recent_zooms[:-4]
        span = TILE * self.cols_per_px
        t_first = int(self.first // span)
        t_last = int((self.first + self.width() * self.cols_per_px) // span)
        missing = []
        for index in range(t_first, t_last + 1):
            image = self._request(index)
            if image is None:
                missing.append(index)
            else:
                # On a whole device pixel: a tile drawn between pixels would be resampled, and blur.
                p.drawImage(QPointF(round(self.x_of(index * span) * ratio) / ratio, 0), image)
        for index in missing:
            self._stand_in(p, index, height)
        for index in (t_first - 1, t_last + 1):         # a margin each side, made ahead
            if index >= 0 and index * span < self.data.n:
                self._request(index)
        self._overlay(p)
        if self.cursor_col is not None:
            x = self.x_of(self.cursor_col + 0.5)
            p.setPen(QPen(QColor(theme.INK), 1, Qt.DashLine))
            p.drawLine(QPointF(x, 0), QPointF(x, self.height()))
        p.end()
        if (self.first, self.cols_per_px, id(self.data)) != getattr(self, "_drawn_at", None):
            self._drawn_at = (self.first, self.cols_per_px, id(self.data))
            self.moved.emit()
        self.frame_ms.append((time.perf_counter() - t0) * 1000)
        del self.frame_ms[:-200]

    def delta_scale(self):
        """(largest |delta| in view, what the delta row's height stands for: that, at least 0.01)."""
        a, b = self.visible()
        top = self.data.max_abs_delta(a, b)
        return top, max(top, 0.01)

    def _overlay(self, p):
        """What isn't tiled: the legend, the event bands and labels, the delta row and the row names."""
        p.setFont(view_font())
        rows = self.rows()
        top, scale = self.delta_scale()
        self.label_boxes = []
        y = 2
        for kind, h, payload in rows:
            rect = QRectF(0, y, self.width(), h)
            if kind in (LEGEND, EVENTS):
                p.fillRect(rect, QColor(theme.SURFACE))
            if kind == LEGEND:
                legend(p, rect, self.data)
            elif kind == DELTA:
                p.fillRect(rect, QColor(GROUP_BACKGROUNDS["pat"]))
                c0, c1 = self.visible()
                painter = TrackPainter(self.data, rows, self.first, self.cols_per_px)
                painter.background = GROUP_BACKGROUNDS["pat"]
                painter.delta(p, rect, c0, c1, scale)
            y += h + ROW_GAP
        event_bands(p, self.data, rows, self.x_of, self.width(), self.show_origin)
        insertion_gaps(p, self.data, rows, self.x_of, self.width())
        TrackPainter(self.data, rows, self.first, self.cols_per_px).text(p, self.width())
        y = 2
        for kind, h, payload in rows:
            rect = QRectF(0, y, self.width(), h)
            if kind == EVENTS:
                event_labels(p, rect, self.data, self.x_of, self.label_boxes, self.show_origin)
            row_label(p, rect, kind, payload, self.data, menu=True)
            y += h + ROW_GAP

    def tiles_ready(self):
        """Every tile in view is drawn (for tests and the speed report)."""
        span = TILE * self.cols_per_px
        first = int(self.first // span)
        last = int((self.first + self.width() * self.cols_per_px) // span)
        return all((self._layout(), self.cols_per_px, i) in self._tiles for i in range(first, last + 1))

    def paint(self, p, width):
        """Draws every row into painter p, `width` pixels wide, at once (exports)."""
        TrackPainter(self.data, self.rows(), self.first, self.cols_per_px, delta_scale=self.delta_scale()[1]).paint(
            p, width)

    # --- interaction ------------------------------------------------------------------------
    def wheelEvent(self, event):
        steps = event.angleDelta().y() / 120
        if not steps:
            return
        x = event.position().x()
        anchor = self.first + x * self.cols_per_px
        self.cols_per_px *= 0.8 ** steps
        self._clamp()
        self.first = anchor - x * self.cols_per_px
        self._clamp()
        self.update()

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            self._drag = (event.position().x(), self.first)
            self._moved = False

    def mouseReleaseEvent(self, event):
        if self._drag is not None and not self._moved and self._pick_sample(event):
            self._drag = None
            return
        if self._drag is not None and not self._moved:
            c = min(max(self.col_at(event.position().x()), 0), self.data.n - 1)
            self.details.emit(self.full_details(c, event.position().y()))
        self._drag = None

    def mouseDoubleClickEvent(self, event):
        """A gentle zoom in (x2) around the pointer; the whole gene is the corner button's."""
        if event.button() == Qt.LeftButton:
            x = event.position().x()
            anchor = self.first + x * self.cols_per_px
            self.cols_per_px /= 2
            self._clamp()
            self.first = anchor - x * self.cols_per_px
            self._clamp()
            self.update()

    def is_whole(self):
        return self.first <= 0 and self.first + self.width() * self.cols_per_px >= self.data.n - 1

    def toggle_whole(self):
        """The whole gene, or back to where the view was before it."""
        if self.is_whole() and self._before_whole:
            self.first, self.cols_per_px = self._before_whole
            self._before_whole = None
            self._clamp()
            self.update()
        else:
            self._before_whole = (self.first, self.cols_per_px)
            self.fit()

    def patients(self):
        return [r for r in self.data.rows if r.group == "patient"]

    def set_sample(self, rank):
        if rank != self.sample_rank:
            self.sample_rank = rank
            self.update()
            self.sample_changed.emit(rank)

    def _pick_sample(self, event):
        """A click on a sample candidate's name: a list of every candidate. True if it was one."""
        kind, payload, _ = self._row_at(event.position().y())
        if kind != TRANSCRIPT or payload.group != "patient" or len(self.patients()) < 2:
            return False
        if event.position().x() > QFontMetrics(view_font()).horizontalAdvance(payload.name) + 160:
            return False
        menu = QMenu(self)
        for r in self.patients():
            action = menu.addAction(r.name + (f"  ({r.relative:.3f})" if r.relative is not None else ""))
            action.setCheckable(True)
            action.setChecked(r.rank == self.sample_rank)
            action.triggered.connect(lambda _=False, k=r.rank: self.set_sample(k))
        menu.exec(event.globalPosition().toPoint())
        return True

    def mouseMoveEvent(self, event):
        x, y = event.position().x(), event.position().y()
        if self._drag is not None:
            if abs(x - self._drag[0]) > 2:
                self._moved = True
            self.first = self._drag[1] - (x - self._drag[0]) * self.cols_per_px
            self._clamp()
        c = min(max(self.col_at(x), 0), self.data.n - 1)
        self.cursor_col = c
        self.hovered.emit(self.to_genomic(c))
        QToolTip.showText(event.globalPosition().toPoint(), self.tooltip(c, y), self)
        if event.modifiers() & Qt.ShiftModifier:
            self.details.emit(self.full_details(c, y))
        self.update()

    def _row_at(self, y):
        top = 2
        for kind, h, payload in self.rows():
            if top <= y < top + h + ROW_GAP:
                return kind, payload, y - top
            top += h + ROW_GAP
        return None, None, 0

    def tooltip(self, c, y):
        """Spec 15 item 3: at most 4 short lines -- an exon's if the pointer is on an exon block, else the
        position (genomic and c.), the bases and the values of the row under the pointer."""
        d = self.data
        c = min(max(c, 0), d.n - 1)
        kind, payload, dy = self._row_at(y)
        if kind == EVENTS:
            info = self.event_at_label(self.x_of(c), y)
            if info is not None:
                return "\n".join(x for x in (info["label"], info["origin"] and t("origin: {origin}",
                                                                                  origin=info["origin"]),
                                              info["effect"]) if x)
        if kind == TRANSCRIPT and 12 <= dy <= 26:
            blk = next((b for b in payload.blocks if b.first <= c <= b.last), None)
            if blk is not None:
                same = t("reference exon {n}", n=blk.identity) if blk.identity else t("no reference exon (novel)")
                same = t("duplicate of {what}", what=same) if blk.duplicate else same
                return (t("{row}, exon {n} ({same})", row=payload.name, n=blk.number, same=same)
                        + f"\n{blk.end - blk.start + 1:,} nt" + (f", {t(blk.status)}" if blk.status else ""))
        f = lambda v: "--" if np.isnan(v) else f"{v:.3f}"
        lines = [t(d.position_text(c)), t("reference {ref}, sample {sample}", ref=d.ref_base[c], sample=d.pat_base[c])]
        tracks = d.tracks
        if kind == SCORES:
            which = "ref" if payload == "ref" else "pat"
            lines.append(t("reference: acceptor {acc}, donor {don}", acc=f(tracks["ref_acc"][c]),
                           don=f(tracks["ref_don"][c])) if which == "ref" else
                         t("sample: acceptor {acc}, donor {don}", acc=f(tracks["pat_acc"][c]),
                           don=f(tracks["pat_don"][c])))
        elif kind == DELTA:
            lines.append(t("delta: acceptor {acc}, donor {don}", acc=f(tracks["delta_acc"][c]),
                           don=f(tracks["delta_don"][c])))
        elif kind == TRANSCRIPT:
            m = int(payload.membership_range(c, c)[0])
            lines.append(f"{payload.name}: " + (t("exon {n}", n=m) if m > 0 else t("intron") if m < 0 else
                                                 t("outside")))
        else:
            lines.append(t("sample: acceptor {acc}, donor {don}", acc=f(tracks["pat_acc"][c]),
                           don=f(tracks["pat_don"][c])))
        info = self.event_at(c)
        if info is not None:        # spec 16 item 2: inside an event, what it is and where its bases come from
            name = info["label"].split(" · ")[0]
            lines = lines[:2] + [f"{name}, origin: {info['origin']}" if info["origin"] else name] + lines[2:3]
        return "\n".join(lines[:4])

    def event_at(self, c):
        """The info (cdna_analysis.event_summary) of the event over column c, else None."""
        return next((info for (a, b, _), info in zip(self.data.events, getattr(self.data, "event_info", []))
                     if a <= c <= b), None)

    def event_at_label(self, x, y):
        return next((info for rect, info in self.label_boxes if rect.contains(QPointF(x, y))), None)

    def full_details(self, c, y):
        """Everything known about the base (and the exon or event) under the pointer, for the side panel:
        [(section, [(key, value)])] (spec 16 item 3)."""
        d = self.data
        c = min(max(c, 0), d.n - 1)
        sections = d.details(c)
        kind, payload, dy = self._row_at(y)
        if kind == TRANSCRIPT:
            blk = next((b for b in payload.blocks if b.first <= c <= b.last), None)
            if blk is not None:
                acc, don = d.tracks["pat_acc" if payload.on_haplotype else "ref_acc"], \
                    d.tracks["pat_don" if payload.on_haplotype else "ref_don"]
                first, last = (blk.first, blk.last) if rightwards(d) else (blk.last, blk.first)
                f = lambda v: "--" if np.isnan(v) else f"{v:.3f}"
                same = (t("duplicate of reference exon {n}", n=blk.identity) if blk.duplicate else
                        t("reference exon {n}", n=blk.identity) if blk.identity else t("novel"))
                pairs = [(N_("Transcript"), payload.name),
                         (N_("Its exon"), f"{blk.number} ({same})"),
                         (N_("On the haplotype") if payload.on_haplotype else N_("On the reference"),
                          f"{blk.start:,}..{blk.end:,} ({blk.end - blk.start + 1:,} nt)"),
                         (N_("Acceptor at start"), f(acc[first])), (N_("Donor at end"), f(don[last]))]
                if blk.status:
                    pairs.append((N_("vs the reference"), blk.status))
                sections.insert(0, (N_("Exon"), pairs))
        info = self.event_at_label(self.x_of(c), y) if kind == EVENTS else self.event_at(c)
        if info is not None:
            sections.insert(0, (N_("Event"), event_pairs(info)))
        return sections

    def leaveEvent(self, _event):
        self.cursor_col = None
        self.update()

    def resizeEvent(self, _event):
        self._clamp()
        if self.minimumHeight() != self.total_height():     # the legend flows to the width
            self.setMinimumHeight(self.total_height())

    def closeEvent(self, event):
        self._pool.shutdown(wait=False, cancel_futures=True)
        super().closeEvent(event)

    # --- export -------------------------------------------------------------------------------
    def _footer(self, p, width):
        """The research-use note under an exported image, which travels without the startup pop-up (spec 14 C7)."""
        p.setPen(QColor(theme.INK_2))
        arrow = "5' ──▶ 3'" if rightwards(self.data) else "3' ◀── 5'"
        p.drawText(QRectF(4, self.total_height(), width - 8, FOOTER), Qt.AlignLeft | Qt.AlignVCenter | Qt.TextWordWrap,
                   f"{orientation_text(self.data)} {arrow}. {cdna_analysis.RESEARCH_NOTE}")

    def render_png(self, path, width=None, pixel_ratio=2):
        """The view as a PNG, `width` logical pixels wide, drawn at pixel_ratio image pixels per logical pixel."""
        width = width or self.width()
        image = QImage(int(width * pixel_ratio), int((self.total_height() + FOOTER) * pixel_ratio), QImage.Format_ARGB32)
        image.setDevicePixelRatio(pixel_ratio)
        image.fill(QColor(theme.SURFACE))
        p = QPainter(image)
        self.paint(p, width)
        self._footer(p, width)
        p.end()
        return image.save(path)

    def render_svg(self, path, width=None):
        from PySide6.QtCore import QSize
        from PySide6.QtSvg import QSvgGenerator

        width = width or self.width()
        gen = QSvgGenerator()
        gen.setFileName(path)
        gen.setSize(QSize(width, self.total_height() + FOOTER))
        gen.setViewBox(QRectF(0, 0, width, self.total_height() + FOOTER))
        gen.setTitle(f"{self.data.gene} {self.data.build}")
        gen.setDescription(cdna_analysis.RESEARCH_NOTE)
        p = QPainter(gen)
        self.paint(p, width)
        self._footer(p, width)
        p.end()
        return True
