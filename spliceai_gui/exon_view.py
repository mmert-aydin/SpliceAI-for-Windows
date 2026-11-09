"""The viewer's exon schematic (spec 8.2 B, spec 15 items 1 and 9): every
transcript as its exons on one shared axis -- exons to scale, introns
compressed -- so the same exon lines up in every row, coloured by the
reference exon each one matches (novel exons grey and hatched).

- The gene's strand: chevrons on the intron lines and a 5'->3' arrow. In
  transcript orientation (the default) a minus-strand gene is mirrored, exon
  1 on the left; the header says which orientation is shown.
- The applied events, at their true place (in an exon or an intron), with
  their type and size ("INV 70 bp"); a click jumps to them in the stacked
  view. A "!" marks an event in or next to a region where the WT prediction
  differs from the reference (cdna_analysis.local_concordance, as in the
  results).
- Between the reference and the WT prediction, and between the WT
  prediction and each sample candidate, what differs, at the exon or
  junction it concerns ("exon 8 skipped", "cryptic donor +37", "pseudoexon
  72 nt in intron 1"), or "no change".
- Every exon and intron where the WT prediction differs from the reference
  is shaded (cdna_analysis.wt_regions).

Hovering shows a short tooltip; a click sends the full details to the
viewer's side panel and jumps the stacked view there.
"""
import bisect
import re

from PySide6.QtCore import QRectF, QSize, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QImage, QPainter, QPen
from PySide6.QtWidgets import QToolTip, QWidget

from spliceai_pipeline import cdna_analysis

from . import theme
from .i18n import t
from .stacked_view import FOOTER, event_color, exon_color, view_font

NAME_W = 210
GAP = 22                # px of a compressed intron
HEADER, EVENTS, NOTE, LINE, SPACE = 18, 34, 15, 20, 6
CAUTION = t("WT prediction differs from the reference in this region; interpret the sample prediction here with "
           "caution.")
ABBREVIATION = cdna_analysis.ABBREVIATION


class Axis:
    """Reference positions -> x: exons (merged over every row) to scale, the gaps between them compressed to
    at most GAP pixels; mirrored for transcript orientation. room {position: bases}: space kept at a position for
    exons in inserted bases, to scale too, in .room as (left, right)."""

    def __init__(self, intervals, x0, x1, mirror, room=None):
        segs = []
        for a, b in sorted(intervals):
            if segs and a <= segs[-1][1] + 1:
                segs[-1][1] = max(segs[-1][1], b)
            else:
                segs.append([a, b])
        self.x0, self.x1, self.mirror = x0, x1, mirror
        width = max(x1 - x0, 50)
        room = room or {}
        exonic = sum(b - a + 1 for a, b in segs) + sum(room.values()) or 1
        gaps = [segs[i + 1][0] - segs[i][1] - 1 for i in range(len(segs) - 1)]
        most = min(GAP, 0.5 * width / max(len(gaps), 1))       # in a narrow view the introns give way, never the width
        scale = (width - most * len(gaps)) / exonic
        for _ in range(3):
            gap_px = [min(most, g * scale) for g in gaps]
            scale = (width - sum(gap_px)) / exonic
        gap_px = [min(most, g * scale) for g in gaps]
        self.scale = scale
        self.room = {}
        self.points = []                # (position, x) breakpoints of a piecewise-linear map
        x = x0
        for i, (a, b) in enumerate(segs):
            self.points.append((a, x))
            at = a
            for p in sorted(p for p in room if a <= p <= b + 1):
                x += (p - at) * scale
                self.points.append((p, x))
                left, x = x, x + room[p] * scale
                self.points.append((p, x))
                self.room[p] = (left, x) if not mirror else (x0 + x1 - x, x0 + x1 - left)
                at = p
            x += (b + 1 - at) * scale
            self.points.append((b + 1, x))
            if i < len(gaps):
                x += gap_px[i]
        self._pos = [p for p, _ in self.points]

    def x(self, pos, near=False):
        """near: at a position with room, the side before the room, where an exon ending there stops."""
        pts = self.points
        if not pts:
            return self.x0
        i = bisect.bisect_right(self._pos, pos) - 1
        if near and i >= 0 and self._pos[i] == pos:
            i = bisect.bisect_left(self._pos, pos)
        if i < 0:
            x = pts[0][1] - min(GAP, (pts[0][0] - pos) * 0.01)
        elif i >= len(pts) - 1:
            x = pts[-1][1] + min(GAP, (pos - pts[-1][0]) * 0.01)
        else:
            (p0, x0), (p1, x1) = pts[i], pts[i + 1]
            x = x0 + (x1 - x0) * (pos - p0) / max(p1 - p0, 1)
        return self.x0 + self.x1 - x if self.mirror else x


def place_notes(fm, marks, lo, hi):
    """[(left, width, text, tooltip)]: the notes in a row between x = lo and hi, each at its mark as far as they
    allow, none over another or past the edge. What doesn't fit becomes "+n more" (the rest in its tooltip)."""
    marks = sorted(marks, key=lambda m: m[0] if m[0] is not None else lo)
    width = lambda text: fm.horizontalAdvance("▼ " + text) + 6
    for keep in range(len(marks), 0, -1):
        rest = [text for _, text in marks[keep:]]
        items = [(lo + 4 if x is None else x, text, text) for x, text in marks[:keep]]
        if rest:
            items.append((items[-1][0], f"+{len(rest)} more", "\n".join(rest)))
        out, edge = [], lo - 4
        for x, text, tip in items:
            left = max(x - 4, edge + 4)
            out.append([left, width(text), text, tip])
            edge = left + width(text)
        edge = hi
        for item in reversed(out):                  # pulled back from the right edge
            item[0] = min(item[0], edge - item[1])
            edge = item[0] - 4
        if out[0][0] >= lo:
            return out
    return []


def _identity(e, canonical):
    return next((k for k, c in enumerate(canonical, 1) if e and e[0] <= c[1] and c[0] <= e[1]), 0)


class ExonView(QWidget):
    exon_clicked = Signal(int, int)     # first and last (genomic) column of the exon, in the stacked view
    event_clicked = Signal(int)         # the index of the event in data.events
    details = Signal(str)               # everything about what was clicked, for the side panel

    def __init__(self, data, parent=None):
        super().__init__(parent)
        self.data = data
        self.result = getattr(data, "result", None)
        self.show_alternatives = False
        self.show_spliceai_alternatives = False
        self.transcript_orientation = True
        self._hit = []                  # (rect, kind, payload)
        self.setMouseTracking(True)
        self.setMinimumHeight(200)

    # --- what is drawn --------------------------------------------------------------------------
    def mirrored(self):
        return self.transcript_orientation and self.data.strand == "-"

    def rightwards(self):
        return (self.data.strand == "+") != self.mirrored()

    def entries(self):
        """([dict] per row: name, group, exons (reference coordinates; None in inserted sequence), the
        transcript (cdna_analysis.TranscriptResult) or None and what it is compared with; the canonical exons
        in transcript order)."""
        r = self.result
        c = r.gene.canonical
        canonical = list(c.exons) if c.strand == "+" else list(reversed(c.exons))
        by = {(t.kind, t.rank): t for t in r.transcripts}
        out = [{"name": f"Reference ({c.accession})", "group": "reference", "t": by[("Reference", 0)],
                "exons": canonical, "compare": None, "cds": c.cds}]
        if self.show_alternatives:
            for alt in r.gene.transcripts[1:]:
                exons = list(alt.exons) if alt.strand == "+" else list(reversed(alt.exons))
                out.append({"name": alt.accession, "group": "alternative", "t": None, "exons": exons,
                            "compare": None, "cds": alt.cds})
        for t in r.transcripts:
            if t.kind == "WT predicted" and (t.rank == 1 or self.show_spliceai_alternatives):
                out.append({"name": f"WT predicted {t.rank}", "group": "wt", "t": t, "exons": t.ref_exons,
                            "compare": "Reference" if t.rank == 1 else None})
        for t in r.transcripts:
            if t.kind == "Patient":
                out.append({"name": f"Sample {t.rank}", "group": "patient", "t": t, "exons": t.ref_exons,
                            "compare": "WT"})
        for e in out:
            e["identities"] = [_identity(x, canonical) for x in e["exons"]]
            e["copies"] = cdna_analysis.in_insertion(r, e["t"]) if e["group"] == "patient" else {}
        return out, canonical

    def content_height(self):
        rows, _ = self.entries() if self.result else ([], [])
        return HEADER + EVENTS + sum(LINE + SPACE + (NOTE if e["compare"] else 0) for e in rows) + 10

    def sizeHint(self):
        return QSize(900, self.content_height())

    def paintEvent(self, _event):
        p = QPainter(self)
        p.fillRect(self.rect(), QColor(theme.SURFACE))
        if self.result is not None:
            self.paint(p, self.width())
        p.end()

    def axis(self, rows=None, width=None):
        rows = rows if rows is not None else self.entries()[0]
        intervals = [tuple(e) for row in rows for e in row["exons"] if e]
        intervals += [(min(x.start, x.end), max(x.start, x.end)) for x in self.result.edits]
        room = {}
        for row in rows:
            used = {}
            for _copy, at, n in row["copies"].values():
                used[at] = used.get(at, 0) + n
            for at, n in used.items():
                room[at] = max(room.get(at, 0), n)
        return Axis(intervals, NAME_W, (width or self.width()) - 12, self.mirrored(), room)

    def paint(self, p, width):
        rows, canonical = self.entries()
        axis = self.axis(rows, width)
        self._hit = []
        self._width = width
        self.notes = {}                 # row name -> what its arrows say (for tests and the side panel)
        p.setFont(view_font())
        self._header(p, width)
        bottom = self.content_height() - 10
        self._shading(p, axis, canonical, HEADER + EVENTS, bottom)
        self._events(p, axis, HEADER)
        for x0, x1, color in self.event_spans:      # each event as a band through every row, behind the exons
            band = QColor(color)
            band.setAlpha(55)
            top = HEADER + EVENTS - 3       # under the labels, not through them
            p.fillRect(QRectF(x0 - 1, top, max(x1 - x0, 1) + 2, bottom - top), band)
            p.setPen(QPen(color, 1, Qt.DashLine))
            p.drawLine(int(x0 - 1), top, int(x0 - 1), int(bottom))
            p.drawLine(int(x1 + 1), top, int(x1 + 1), int(bottom))
        y = HEADER + EVENTS
        wt1 = next((e for e in rows if e["group"] == "wt" and e["t"].rank == 1), None)
        for row in rows:
            if row["compare"]:
                base = canonical if row["compare"] == "Reference" else (wt1["exons"] if wt1 else canonical)
                self._notes(p, axis, row, base, canonical, wt1, y)
                y += NOTE
            self._row(p, axis, row, y)
            y += LINE + SPACE
        if self.minimumHeight() != int(y + 10):
            self.setMinimumHeight(int(y + 10))

    def _header(self, p, width):
        strand = t("+ strand") if self.data.strand == "+" else t("− strand")
        if self.mirrored():
            text = t("Transcript orientation (5'→3'), {strand}: exon 1 on the left", strand=strand)
        elif self.transcript_orientation:
            text = t("Transcript orientation (5'→3'), {strand}", strand=strand)
        elif self.data.strand == "+":
            text = t("Genomic orientation, {strand}", strand=strand)
        else:
            text = t("Genomic orientation, {strand}: exon 1 on the right", strand=strand)
        arrow = "5' ──▶ 3'" if self.rightwards() else "3' ◀── 5'"
        self.header_text = f"{text}    {arrow}"
        p.setPen(QColor(theme.INK))
        p.drawText(QRectF(4, 0, width - 8, HEADER), Qt.AlignLeft | Qt.AlignVCenter, self.header_text)

    def _shading(self, p, axis, canonical, top, bottom):
        """Spec 15 item 9: where the WT prediction differs from the reference (the regions C2 reads)."""
        regions = cdna_analysis.wt_regions(self.result)
        color = QColor(theme.WARN_SOFT)
        color.setAlpha(170)
        self.shaded = []
        for (kind, k), what in sorted(regions.items()):
            if kind == "exon":
                a, b = canonical[k - 1]
            else:
                e1, e2 = canonical[k - 1], canonical[k]
                a, b = min(e1[1], e2[1]) + 1, max(e1[0], e2[0]) - 1
            x0, x1 = sorted((axis.x(a), axis.x(b + 1)))
            rect = QRectF(x0 - 2, top, max(x1 - x0, 2) + 4, bottom - top)
            p.fillRect(rect, color)
            self.shaded.append((kind, k))
            sites = cdna_analysis.region_sites(self.result, kind, k)
            scores = ("\n" + t("SpliceAI (WT) at the reference's splice sites: {sites} (site threshold {threshold:g})",
                                sites=", ".join(f"{t(name)} {p:.2f}" for name, p in sites),
                                threshold=self.result.settings.threshold)) if sites else ""
            self._hit.append((rect, "region", t("WT prediction differs from the reference here: {kind} {k} ({what})",
                                                kind=t(kind), k=k, what=t(what)) + scores))

    def _events(self, p, axis, top):
        local = cdna_analysis.local_concordance(self.result)
        labels = [e for e in self.result.compat.alignment.events if e.apply]
        self.event_marks = []
        self.event_spans = []           # (x0, x1, colour): drawn through every row once the rows are drawn
        named, marks, spots = {}, [], []
        for i, e in enumerate(self.result.edits):
            lo, hi = min(e.start, e.end), max(e.start, e.end)
            insertion = e.end < e.start
            x0, x1 = sorted((axis.x(lo), axis.x(hi + 1)))
            size = max(e.end - e.start + 1, len(e.alt))
            text = f"{ABBREVIATION.get(e.kind, e.kind.upper())} {size:,} bp"
            color = event_color(e.kind)      # spec 16: the event type's colour, as in the stacked view
            p.setPen(QPen(color, 2))
            if insertion or x1 - x0 < 3:
                x = axis.x(e.start - 0.5 if insertion else lo)
                p.drawLine(int(x), top + 14, int(x), top + EVENTS - 2)
                p.drawText(QRectF(x - 4, top + 4, 8, 10), Qt.AlignCenter, "▼")
                x0 = x1 = x
            else:
                p.fillRect(QRectF(x0, top + 16, x1 - x0, 6), color)
            self.event_spans.append((x0, x1, color))
            caution = bool(local[i]["differing"] or local[i]["introns"]) if i < len(local) else False
            name = labels[i].label if i < len(labels) and labels[i].label else e.kind
            # spec 16: the event's full description, as everywhere else
            info = self.data.event_info[i]["line"] if i < len(getattr(self.data, "event_info", [])) else \
                f"at {e.start:,} ({size:,} bp)"
            tip = f"{name}: {info}" + (f"\n! {CAUTION}" if caution else "")
            self.event_marks.append(((x0 + x1) / 2, text, caution))
            spots.append((QRectF(x0 - 4, top, x1 - x0 + 8, 22), "event", (i, tip)))
            named.setdefault(text + ("  !" if caution else ""), []).append((color, i, tip))
            marks.append((x0 + 7 if x0 == x1 else x0, text + ("  !" if caution else "")))    # beside a line, not over it
        for left, w, text, tip in place_notes(p.fontMetrics(), marks, NAME_W, self._width - 6):
            if text in named:
                color, i, tip = named[text].pop(0)
                p.setPen(color.darker(140))
                self._hit.append((QRectF(left, top + 20, w, 16), "event", (i, tip)))
            else:                       # "+n more": they are in its tooltip
                p.setPen(QColor(theme.INK_2))
                self._hit.append((QRectF(left, top + 20, w, 16), "note", tip))
            p.drawText(QRectF(left, top + 22, w, 12), Qt.AlignLeft | Qt.AlignVCenter, text)
        self._hit += spots              # last in, first found: a mark before a label moved over it

    def _notes(self, p, axis, row, base, canonical, wt1, y):
        """Small arrows over the row: what differs from the row it is compared with."""
        marks = self.differences(axis, row, base, canonical, wt1)
        self.notes[row["name"]] = [text for _, text in marks] or ["no change"]
        if not marks:
            p.setPen(QColor(theme.INK_3))
            text = "▼ " + (t("no change vs the reference") if row["compare"] == "Reference" else
                           t("no change vs WT predicted 1"))
            p.drawText(QRectF(NAME_W, y, 320, NOTE), Qt.AlignLeft | Qt.AlignVCenter, text)
            return
        p.setPen(QColor(theme.WARN))
        for left, w, text, tip in place_notes(p.fontMetrics(), marks, NAME_W, self._width - 6):
            p.drawText(QRectF(left, y, w, NOTE), Qt.AlignLeft | Qt.AlignVCenter, "▼ " + text)
            self._hit.append((QRectF(left, y, w, NOTE), "note", tip))

    def differences(self, axis, row, base, canonical, wt1):
        """[(x or None, text)] of what differs between row and what it is compared with."""
        t = row["t"]
        rows = t.vs_ref[1] if row["compare"] == "Reference" else t.vs_wt[1]
        labels = (cdna_analysis.exon_labels(wt1["exons"], canonical) if row["compare"] == "WT" and wt1 else
                  [f"Exon {k}" for k in range(1, len(canonical) + 1)])
        plus = self.data.strand == "+"
        marks = {}          # as a dict: both exons of a retained intron name it, once is enough
        add = lambda x, text: marks.setdefault((x, text))
        for r in rows:
            status = r["status"]
            if status == "identical":
                continue
            if status == "missing":
                a, b = base[r["ref"] - 1]
                label = labels[r["ref"] - 1]
                add((axis.x(a) + axis.x(b + 1)) / 2, f"{label[0].lower() + label[1:]} skipped")
                continue
            if r["ref"] is None:
                e = row["exons"][r["other"] - 1]
                if e is None:
                    copy = row["copies"][r["other"]][0]
                    x0, x1 = self.inserted_boxes(axis, row)[r["other"]]
                    add(x0 if x1 - x0 < 8 else (x0 + x1) / 2, f"exon {copy} (dup)" if copy else "an exon in inserted sequence")
                    continue
                k = cdna_analysis.intron_of(canonical, min(e), max(e))
                add((axis.x(e[0]) + axis.x(e[1] + 1, near=True)) / 2,
                    f"{'pseudoexon' if k else 'novel exon'} {e[1] - e[0] + 1:,} nt" + (f" in intron {k}" if k else ""))
                continue
            e = row["exons"][r["other"] - 1] if r["other"] else None
            for part in status.split(", "):
                m = re.match(r"(start|end) shifted ([+-]\d+)", part)
                if m and e:
                    first = (e[0] if plus else e[1]) if m.group(1) == "start" else (e[1] if plus else e[0])
                    kind = "acceptor" if m.group(1) == "start" else "donor"
                    add(axis.x(first + 0.5), f"cryptic {kind} {m.group(2)}")
            m = re.search(r"one exon with reference exon (\d+)", status)
            if m:
                k = min(int(m.group(1)), r["ref"])
                e1, e2 = canonical[k - 1], canonical[k]
                add((axis.x(min(e1[1], e2[1])) + axis.x(max(e1[0], e2[0]))) / 2, f"intron {k} retained")
        return list(marks)

    def inserted_boxes(self, axis, row):
        """{exon number: (left, right)} of a row's exons in inserted bases, in the room kept at their insertion, in
        transcript order."""
        out, used = {}, {}
        for number, (_copy, at, n) in sorted(row["copies"].items()):
            left, right = axis.room[at]
            w, before = n * axis.scale, used.get(at, 0)
            out[number] = (left + before, left + before + w) if self.rightwards() else (right - before - w, right - before)
            used[at] = before + w
        return out

    def _row(self, p, axis, row, y):
        tr_ = row["t"]
        label = row["name"]
        if tr_ is not None and tr_.relative is not None:
            weak = "weak evidence" in cdna_analysis.evidence(tr_, self.result.settings)
            label += f"  ({tr_.relative:.3f}{', weak' if weak else ''})"
        p.setPen(QColor(theme.INK))
        p.drawText(QRectF(4, y, NAME_W - 8, LINE), Qt.AlignLeft | Qt.AlignVCenter, label)
        placed = [e for e in row["exons"] if e]
        cy = y + LINE / 2
        if placed:
            x0, x1 = sorted((axis.x(min(a for a, _ in placed)), axis.x(max(b for _, b in placed) + 1)))
            p.setPen(QPen(QColor(theme.INK_2), 1))
            p.drawLine(int(x0), int(cy), int(x1), int(cy))
            x = x0 + 20
            while x < x1:
                d = 3 if self.rightwards() else -3
                p.drawLine(int(x - d), int(cy - 3), int(x), int(cy))
                p.drawLine(int(x - d), int(cy + 3), int(x), int(cy))
                x += 40
        for number, (e, ident) in enumerate(zip(row["exons"], row["identities"]), 1):
            if e is None:
                continue
            a, b = sorted((axis.x(e[0]), axis.x(e[1] + 1, near=True)))
            rect = QRectF(a, y + 3, max(b - a, 3), LINE - 6)
            cds = row.get("cds")
            if cds:         # an annotated transcript: UTR thin, coding part full height
                p.fillRect(rect.adjusted(0, 4, 0, -4), exon_color(ident))
                lo, hi = max(e[0], cds[0]), min(e[1], cds[1])
                if lo <= hi:
                    ca, cb = sorted((axis.x(lo), axis.x(hi + 1, near=True)))
                    p.fillRect(QRectF(ca, y + 3, max(cb - ca, 2), LINE - 6), exon_color(ident))
            else:
                p.fillRect(rect, exon_color(ident))
            if not ident:
                p.fillRect(rect, QBrush(QColor("white"), Qt.FDiagPattern))
            text = str(ident) if ident else "new"
            if p.fontMetrics().horizontalAdvance(text) + 4 <= rect.width():
                # Dark on the light exon colours, as in the stacked view; centred on a whole pixel, or it blurs.
                p.setPen(QColor(theme.INK) if exon_color(ident).lightness() > 130 else QColor("white"))
                p.drawText(QRectF(round(rect.center().x()) - 40, y + 3, 80, LINE - 6), Qt.AlignCenter, text)
            self._hit.append((rect, "exon", (row, number, e, ident)))
        for number, (a, b) in self.inserted_boxes(axis, row).items():
            copy, _at, n = row["copies"][number]
            rect = QRectF(a, y + 3, max(b - a, 3), LINE - 6)
            p.fillRect(rect, exon_color(copy))
            if not copy:
                p.fillRect(rect, QBrush(QColor("white"), Qt.FDiagPattern))
            text = str(copy) if copy else "new"
            if p.fontMetrics().horizontalAdvance(text) + 4 <= rect.width():
                p.setPen(QColor(theme.INK) if exon_color(copy).lightness() > 130 else QColor("white"))
                p.drawText(QRectF(round(rect.center().x()) - 40, y + 3, 80, LINE - 6), Qt.AlignCenter, text)
            what = (t("exon {n} (dup): a copy of reference exon {n}", n=copy) if copy else
                    t("an exon in inserted sequence"))
            self._hit.append((rect, "note", t("{row}: {what}, {n:,} nt, in the inserted bases", row=row["name"],
                                              what=what, n=n)))

    def render_png(self, path, width=None, pixel_ratio=2):
        """The schematic as a PNG, `width` logical pixels wide, at pixel_ratio image pixels per logical pixel, with
        the research-use note under it (spec 14 C7)."""
        width = width or max(self.width(), 900)
        height = self.content_height() + FOOTER
        image = QImage(int(width * pixel_ratio), int(height * pixel_ratio), QImage.Format_ARGB32)
        image.setDevicePixelRatio(pixel_ratio)
        image.fill(QColor(theme.SURFACE))
        p = QPainter(image)
        self.paint(p, width)
        p.setPen(QColor(theme.INK_2))
        p.drawText(QRectF(4, height - FOOTER, width - 8, FOOTER), Qt.AlignLeft | Qt.AlignVCenter | Qt.TextWordWrap,
                   cdna_analysis.RESEARCH_NOTE)
        p.end()
        return image.save(path)

    # --- interaction -------------------------------------------------------------------------
    def _at(self, pos):
        for rect, kind, payload in reversed(self._hit):
            if rect.contains(pos):
                return kind, payload
        return None, None

    def tooltip_at(self, pos):
        """Spec 15 item 3: short -- an exon's number, length and status, an event's type and size."""
        kind, payload = self._at(pos)
        if kind == "exon":
            row, number, e, ident = payload
            same = t("reference exon {n}", n=ident) if ident else t("no reference exon (novel)")
            status = self._status(row, number)
            return (t("{row}, exon {n} ({same})", row=row["name"], n=number, same=same) + f"\n{e[1] - e[0] + 1:,} nt"
                    + (f", {t(status)}" if status else ""))
        if kind == "event":
            return payload[1]
        if kind in ("note", "region"):
            return payload
        return ""

    def _status(self, row, number):
        t = row["t"]
        if t is None:
            return ""
        return next((r["status"] for r in t.vs_ref[1] if r["other"] == number), "")

    def mouseMoveEvent(self, event):
        text = self.tooltip_at(event.position())
        if text:
            QToolTip.showText(event.globalPosition().toPoint(), text, self)
        else:
            QToolTip.hideText()

    def mousePressEvent(self, event):
        self.click(event.position())

    def click(self, pos):
        kind, payload = self._at(pos)
        if kind == "event":
            self.event_clicked.emit(payload[0])
            self.details.emit(payload[1])
        elif kind == "exon":
            row, number, e, ident = payload
            self.details.emit(self._exon_details(row, number, e, ident))
            block = self._block(row, number)
            if block is not None:
                self.exon_clicked.emit(block.first, block.last)
        elif kind in ("note", "region"):
            self.details.emit(payload)

    def _block(self, row, number):
        view_row = next((r for r in self.data.rows if r.name == row["name"]
                         or (row["group"] == "reference" and r.group == "reference")), None)
        return next((b for b in view_row.blocks if b.number == number), None) if view_row else None

    def _exon_details(self, row, number, e, ident):
        tr_ = row["t"]
        lines = [t("{row}: exon {n} of {total}", row=row["name"], n=number, total=len(row["exons"])),
                 t("reference {where} ({n:,} nt)", where=f"{self.data.chrom}:{e[0]:,}-{e[1]:,}", n=e[1] - e[0] + 1),
                 t("matches reference exon {n}", n=ident) if ident else t("no reference exon (novel)")]
        if tr_ is not None:
            if tr_.kind == "Patient":
                h = tr_.exons[number - 1]
                lines.append(t("haplotype {start:,}-{end:,}", start=h[0], end=h[1]))
            vs_ref = next((r["status"] for r in tr_.vs_ref[1] if r["other"] == number), "")
            if vs_ref:
                lines.append(t("vs the reference ({accession}): {status}",
                               accession=self.result.gene.canonical.accession, status=t(vs_ref)))
            if tr_.kind != "Reference":
                vs_wt = next((r["status"] for r in tr_.vs_wt[1] if r["other"] == number), "")
                if vs_wt:
                    lines.append(t("vs WT predicted 1: {status}", status=t(vs_wt)))
                lines.append(t("change vs the reference: {change}", change=tr_.change_vs_ref))
        return "\n".join(lines)
