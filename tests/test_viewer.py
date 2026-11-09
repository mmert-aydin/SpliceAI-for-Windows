"""Phase 8: the viewer (spec 8.2). The shared columns from the coordinate
map (haplotype.columns), what the viewer shows computed without Qt
(view_data.py: tracks, exons, the completeness rule, hover text, CSV out
and back), and the viewer window offscreen: the three views, PNG / SVG,
navigation, a results file opened with neither SpliceAI nor a FASTA nor the
database, and drawing a 2.3 Mb gene.

    .venv\\Scripts\\python tests\\test_viewer.py
"""
import json
import os
import random
import subprocess
import sys
import tempfile
import time

import helpers
from helpers import check

helpers.watchdog(1200)
helpers.isolate_home()          # the viewer keeps its side panel's width in the settings
os.environ["QT_QPA_PLATFORM"] = "offscreen"
helpers.setup()

import numpy as np  # noqa: E402

from spliceai_pipeline import bundle, cdna_analysis as ca, haplotype, refdb  # noqa: E402
from spliceai_pipeline.haplotype import Edit  # noqa: E402
from spliceai_pipeline.view_data import ViewData  # noqa: E402

rng = random.Random(21)
rand_seq = lambda n: "".join(rng.choice("ACGT") for _ in range(n))
tmp = tempfile.mkdtemp(prefix="spliceai_viewer_")

# --- the shared columns -------------------------------------------------------------------------
REF = rand_seq(2000)
START = 5001


def column_problems(ref_seq, edits, a=START, b=START + 1999):
    hap, cmap = haplotype.build(ref_seq, START, edits)
    r, h, k = haplotype.columns(cmap, a, b)
    problems = []
    placed = r[r > 0]
    if list(placed) != list(range(a, b + 1)):
        problems.append("reference positions aren't each start..end once, in order")
    if (h[k == haplotype.DELETED] != 0).any() or (r[k == haplotype.INSERTED] != 0).any():
        problems.append("a deletion must be a gap in the patient rows, an insertion one in the reference rows")
    for c in np.nonzero(k == haplotype.INVERTED)[0]:
        comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
        if hap[h[c] - 1] != comp[ref_seq[r[c] - START]]:
            problems.append(f"column {c}: inverted base isn't the complement")
            break
    for c in np.nonzero(k == haplotype.SAME)[0]:
        if hap[h[c] - 1] != ref_seq[r[c] - START]:
            problems.append(f"column {c}: unchanged base differs")
            break
    # The columns come from the coordinate map: each reference base's column holds
    # the haplotype base the map gives it, no haplotype base twice, inserted ones inserted.
    if len(set(h[h > 0].tolist())) != int((h > 0).sum()):
        problems.append("a haplotype position is in two columns")
    for c in range(len(r)):
        if r[c] and (h[c] or None) != cmap.to_hap(int(r[c])):
            problems.append(f"column {c}: reference {r[c]} -> haplotype {h[c]}, the map says {cmap.to_hap(int(r[c]))}")
            break
        if not r[c] and cmap.to_ref(int(h[c]))[0] is not None:
            problems.append(f"column {c}: haplotype {h[c]} isn't inserted")
            break
    return problems, (r, h, k)


problems, (r, h, k) = column_problems(REF, [])
check("columns: no edits is one column per base, both sides", not problems and (r == np.arange(START, START + 2000)).all()
      and (h == np.arange(1, 2001)).all())
problems, (r, h, k) = column_problems(REF, [Edit(START + 100, START + 99, "ACGTA", "insertion")])
ins = np.nonzero(k == haplotype.INSERTED)[0]
check("columns: an insertion opens 5 gap columns in the reference rows, right after its anchor base",
      not problems and len(ins) == 5 and r[ins[0] - 1] == START + 99 and r[ins[-1] + 1] == START + 100, str(problems))
problems, (r, h, k) = column_problems(REF, [Edit(START + 300, START + 349, "", "deletion")])
check("columns: a deletion opens 50 gap columns in the patient rows", not problems
      and (k == haplotype.DELETED).sum() == 50 and (h[r == START + 320] == 0).all(), str(problems))
problems, (r, h, k) = column_problems(REF, [Edit(START + 500, START + 599, "", "inversion")])
inv = np.nonzero(k == haplotype.INVERTED)[0]
check("columns: an inversion in reference order, the patient side reversed (and complementary)",
      not problems and len(inv) == 100 and (np.diff(r[inv]) == 1).all() and (np.diff(h[inv]) == -1).all(), str(problems))
bad = []
for trial in range(300):
    edits, pos = [], START + rng.randrange(0, 30)
    while pos < START + 1950:
        kind = rng.choice(["SNV", "insertion", "deletion", "inversion", "complex", "duplication"])
        n = rng.randrange(1, 40)
        end = min(pos + n - 1, START + 1999)
        e = {"SNV": Edit(pos, pos, "A" if REF[pos - START] != "A" else "G", kind),
             "insertion": Edit(pos, pos - 1, rand_seq(n), kind), "deletion": Edit(pos, end, "", kind),
             "inversion": Edit(pos, end, "", kind), "complex": Edit(pos, end, rand_seq(rng.randrange(0, 12)), kind),
             "duplication": Edit(pos, pos - 1, REF[max(0, pos - START - n):pos - START], kind)}[kind]
        edits.append(e)
        pos = (e.end if e.end >= e.start else e.start - 1) + 1 + rng.randrange(1, 200)
    a, b = START + rng.randrange(0, 100), START + 1999 - rng.randrange(0, 100)
    p, _ = column_problems(REF, edits, a, b)
    if p:
        bad.append((trial, p[0]))
check("columns: 300 random edit sets and stretches -- every base once, in order, gaps on the right side", not bad,
      str(bad[:2]))

# --- view data on a real analysis ---------------------------------------------------------------------
REFSEQ, MANE = str(helpers.THIRDPARTY / "refseq"), str(helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz")
real_db = os.path.join(helpers.REAL_HOME, "SpliceAI_reference_data", "spliceai-refdb", "hg38")
reader = refdb.Reader(real_db) if os.path.exists(os.path.join(real_db, "manifest.json")) else None
scores = refdb.ReferenceScores("hg38", str(helpers.HG38), reader=reader, cache_dir=os.path.join(tmp, "cache"),
                               refseq_dir=REFSEQ, mane_dir=str(helpers.MANE_DIR))
fasta = scores.predictor().annotator.ref_fasta
table = ca.GeneTable("hg38", REFSEQ, MANE)


def analysis(name, edits_for):
    gene = table.find(name)[0]
    record = fasta["chr" + gene.chrom]
    fetch = lambda a, b: record[a - 1:b].seq
    region = (gene.mask[0] - 5000, gene.mask[1] + 5000)
    ref = fetch(*region).upper()
    patient, _ = haplotype.build(ref, region[0], edits_for(gene, region, ref))
    compat = ca.check_compatibility(patient[3000:len(patient) - 3000], gene, fetch, len(record))
    for e in compat.alignment.events:
        e.apply = True
    return ca.run_analysis("hg38", gene, compat, scores, fetch, table)


def mrpl39_edits(gene, region, ref):
    exons = list(reversed(gene.canonical.exons))
    site = exons[2][0] - 1
    return [Edit(site, site, "A" if ref[site - region[0]] != "A" else "C", "SNV"),
            Edit(exons[6][1] + 700, exons[6][1] + 699, rand_seq(70), "insertion"),
            Edit(exons[8][1] + 300, exons[8][1] + 379, "", "deletion")]


result = analysis("MRPL39", mrpl39_edits)
d = ViewData.from_analysis(result)
plan, pred = result.plan, result.patient
ref_cols = (d.ref >= plan["start"]) & (d.ref <= plan["end"])
check("view data: reference scores on their columns are the database's (or live) values",
      np.array_equal(d.tracks["ref_acc"][ref_cols], result.wt_values[d.ref[ref_cols] - plan["start"], 0]))
placed = d.ref[d.ref > 0]
outside = (d.ref > 0) & ~ref_cols
check("view data: up to 3 kb of reference beyond the gene on each side (as far as the region goes), not scored",
      placed.min() <= plan["start"] - 1 and placed.max() >= plan["end"] + 1
      and plan["start"] - placed.min() <= 3000 and placed.max() - plan["end"] <= 3000
      and np.isnan(d.tracks["ref_acc"][outside]).all(), f"{placed.min()}..{placed.max()} vs {plan['start']}..{plan['end']}")
from spliceai_pipeline.view_data import FLANK_WHY, TranscriptRow  # noqa: E402

wide = TranscriptRow("wide annotated", "alternative", 0, given=np.ones(d.n, np.int64))
d.rows.append(wide)
check("...the flank is hatched as reference only, and never reported as transcribed positions without a score",
      any(why == FLANK_WHY for _, _, why in d.unscored()) and d.unscored_transcribed() == [])
d.rows.remove(wide)
hap_cols = (d.hap >= pred.start) & (d.hap <= pred.end)
check("view data: patient scores on their columns are the patient's prediction; delta = patient - reference",
      np.array_equal(d.tracks["pat_don"][hap_cols], pred.values[d.hap[hap_cols] - pred.start, 1])
      and np.allclose(d.tracks["delta_acc"], d.tracks["pat_acc"] - d.tracks["ref_acc"], equal_nan=True))
canonical = next(r for r in d.rows if r.group == "reference")
c = result.gene.canonical
check("view data: the canonical transcript's exons, in transcript order, on the right columns",
      [(b.start, b.end) for b in canonical.blocks] == list(reversed(c.exons))
      and all(d.ref[b.first] == min(b.start, b.end) for b in canonical.blocks))
check("view data: the patient's 70 bp insertion is 70 inserted columns, labelled INS1+n",
      (d.kind == haplotype.INSERTED).sum() == 70
      and d.label(int(np.nonzero(d.kind == haplotype.INSERTED)[0][4])) == "INS1+5")
check("completeness: every exonic and intronic column of the canonical and curated alternative transcripts "
      "has scores", d.unscored_transcribed() == [] and not np.isnan(d.tracks["ref_acc"][ref_cols]).any())
check("view data: the three differences are located on the columns",
      len(d.events) == 3 and all(0 <= a <= b < d.n for a, b, _ in d.events), str(d.events))
mid = canonical.blocks[3]
hover = d.hover(mid.first + 2)
check("details (spec 16 item 3): sections Position, Bases, Scores, Transcripts -- genomic and c. positions, both "
      "bases, both scores, delta, source, exon membership",
      [title for title, _ in d.details(mid.first + 2)] == ["Position", "Bases", "Scores", "Transcripts"]
      and "Genomic: chr21:" in hover and "c.: c." in hover and "Reference: " in hover and "Sample acceptor" in hover
      and "Delta donor" in hover
      and f"exon {mid.identity}" in hover, hover)
intron_col = canonical.blocks[3].last + 30
intron_hover = d.hover(intron_col)
merged = dict((what, rows) for rows, what in d.memberships(intron_col))
check("details: an intronic base gets its intronic c. position (e.g. c.123-30), and the transcripts that agree "
      "share one line ('WT predicted 1-5, Sample 1-5: intron 4')",
      any(line.strip().startswith("c.: c.") and ("+" in line or "-" in line[8:]) for line in intron_hover.splitlines())
      and any("WT predicted 1-5" in rows for rows in merged.values())
      and all(what.startswith("intron") for what in merged), intron_hover)

# Spec 16 items 1 and 2: the events described.
info = d.event_info
snv, ins, dele = (next(x for x in info if x["kind"] == kind) for kind in ("SNV", "insertion", "deletion"))
check("events (spec 16): one description per applied event -- type, size, coordinates, place in the canonical "
      "transcript with c., nearest donor and acceptor",
      len(info) == 3 and snv["label"].startswith("SNV 1 bp · chr21:") and snv["where"] == "intron 3"
      and "donor of exon 3: 1 bp" in snv["sites"] and snv["c"].startswith("c.") and "+1" in snv["c"]
      and dele["label"].startswith("DEL 80 bp") and dele["where"] == "intron 8" and ins["label"].startswith("INS 70 bp")
      and ins["where"] == "intron 6" and "acceptor of exon 9: 299 bp" in dele["sites"], " | ".join(x["label"] for x in info))
check("events: where inserted bases come from -- random bases are novel; a deletion or SNV puts none in; what "
      "sample 1 does there (exon 3 skipped at the SNV)", ins["origin"] == "novel (not found in the gene region)"
      and dele["origin"] == "" == snv["origin"] and "sample 1: exon 3 skipped" in snv["effect"]
      and all(x["line"].startswith(x["label"]) and x["effect"] in x["line"] for x in info),
      " | ".join(x["line"] for x in info))
csv_path = os.path.join(tmp, "view.csv")
d.write_csv(csv_path, 100, 4000)
back = ViewData.from_csv(csv_path)
check("CSV of a stretch, read back: the same bases, scores and exon membership (per-base view only, said so)",
      back.n == 3901 and np.array_equal(back.ref, d.ref[100:4001]) and "".join(back.pat_base) == "".join(d.pat_base[100:4001])
      and np.allclose(back.tracks["pat_acc"], d.tracks["pat_acc"][100:4001], atol=5e-7, equal_nan=True)
      and np.array_equal(back.rows[0].membership, d.rows[0].membership[100:4001]) and back.limited)
text = open(csv_path, encoding="utf-8").read()
old = os.path.join(tmp, "old.csv")
with open(old, "w", encoding="utf-8", newline="") as fh:
    fh.write(text.replace("sample", "patient"))
check("the CSV's header line carries the research-use note and says 'sample'; an export from before (saying "
      "'patient') still reads", ca.RESEARCH_NOTE in text.splitlines()[0] and "sample base" in text.splitlines()[1]
      and ViewData.from_csv(old).n == 3901)
try:
    ViewData.from_csv(os.path.join(helpers.REPO, "README.md"))
    check("a CSV that isn't the viewer's is refused", False)
except ValueError as exc:
    check("a CSV that isn't the viewer's is refused", "isn't a CSV exported" in str(exc))

# --- the window -----------------------------------------------------------------------------------
from PySide6.QtCore import QPoint, QPointF, QRectF, Qt  # noqa: E402
from PySide6.QtGui import QImage, QPainter  # noqa: E402
from PySide6.QtTest import QTest  # noqa: E402
from PySide6.QtWidgets import QApplication, QMessageBox  # noqa: E402

app = QApplication(sys.argv)
QMessageBox.warning = staticmethod(lambda *a, **k: QMessageBox.Ok)
from PySide6.QtGui import QFont, QFontMetrics  # noqa: E402

from spliceai_gui import config, theme  # noqa: E402
from spliceai_gui.stacked_view import DELTA, DNA, EVENTS, FOOTER, LEGEND, RULER, SCORES, TRANSCRIPT, exon_label, \
    orientation_text  # noqa: E402
from spliceai_gui.viewer_window import ViewerWindow  # noqa: E402


def wait_tiles(view, seconds=30):
    t0 = time.time()
    while time.time() - t0 < seconds:
        view.repaint()
        app.processEvents()
        if view.tiles_ready():
            return True
        time.sleep(0.02)
    return False


def schematic(view, width=1200):
    """Paints the exon schematic once (it records what it drew)."""
    image = QImage(width, max(view.content_height(), 100), QImage.Format_ARGB32)
    p = QPainter(image)
    view.paint(p, width)
    p.end()


t0 = time.perf_counter()
w = ViewerWindow(result, candidate=1)
w.resize(1300, 900)
w.show()
app.processEvents()
opened_ms = (time.perf_counter() - t0) * 1000
check("viewer (spec 15 item 8): exon schematic, stacked view, sequence comparison -- it opens on the schematic",
      [w.tabs.tabText(i) for i in range(w.tabs.count())] == ["Exon schematic", "Stacked view", "Sequence comparison"]
      and w.tabs.currentIndex() == 0)
check("viewer (item 2): the stacked view draws nothing and the sequence comparison isn't built until shown; open in "
      "under 1 s", not w.stacked._tiles and w.sequences is None and opened_ms < 1000, f"{opened_ms:.0f} ms")
check("summary card (spec 16 item 1): one line per applied event at the top of the viewer; the event list next to "
      "'Zoom to event'; Focus on (there are events)", w.summary.isVisibleTo(w) and w.summary_lines == [x["line"] for x in info]
      and w.event_combo.count() == 3 and w.zoom_event_button.isEnabled() and w.focus.isChecked()
      and not w.origin_button.isEnabled())
w.tabs.setCurrentIndex(1)
ready = wait_tiles(w.stacked)
first_event = d.events[0]
e0 = info[0]["event"]
k1, k2 = ca.flanking_exons(result, min(e0.start, e0.end), max(e0.start, e0.end))
flank = [b for b in canonical.blocks if k1 <= b.number <= k2]
a, b = w.stacked.visible_genomic()
centre = (a + b) / 2
check("the stacked view opens centred on the first event, its flanking exons in view (not the whole gene)",
      a <= min(bl.first for bl in flank) and max(bl.last for bl in flank) <= b and b - a < d.n / 2
      and abs(centre - (first_event[0] + first_event[1]) / 2) <= (b - a) * 0.06 + 2,
      f"{a}..{b}, event {first_event[:2]}, exons {k1}-{k2}")
kinds = [k for k, _, _ in w.stacked.rows()]
check("labels above the ruler: each event in view labelled with its description, the legend names the event "
      "types", kinds[:3] == [LEGEND, EVENTS, RULER] and w.stacked.label_boxes
      and w.stacked.label_boxes[0][1]["label"] == info[0]["label"])
focused = [(k, p.name if k == TRANSCRIPT else p) for k, _, p in w.stacked.rows() if k not in (LEGEND, EVENTS, RULER)]
check("Focus: only Reference DNA, Reference scores, the reference transcript, WT predicted 1, the sample's "
      "sequence, scores and delta, and Sample 1", focused == [
          (DNA, "ref"), (SCORES, "ref"), (TRANSCRIPT, canonical.name), (TRANSCRIPT, "WT predicted 1"), (DNA, "pat"),
          (SCORES, "pat"), (DELTA, None), (TRANSCRIPT, "Sample 1")], str(focused))
w.stacked.show_genomic(first_event[0] - 60, first_event[1] + 60)
top_delta, scale = w.stacked.delta_scale()
va, vb = w.stacked.visible()
image = QImage(w.stacked.width(), 60, QImage.Format_ARGB32)
painter = QPainter(image)
from spliceai_gui.exon_view import Axis, place_notes  # noqa: E402
from spliceai_gui.stacked_view import ROW_GAP, TrackPainter, exon_color, insertion_gaps, view_font  # noqa: E402

marks = TrackPainter(w.stacked.data, [], w.stacked.first, w.stacked.cols_per_px).delta(
    painter, QRectF(0, 0, w.stacked.width(), 40), va, vb, scale)
painter.end()
big = np.nanmax(np.abs(np.r_[d.tracks["delta_acc"][va:vb + 1], d.tracks["delta_don"][va:vb + 1]]))
check("delta row: scaled to the largest |delta| in view (that maximum is what the row says), positions with "
      "|delta| >= 0.1 marked", abs(top_delta - big) < 1e-6 and scale == max(big, 0.01) and (len(marks) > 0) == (big >= 0.1),
      f"max {top_delta:.3f}, {len(marks)} marks")
font = QFont(theme.mono_font())
font.setPointSize(8)
fm = QFontMetrics(font)
check("exon numbers on every exon block with room: 'Exon 14', else '14', else none",
      exon_label(fm, 14, 200) == "Exon 14" and exon_label(fm, 14, fm.horizontalAdvance("14") + 6) == "14"
      and exon_label(fm, 14, 6) == "" and exon_label(fm, 0, 200) == "novel")
check("the orientation is said above the stacked view, not on the ruler",
      "genomic orientation (− strand" in w.orientation_label.text() and "3' ◀── 5'" in w.orientation_label.text())
w.spliceai_alternatives.setChecked(True)
check("ticking an alternatives box turns Focus off, so what was asked for is shown", not w.focus.isChecked()
      and any(k == TRANSCRIPT and p.name == "Sample 2" for k, _, p in w.stacked.rows()))
w.spliceai_alternatives.setChecked(False)
w.focus.setChecked(True)
w.stacked.fit()
ready = ready and wait_tiles(w.stacked)
check("viewer: shown, the stacked view's tiles are drawn in the background and cached", ready
      and len(w.stacked._tiles) >= 1)
zoom_a = w.stacked.cols_per_px * 0.5
for cpp in (zoom_a, zoom_a * 0.7, zoom_a):          # leave a zoom before its tiles arrive, then come back
    w.stacked.cols_per_px = cpp
    w.stacked._clamp()
    w.stacked.repaint()
check("tiles asked for at a zoom left before they arrived are drawn when that zoom is shown again",
      wait_tiles(w.stacked, 10))
w.stacked.fit()
png, svg = os.path.join(tmp, "view.png"), os.path.join(tmp, "view.svg")
check("viewer: PNG and SVG from the same drawing", w.stacked.render_png(png) and os.path.getsize(png) > 5000
      and w.stacked.render_svg(svg) and "<svg" in open(svg, encoding="utf-8").read(2000))
check("PNG and SVG carry the research-use note in a footer under the drawing",
      QImage(png).height() == 2 * (w.stacked.total_height() + FOOTER)
      and ca.RESEARCH_NOTE.split(":")[0] in open(svg, encoding="utf-8").read())
a, b = w.stacked.visible()
check("viewer: the whole gene fits", a == 0 and b >= d.n - 2, f"{a}..{b} of {d.n}")
events = sorted(d.events)
seen = []
for _ in range(len(events) + 1):
    w.step_difference(1)
    a, b = w.stacked.visible_genomic()
    seen.append(next((i for i, e in enumerate(events) if a <= e[0] and e[1] <= b and b - a < d.n / 2), None))
w.step_difference(-1)
a, b = w.stacked.visible_genomic()
check("viewer: 'Next difference' goes through the differences in order, round again; 'Previous' goes back",
      seen == list(range(len(events))) + [0] and a <= events[-1][0] and events[-1][1] <= b, str(seen))
exon3 = next(blk for blk in canonical.blocks if blk.number == 3)
check("viewer: a c. position in the jump box goes there (c.1)",
      w.jump("c.1") and w.stacked.cursor_col == d.column_of_ref(c.cds[1]))
w.tabs.setCurrentIndex(0)
w.exons.exon_clicked.emit(exon3.first, exon3.last)
a, b = w.stacked.visible_genomic()
check("viewer: clicking an exon in the schematic zooms the stacked view to it", a <= exon3.first and exon3.last <= b
      and w.tabs.currentIndex() == 1)
w.stacked.zoom_to(exon3.first, exon3.first + 30)
w.stacked.render_png(os.path.join(tmp, "bases.png"))
check("viewer: at base level it still draws (letters and values)", 1 / w.stacked.cols_per_px >= 8)

# Item 3: short tooltips, the details in the side panel.
lines, worst = 0, ""
for x in range(5, 1290, 37):
    for y in range(5, w.stacked.total_height(), 9):
        text = w.stacked.tooltip(w.stacked.col_at(x), y)
        if text.count("\n") + 1 > lines:
            lines, worst = text.count("\n") + 1, text
top = 2
for kind, h, payload in w.stacked.rows():
    if kind == TRANSCRIPT and payload.group == "reference":
        break
    top += h + ROW_GAP
w.stacked.show_genomic(exon3.first, exon3.last)
block_tip = w.stacked.tooltip(w.stacked.to_display((exon3.first + exon3.last) // 2), top + 18)
check("tooltips (item 3): at most 4 short lines anywhere in the stacked view; on an exon: its number, length and "
      "status", lines <= 4 and block_tip.startswith(f"{canonical.name}, exon 3 (reference exon 3)")
      and "nt" in block_tip and block_tip.count("\n") <= 2, f"{lines} lines: {worst!r} | {block_tip!r}")
QTest.mouseClick(w.stacked, Qt.LeftButton, Qt.NoModifier, QPoint(300, top + 18))
app.processEvents()
panel = w.details_panel.toPlainText()
check("a click puts everything about the base or exon in the side panel, as a key/value table in sections: "
      "Exon, Position, Bases, Scores, Transcripts", all(s in panel for s in ("Exon", "Position", "Genomic", "chr21:",
                                                                             "Bases", "Scores", "Ref. acceptor",
                                                                             "Transcripts"))
      and "<table" in w.details_panel.toHtml() and len(panel.splitlines()) > 8, panel[:200])
check("the side panel is a splitter, about 320 px by default", abs(w.splitter.sizes()[1] - 320) <= 2,
      str(w.splitter.sizes()))
QTest.mouseMove(w.stacked, QPoint(500, top + 60))
w.details_panel.clear()
w.stacked.mouseMoveEvent(type("E", (), {"position": lambda self: QPointF(510, top + 60),
                                        "globalPosition": lambda self: QPointF(510, top + 60),
                                        "modifiers": lambda self: Qt.ShiftModifier})())
check("Shift + hover updates the side panel too", "chr21:" in w.details_panel.toPlainText())
w.stacked.show_genomic(exon3.first, exon3.last)
first = w.stacked.first
QTest.mousePress(w.stacked, Qt.LeftButton, Qt.NoModifier, QPoint(400, top + 60))
QTest.mouseMove(w.stacked, QPoint(300, top + 60))
QTest.mouseRelease(w.stacked, Qt.LeftButton, Qt.NoModifier, QPoint(300, top + 60))
panned = w.stacked.first > first
QTest.mouseDClick(w.stacked, Qt.LeftButton, Qt.NoModifier, QPoint(300, top + 60))
check("drag pans; double-click shows the whole gene again", panned and w.stacked.visible() == (0, w.stacked.data.n - 1),
      f"{first} -> {w.stacked.first}")
ref_row = next(r for r in w.stacked.data.rows if r.group == "reference")
cds_cols = sorted(w.stacked.data.columns_of(*result.gene.canonical.cds))
check("the reference transcript knows its coding part (UTRs drawn thinner)",
      ref_row.cds is not None and tuple(sorted(ref_row.cds)) == tuple(cds_cols)
      and min(b.first for b in ref_row.blocks) <= cds_cols[0] <= cds_cols[1] <= max(b.last for b in ref_row.blocks),
      f"{ref_row.cds} {cds_cols}")

# Item 1: the strand, and transcript orientation for the stacked view (MRPL39 is on the minus strand).
check("strand: the window and the orientation line say '− strand'; genomic orientation by default in the stacked view",
      "− strand" in w.windowTitle() and "genomic orientation (− strand" in orientation_text(w.stacked.data)
      and not w.stacked_orientation.isChecked() and w.stacked_orientation.isEnabled())
g_first = next(blk for blk in canonical.blocks if blk.number == 1)
g_second = next(blk for blk in canonical.blocks if blk.number == 2)
probe = (exon3.first + exon3.last) // 2
genomic_text = d.position_text(probe)
w.stacked_orientation.setChecked(True)
f = w.stacked.data
f_canon = next(r for r in f.rows if r.group == "reference")
f1 = next(blk for blk in f_canon.blocks if blk.number == 1)
f2 = next(blk for blk in f_canon.blocks if blk.number == 2)
comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N", "-": "-"}
event = d.events[0]
w.stacked.show_genomic(event[0], event[1])
va, vb = w.stacked.visible()
check("transcript orientation (minus strand): exon 1 now on the left, bases complemented, the line above says so",
      f.reversed and g_first.first > g_second.first and f1.first < f2.first
      and "transcript orientation 5'→3'" in orientation_text(f) and "transcript orientation" in w.orientation_label.text()
      and all(f.ref_base[w.stacked.to_display(k)] == comp[d.ref_base[k]] for k in range(probe, probe + 50)))
check("transcript orientation: tooltips give the same genomic and c. positions, exon numbers stay, the event is "
      "where it was", f.position_text(w.stacked.to_display(probe)) == genomic_text
      and f1.number == 1 and f1.identity == 1
      and va <= w.stacked.to_display(event[1]) and w.stacked.to_display(event[0]) <= vb
      and f.events[0][:2] == (d.n - 1 - event[1], d.n - 1 - event[0]))
w.stacked_orientation.setChecked(False)
check("...and back to genomic", not w.stacked.data.reversed and w.stacked.data is w.data)

# Item 1 and 9 on the exon schematic.
schematic(w.exons)
rows, canonical_exons = w.exons.entries()
axis = w.exons.axis(rows, 1200)
x1, x2 = axis.x(canonical_exons[0][0]), axis.x(canonical_exons[1][0])
check("exon schematic: transcript orientation by default -- exon 1 left of exon 2 for this minus-strand gene; the "
      "header says so", w.exon_orientation.isChecked() and x1 < x2 and "exon 1 on the left" in w.exons.header_text
      and "5' ──▶ 3'" in w.exons.header_text)
w.exon_orientation.setChecked(False)
schematic(w.exons)
axis_g = w.exons.axis(rows, 1200)
check("...genomic orientation: exon 1 on the right, the arrow points left",
      axis_g.x(canonical_exons[0][0]) > axis_g.x(canonical_exons[1][0]) and "3' ◀" in w.exons.header_text)
w.exon_orientation.setChecked(True)
schematic(w.exons)
kinds = [text for _, text, _ in w.exons.event_marks]
check("exon schematic (item 9): every applied event marked with its type and size",
      sorted(kinds) == sorted(["DEL 80 bp", "INS 70 bp", "SNV 1 bp"]), str(kinds))
x_event = w.exons.event_marks[0][0]
w.tabs.setCurrentIndex(0)
w.exons.click(QPointF(x_event, 18 + 20))
a, b = w.stacked.visible_genomic()
first_event = d.events[0]
check("clicking an event mark jumps to it in the stacked view", w.tabs.currentIndex() == 1
      and a <= first_event[0] and first_event[1] <= b)
notes = w.exons.notes
check("between the reference and WT, and WT and each sample candidate: what differs, or 'no change'",
      notes["WT predicted 1"] == ["no change"] and notes["Sample 1"] and notes["Sample 1"] != ["no change"]
      and any("skipped" in n or "cryptic" in n or "pseudoexon" in n for n in notes["Sample 1"]), str(notes)[:300])
check("no shading where the WT prediction matches the reference (MRPL39), no '!' on the events",
      w.exons.shaded == [] and not any(c_ for _, _, c_ in w.exons.event_marks))
tip = ""
for rect, kind, payload in w.exons._hit:
    if kind == "exon":
        tip = w.exons.tooltip_at(rect.center())
        break
check("exon schematic tooltips are short: exon number, length, status", tip.count("\n") <= 2 and "nt" in tip
      and ", exon " in tip, tip)

# A WT prediction that differs from the reference next to an event: shaded, and the event gets its "!".
import copy  # noqa: E402

fake = copy.copy(result)
wt1 = copy.deepcopy(next(t for t in result.transcripts if t.kind == "WT predicted" and t.rank == 1))
target = ca.locate(result, result.edits[0].start, result.edits[0].end)[1]
wt1.vs_ref = (wt1.vs_ref[0], [dict(r, status="end shifted +12") if r["ref"] == target else r for r in wt1.vs_ref[1]],
              *wt1.vs_ref[2:])
fake.transcripts = [wt1 if t.kind == "WT predicted" and t.rank == 1 else t for t in result.transcripts]
odd = ViewerWindow(fake)
schematic(odd.exons)
local = ca.local_concordance(fake)
check("the region where the WT prediction differs is shaded, and the event next to it gets '!' with the caution -- "
      "the same regions as the results' local concordance", odd.exons.shaded == [("exon", target)]
      and odd.exons.event_marks[0][2] and bool(local[0]["differing"])
      and "WT prediction differs from the reference in this region" in
      next(p[1] for r, k, p in odd.exons._hit if k == "event" and p[0] == 0), str(odd.exons.shaded))
region_tip = next(p for r, k, p in odd.exons._hit if k == "region")
c_ = result.gene.canonical
a_, b_ = (c_.exons if c_.strand == "+" else list(reversed(c_.exons)))[target - 1]
donor = float(result.wt_values[(b_ if c_.strand == "+" else a_) - result.plan["start"], 1])
check("hovering the shaded region gives SpliceAI's WT probability at the reference's own sites there",
      f"donor of exon {target} {donor:.2f}" in region_tip and "site threshold" in region_tip, region_tip)
odd.close()

w.tabs.setCurrentIndex(2)
app.processEvents()
html = w.sequences.text.toHtml() if w.sequences else ""
check("sequence view (built when first shown): cDNA alignment of reference, WT and the candidate; protein "
      "alignment with its change", "Reference" in html and "WT predicted 1" in html and "Sample 1" in html
      and "Protein" in html and "p.(" in html)
w.spliceai_alternatives.setChecked(True)
w.alternatives.setChecked(True)
check("viewer: 'SpliceAI alternatives' adds the WT and sample candidates 2 to 5",
      sum(1 for kind, _, p in w.stacked.rows() if kind == 4 and p.group in ("wt", "patient")) == 10)

# A plus-strand gene: the orientations agree. With a tandem duplication: where its bases come from (spec 16 item 2).
def usp25_edits(gene, region, ref):
    exons = list(gene.canonical.exons)
    x = exons[5][1] + 900
    return [Edit(exons[3][1] + 400, exons[3][1] + 479, "", "deletion"),
            Edit(x, x - 1, ref[x - 150 - region[0]:x - region[0]], "insertion")]


plus = analysis("USP25", usp25_edits)
pw = ViewerWindow(plus)
pd = pw.data
schematic(pw.exons)
p_rows, p_canon = pw.exons.entries()
p_axis = pw.exons.axis(p_rows, 1200)
p_ref = next(r for r in pd.rows if r.group == "reference")
pb1 = next(blk for blk in p_ref.blocks if blk.number == 1)
pb2 = next(blk for blk in p_ref.blocks if blk.number == 2)
pw.stacked.set_transcript_orientation(True)
check("plus strand: '+ strand', exon 1 on the left in both views; transcript orientation changes nothing",
      "+ strand" in pw.windowTitle() and not pw.stacked_orientation.isEnabled() and not pw.stacked.data.reversed
      and pb1.first < pb2.first and p_axis.x(p_canon[0][0]) < p_axis.x(p_canon[1][0])
      and "5' ──▶ 3'" in pw.exons.header_text)
dup = next((x for x in pd.event_info if x["kind"] == "duplication"), None)
src = dup["sources"][0] if dup and dup["sources"] else None
check("a tandem duplication's origin: 'tandem copy of chr:start-end', the 150 bases before it", dup is not None
      and dup["origin"].startswith("tandem copy of chr") and src is not None and src[1] - src[0] + 1 == 150
      and src[2] == 1 and dup["origin"] in dup["line"], dup["line"] if dup else str(pd.event_info))
pw.show()
pw.tabs.setCurrentIndex(1)
app.processEvents()
pw.event_combo.setCurrentIndex(pd.event_info.index(dup))
went = pw.go_to_origin()
sa, sb = pw.data.columns_of(src[0], src[1])
va, vb = pw.stacked.visible_genomic()
check("'Go to origin' zooms to the source region and outlines it ('Show origin' on)", went and va <= sa and sb <= vb
      and pw.stacked.show_origin and pw.origin_button.isChecked() and pw.origin_button.isEnabled())
d_first, d_last, _ = pd.events[pd.event_info.index(dup)]
pw.stacked.show_genomic(d_first - 20, d_last + 20)
tip = pw.stacked.tooltip(pw.stacked.to_display((d_first + d_last) // 2), 200)
check("hovering the inserted block shows its origin", "origin: tandem copy of" in tip, tip)


class Recorder:
    """A painter that notes what insertion_gaps() draws."""

    def __init__(self):
        self.texts, self.dashed = [], 0

    def fontMetrics(self):
        return QFontMetrics(view_font())

    def fillRect(self, *_):
        pass

    def setPen(self, *_):
        pass

    def drawText(self, _rect, _flags, text):
        self.texts.append(text)

    def drawLine(self, *_):
        self.dashed += 1


view = pw.stacked
first, last = next((a, b) for (a, b, _t), info in zip(view.data.events, pd.event_info) if info is dup)
view.show_genomic(first - 100, last + 100)
gap = Recorder()
insertion_gaps(gap, view.data, view.rows(), view.x_of, view.width())
check("spec 17 item 5: the gap of an insertion is labelled, and the transcripts run across it as a dashed line",
      gap.texts == ["150 bp inserted"] and gap.dashed >= 1, f"{gap.texts}, {gap.dashed} dashed")
twice = os.path.join(tmp, "twice.png")
view.render_png(twice, 800)
wide = QImage(twice).width()
pw.exons.render_png(twice, 800)
check("both exports are at twice the resolution, the schematic's research-use note under it",
      wide == 1600 and QImage(twice).width() == 1600
      and QImage(twice).height() == 2 * (pw.exons.content_height() + FOOTER), str(wide))
pw.exons.render_png(twice, 480)
labels = sorted((r.left(), r.right()) for r, k, _ in pw.exons._hit if k in ("event", "note") and r.height() == 16)
check("in a narrow schematic the two events' labels keep off each other and inside it",
      len(labels) == 2 and labels[0][1] <= labels[1][0] and labels[1][1] <= 480, str(labels))
fm = QFontMetrics(view_font())
notes = place_notes(fm, [(100, "intron 2 retained"), (105, "cryptic donor +139"), (880, "pseudoexon 117 nt")], 0, 900)
crowded = place_notes(fm, [(x, f"cryptic donor {x}") for x in range(0, 300, 10)], 0, 500)
check("notes over a row: none over another or past the edge; what doesn't fit is '+n more'",
      len(notes) == 3 and all(a[0] + a[1] <= b[0] for a, b in zip(notes, notes[1:]))
      and notes[-1][0] + notes[-1][1] <= 900 and crowded[-1][2].endswith("more")
      and crowded[-1][0] + crowded[-1][1] <= 500, str([n[:3] for n in crowded]))
narrow = [x for _, x in Axis([(1000 * i, 1000 * i + 100) for i in range(60)], 100, 300, False).points]
check("in a narrow schematic the introns give way: every exon keeps some width, all inside the view",
      narrow[-1] <= 300.5 and all(b > a for a, b in zip(narrow[::2], narrow[1::2])), f"{narrow[:4]} .. {narrow[-1]}")
saved_plus = bundle.save(plus, os.path.join(tmp, "plus"))
pw.close()
saved = bundle.save(result, os.path.join(tmp, "result"))
check("viewer: a results file dropped onto it opens", w.open_file(saved) and w.data.gene == "MRPL39")
check("viewer: the CSV export opens, per-base view only, with the limitation shown (no schematic)",
      w.open_file(csv_path) and w.data.limited and not w.tabs.isTabEnabled(0) and w.tabs.currentIndex() == 1)
w.resize(1800, 900)              # room for it with the offscreen platform's wide letters
app.processEvents()
w.splitter.setSizes([w.splitter.width() - 410, 410])
w.close()
again = ViewerWindow(saved)
again.resize(1300, 850)
again.show()
app.processEvents()
check("the side panel's width is kept in the settings (config.py DEFAULTS: viewer_details_width) and the next "
      "viewer opens with it", config.DEFAULTS["viewer_details_width"] == 320
      and config.load()["viewer_details_width"] == 410 and abs(again.splitter.sizes()[1] - 410) <= 2,
      f"{config.load()['viewer_details_width']}, {again.splitter.sizes()}")
again.close()

# Spec 16 item 3: the toolbar fits at 1366 px with 100, 125 and 150 % scaling -- every control shown whole, or in
# the toolbar's overflow menu, none squeezed. Laid out by the Windows platform (its fonts), never shown.
fit_script = r"""
import os, sys
os.environ["QT_QPA_PLATFORM"] = "windows" if sys.platform == "win32" else "offscreen"
from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QLineEdit, QToolButton
app = QApplication(sys.argv)
from spliceai_gui.viewer_window import ViewerWindow
w = ViewerWindow(sys.argv[1])
w.setAttribute(Qt.WA_DontShowOnScreen)
w.show()
for width in (1366, 1093, 911):
    w.resize(width, 760)
    for _ in range(4):
        app.processEvents()
    for bar in w.toolbars:
        ext = next((b for b in bar.findChildren(QToolButton) if b.objectName() == "qt_toolbar_ext_button"), None)
        overflow = ext is not None and ext.isVisible()
        for action in bar.actions():
            widget = bar.widgetForAction(action)
            if widget is None:
                continue
            if widget.isVisible():
                inside = widget.geometry().right() <= bar.width() and widget.geometry().left() >= 0
                least = widget.minimumWidth() if isinstance(widget, QLineEdit) else widget.sizeHint().width()
                if not inside or widget.width() < least:
                    print(f"CUT {width} {widget.metaObject().className()} {getattr(widget, 'text', lambda: '')()!r}"
                          f" {widget.geometry().getRect()} bar {bar.width()}")
            elif not overflow:
                print(f"LOST {width} {getattr(widget, 'text', lambda: '')()!r}")
        print(f"WIDTH {width} {bar.windowTitle()}: overflow {overflow}")
print("FIT-DONE")
"""
run = subprocess.run([sys.executable, "-W", "ignore", "-c", fit_script, saved_plus], env=dict(os.environ,
                     PYTHONPATH=str(helpers.REPO)), cwd=str(helpers.REPO), capture_output=True, text=True, timeout=300)
fit_lines = run.stdout.splitlines()
check("viewer toolbar (spec 16 item 3): at 1366 px with 100, 125 and 150 % scaling every control is shown whole or "
      "in the overflow menu, nothing squeezed; at 100 % nothing overflows",
      "FIT-DONE" in fit_lines and not any(x.startswith(("CUT", "LOST")) for x in fit_lines)
      and all(x.endswith("overflow False") for x in fit_lines if x.startswith("WIDTH 1366")),
      " | ".join(x for x in fit_lines if not x.startswith("FIT")) + run.stderr[-300:])



# Spec 17: a duplicated exon in the drawings, on both strands, and after reloading the results file.
def duplicate(k, flank):
    def edits(gene, region, ref):
        exons = list(gene.canonical.exons) if gene.strand == "+" else list(reversed(gene.canonical.exons))
        start, end = exons[k - 1][0] - flank, exons[k - 1][1] + flank
        return [Edit(end + 1, end, ref[start - region[0]:end + 1 - region[0]], "duplication")]

    return edits


for name, k, flank in (("VHL", 2, 40), ("TP53", 7, 87)):
    twin = analysis(name, duplicate(k, flank))
    tw = ViewerWindow(twin)
    schematic(tw.exons)
    sample_row = next(r for r in tw.data.rows if r.name == "Sample 1")
    block = next((b for b in sample_row.blocks if b.duplicate), None)
    reference = next(b for b in tw.data.rows[0].blocks if b.number == k)
    check(f"{name}: the stacked view has the copy in the colour of exon {k}, 'duplicate of exon {k}', not 'novel'",
          block is not None and block.identity == k and exon_color(block.identity) == exon_color(reference.identity)
          and exon_label(QFontMetrics(view_font()), k, 400, True) == f"duplicate of exon {k}"
          and exon_label(QFontMetrics(view_font()), k, QFontMetrics(view_font()).horizontalAdvance(f"dup {k}") + 4,
                         True) == f"dup {k}", str(block))
    png = os.path.join(tmp, f"{name}.png")
    tw.exons.render_png(png, 1200)
    boxes = [r for r, kind, tip in tw.exons._hit if kind == "note" and f"exon {k} (dup)" in tip and "Sample 1" in tip]
    drawn = QImage(png)
    corner = drawn.pixelColor(int(2 * boxes[0].left() + 3), int(2 * boxes[0].top() + 3)) if boxes else None
    row_y = boxes[0].top() if boxes else 0
    others = [r for r, kind, p_ in tw.exons._hit if kind == "exon" and r.top() == row_y]
    check(f"{name}: the schematic has a second box in that colour, labelled 'exon {k} (dup)', over no other exon",
          len(boxes) == 1 and corner == exon_color(k) and f"exon {k} (dup)" in tw.exons.notes["Sample 1"]
          and not any(r.intersects(boxes[0]) for r in others), f"{boxes}, {corner}, {tw.exons.notes['Sample 1']}")
    tw.close()
    back = ViewerWindow(bundle.save(twin, os.path.join(tmp, f"{name}_twin")))
    back_row = next(r for r in back.data.rows if r.name == "Sample 1")
    check(f"{name}: the same after reloading the results file", any(b.duplicate and b.identity == k for b in back_row.blocks))
    back.close()

# A results file opened with no SpliceAI, no FASTA, no database: a separate process.
script = r"""
import os, sys
os.environ["QT_QPA_PLATFORM"] = "offscreen"
import importlib.util
assert importlib.util.find_spec("spliceai") is None, "spliceai importable"
from PySide6.QtWidgets import QApplication
app = QApplication(sys.argv)
from spliceai_gui.viewer_window import ViewerWindow
w = ViewerWindow(sys.argv[1])
w.show(); app.processEvents()
seq = w.build_sequences()
ok = w.stacked.render_png(sys.argv[2]) and "Protein" in seq.text.toHtml()
print("VIEWER-OK" if ok else "VIEWER-FAILED")
"""
home = tempfile.mkdtemp()
env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH",)}
env.update(USERPROFILE=home, HOME=home, PYTHONPATH=str(helpers.REPO))
run = subprocess.run([sys.executable, "-W", "ignore", "-c", script, saved, os.path.join(tmp, "alone.png")],
                     env=env, cwd=home, capture_output=True, text=True, timeout=300)
check("a results file opens in the viewer with no SpliceAI, no FASTA and no database (another process, empty home)",
      "VIEWER-OK" in run.stdout, (run.stdout + run.stderr)[-300:])

# --- a 2.3 Mb gene ------------------------------------------------------------------------------------
if reader is not None:
    def dmd_edits(gene, region, ref):
        mid = (gene.mask[0] + gene.mask[1]) // 2
        return [Edit(mid, mid + 49, "", "deletion")]

    t0 = time.time()
    big = analysis("DMD", dmd_edits)
    built = time.time() - t0
    t0 = time.perf_counter()
    v = ViewerWindow(big)
    v.resize(1400, 900)
    v.show()
    app.processEvents()
    opened = (time.perf_counter() - t0) * 1000
    bd = v.data
    v.tabs.setCurrentIndex(1)
    v.stacked.fit()
    tiled = wait_tiles(v.stacked, 60)
    frames = []
    for i in range(30):
        v.stacked.first += (1 if i % 10 < 5 else -1) * v.stacked.width() * v.stacked.cols_per_px / 20
        v.stacked._clamp()
        t0 = time.perf_counter()
        v.stacked.repaint()
        frames.append((time.perf_counter() - t0) * 1000)
    frames.sort()
    image = os.path.join(tmp, "dmd.png")
    t0 = time.time()
    v.stacked.render_png(image)
    whole = time.time() - t0
    t0 = time.time()
    for x in range(0, 1400, 7):
        v.stacked.tooltip(v.stacked.col_at(x), 60)
    hover_ms = (time.time() - t0) / 200 * 1000
    check(f"a 2.3 Mb gene (DMD, {bd.n:,} columns, spec 15 item 2): open in under 1 s, its tiles drawn in the "
          "background, panning under 50 ms a frame, hover instant",
          opened < 1000 and tiled and frames[len(frames) // 2] < 50 and hover_ms < 20
          and bd.unscored_transcribed() == [],
          f"analysis {built:.0f} s, open {opened:.0f} ms, pan median {frames[len(frames) // 2]:.1f} ms "
          f"(max {frames[-1]:.1f}), whole-gene export {whole * 1000:.0f} ms, hover {hover_ms:.2f} ms")
    v.close()
else:
    helpers.skip("a 2.3 Mb gene", "the Phase 2 database isn't on this PC (live computation would take minutes)")

helpers.finish()
