"""Spec 15 item 2: the viewer's speed, measured (a report, not a check; the
checks are in test_viewer.py). For the largest gene (DMD, about 2.2 Mb) and a
typical one (MRPL39, 22 kb): the time from "Open viewer" to the first view
painted, and frame times while panning and zooming the stacked view at
several zoom levels. Laid out by the Windows platform, never shown on screen.

    .venv\\Scripts\\python tests\\perf_viewer.py [label]
"""
import os
import random
import statistics
import sys
import time

import helpers

home = helpers.isolate_home()
os.environ["QT_QPA_PLATFORM"] = "windows" if sys.platform == "win32" else "offscreen"
helpers.setup()

from PySide6.QtCore import Qt  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from spliceai_pipeline import cdna_analysis as ca, haplotype, refdb  # noqa: E402
from spliceai_pipeline.haplotype import Edit  # noqa: E402

app = QApplication(sys.argv)
from spliceai_gui import exon_view, stacked_view  # noqa: E402
from spliceai_gui.viewer_window import ViewerWindow  # noqa: E402

REFSEQ, MANE = str(helpers.REPO / "refseq"), str(helpers.MANE_DIR / "MANE.GRCh38.v1.5.summary.txt.gz")
db = os.path.join(helpers.REAL_HOME, "SpliceAI_reference_data", "spliceai-refdb", "hg38")
scores = refdb.ReferenceScores("hg38", str(helpers.HG38), reader=refdb.Reader(db), refseq_dir=REFSEQ,
                               mane_dir=str(helpers.MANE_DIR))
fasta = scores.predictor().annotator.ref_fasta
table = ca.GeneTable("hg38", REFSEQ, MANE)
rng = random.Random(5)

painted = []
for cls in (stacked_view.StackedView, exon_view.ExonView):
    original = cls.paintEvent

    def hooked(self, event, original=original):
        original(self, event)
        painted.append(time.perf_counter())
    cls.paintEvent = hooked


def analysis(name):
    gene = table.find(name)[0]
    record = fasta["chr" + gene.chrom]
    fetch = lambda a, b: record[a - 1:b].seq
    region = (gene.mask[0] - 5000, gene.mask[1] + 5000)
    ref = fetch(*region).upper()
    mid = (gene.mask[0] + gene.mask[1]) // 2
    insert = "".join(rng.choice("ACGT") for _ in range(70))
    sample, _ = haplotype.build(ref, region[0], [Edit(mid, mid + 49, "", "deletion"),
                                                 Edit(mid + 3000, mid + 2999, insert, "insertion")])
    compat = ca.check_compatibility(sample, gene, fetch, len(record))
    for e in compat.alignment.events:
        e.apply = True
    return ca.run_analysis("hg38", gene, compat, scores, fetch, table)


def open_viewer(result):
    """(ms from 'Open viewer' to the first finished paint, the window)."""
    painted.clear()
    t0 = time.perf_counter()
    w = ViewerWindow(result)
    w.setAttribute(Qt.WA_DontShowOnScreen)
    w.resize(1400, 900)
    w.show()
    while not painted and time.perf_counter() - t0 < 120:
        app.processEvents()
    return (painted[0] - t0) * 1000, w


def settle_ms(view, limit=10):
    """After a zoom step: ms until every tile in view is drawn (the background thread's part)."""
    if not hasattr(view, "tiles_ready"):
        return float("nan")
    t0 = time.perf_counter()
    while not view.tiles_ready() and time.perf_counter() - t0 < limit:
        view.repaint()
        app.processEvents()
        time.sleep(0.002)
    return (time.perf_counter() - t0) * 1000


def frames(view, shown, n=30):
    """Frame times (ms) of n pans by a twentieth of the width, then n zoom steps, with `shown` columns in view."""
    view.cols_per_px = max(shown / view.width(), 1 / 40)
    view.first = max(0.0, view.data.n / 2 - shown / 2)
    view._clamp()
    pan, zoom, settle = [], [], []
    for i in range(n):
        view.first += (1 if i % 10 < 5 else -1) * view.width() * view.cols_per_px / 20
        view._clamp()
        t0 = time.perf_counter()
        view.repaint()
        pan.append((time.perf_counter() - t0) * 1000)
    for i in range(n):
        view.cols_per_px *= 0.8 if i % 2 else 1.25
        view._clamp()
        t0 = time.perf_counter()
        view.repaint()
        zoom.append((time.perf_counter() - t0) * 1000)
    for i in range(6):                  # zoom levels not shown before: until their tiles are all drawn
        view.cols_per_px *= 0.93
        view._clamp()
        settle.append(settle_ms(view))
    return pan, zoom, settle


label = sys.argv[1] if len(sys.argv) > 1 else "now"
lines = [f"Viewer speed ({label}; Qt platform {app.platformName()}, window 1400 x 900)", "",
         "| Gene | Columns | Open viewer -> first view | Columns in view | Pan: median / max | Zoom: median / max | "
         "New zoom's tiles all drawn: median |",
         "|---|---|---|---|---|---|---|"]
for name in ("DMD", "MRPL39"):
    result = analysis(name)
    ms, w = open_viewer(result)
    view = w.stacked
    w.tabs.setCurrentIndex(next(i for i in range(w.tabs.count()) if w.tabs.widget(i).isAncestorOf(view)))
    for _ in range(5):
        app.processEvents()
    n, first = view.data.n, True
    for shown in (n, 200_000, 10_000, 1_000, 150):
        if shown > n:
            continue
        view.repaint()                  # warm (tiles or pyramids) at this zoom
        pan, zoom, settle = frames(view, shown)
        settled = f"{statistics.median(settle):.0f} ms" if settle and settle[0] == settle[0] else "--"
        lines.append(f"| {name if first else ''} | {f'{n:,}' if first else ''} | {f'{ms:,.0f} ms' if first else ''} | "
                     f"{shown:,} | {statistics.median(pan):.1f} / {max(pan):.1f} ms | "
                     f"{statistics.median(zoom):.1f} / {max(zoom):.1f} ms | {settled} |")
        first = False
    w.close()
print("\n".join(lines))
if len(sys.argv) > 2:
    with open(sys.argv[2], "a", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n\n")
helpers.finish()
