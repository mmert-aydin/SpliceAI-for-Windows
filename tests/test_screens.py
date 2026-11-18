"""Spec 17 item 5: the whole program and the viewer at 1280x720, 1366x768,
1600x900 and 1920x1080, at 100%, 125%, 150% and 175% Windows scaling, and at
2560x1440 and 3840x2160 up to 200%: nothing cut off, text sharp; the Settings
pane no taller than its contents; every dialog and every "(?)" help popup
fitting the screen. A screenshot of each view at each combination goes to
docs/phase_reports/phase_13_screens/.

Each scaling runs in its own process on the real Windows platform, whose fonts
decide every width (the offscreen platform draws text as boxes), with the
windows laid out but never shown on screen. QT_SCALE_FACTOR, set before Qt
starts, multiplies the monitor's own ratio, so it is divided by it here: the
ratio each process runs at is the one it is named after, and is checked. One real run (TP53 with an insertion
and a deletion, from the database) is saved first and opened at each size.

    .venv\\Scripts\\python tests\\test_screens.py [--no-screens]
"""
import json
import os
import subprocess
import sys
import tempfile

import helpers
from helpers import check

SCREENS = ((1280, 720), (1366, 768), (1600, 900), (1920, 1080), (2560, 1440), (3840, 2160))
SCALES = (1.0, 1.25, 1.5, 1.75, 2.0)


def offered(width, scale):
    """Windows offers 200% only on the large screens."""
    return scale <= 1.75 or width >= 2560
OUT = helpers.REPO / "docs" / "phase_reports" / "phase_13_screens"
TASKBAR = 48            # px of the screen the Windows taskbar takes


def child(bundle_path, scale, save):
    """One scaling: every screen size; prints a JSON list of problems per size."""
    home = helpers.isolate_home()
    os.makedirs(home / ".spliceai_gui")
    with open(home / ".spliceai_gui" / "config.json", "w") as fh:
        json.dump({"fasta_path_hg19": str(helpers.HG19), "fasta_path_hg38": str(helpers.HG38), "build": "hg38",
                   "precomputed_mode": "none", "use_snpeff": False, "hide_reference_note": True,
                   "mane_dir": str(helpers.MANE_DIR), "cdna_show_intro": False}, fh)
    helpers.setup()
    import numpy as np
    from PySide6.QtCore import QPoint, Qt
    from PySide6.QtGui import QColor, QImage, QPainter
    from PySide6.QtWidgets import QAbstractButton, QAbstractScrollArea, QApplication, QLabel, QScrollArea, QToolBar

    from spliceai_gui import theme
    from spliceai_gui.main_window import MainWindow
    from spliceai_gui.stacked_view import ROW_GAP, TRANSCRIPT
    from spliceai_gui.viewer_window import ViewerWindow

    app = QApplication(sys.argv)
    theme.apply_theme(app)
    w = MainWindow()
    w.setAttribute(Qt.WA_DontShowOnScreen)
    w.show()
    tab = w.cdna_tab
    tab._on_continue_without_database()
    result = tab.open_results(bundle_path)
    tab.build_combo.setCurrentText(result.build)
    tab.gene_edit.setText(result.gene.name)
    tab.genes, tab.gene = [result.gene], result.gene
    tab._check_done(result.compat)
    tab.show_result(result)
    viewer = ViewerWindow(bundle_path)
    viewer.setAttribute(Qt.WA_DontShowOnScreen)
    viewer.show()

    def settle(n=4):
        for _ in range(n):
            app.processEvents()

    def cut_off(window, name):
        """Controls a window shows but doesn't give room to: outside it, or narrower than their text."""
        out = []
        frame = window.rect()
        for c in window.findChildren(QAbstractButton) + window.findChildren(QLabel):
            if not c.isVisible() or c.window() is not window or isinstance(c.parent(), QToolBar):
                continue
            text = c.text() or c.toolTip()[:30]
            if not text.strip():
                continue
            area = c.parent()
            while area is not None and area is not window and not isinstance(area, QScrollArea):
                area = area.parent()        # (a dialog's parent is in another window: stop at this one)
            area = None if area is window else area
            box = area.viewport() if area is not None else window
            right = c.mapTo(box, QPoint(c.width(), 0)).x()
            if right > box.width() + 1:
                out.append(f"{name}: '{text[:30]}' ends at {right}, past {box.width()}")
            elif isinstance(c, QLabel):
                if not c.wordWrap() and c.sizePolicy().horizontalPolicy().name != "Ignored" \
                        and c.width() < c.minimumSizeHint().width():
                    out.append(f"{name}: label '{text[:30]}' is {c.width()} px, needs {c.minimumSizeHint().width()}")
            elif c.width() < c.minimumSizeHint().width():
                out.append(f"{name}: button '{text[:30]}' is {c.width()} px, needs {c.minimumSizeHint().width()}")
            elif area is None and not frame.contains(c.mapTo(window, QPoint(0, 0))):
                out.append(f"{name}: '{text[:30]}' outside the window")
        return out

    report = {}
    scroll = w.centralWidget()
    if abs(w.devicePixelRatioF() - scale) > 0.01:
        report["ratio"] = [f"running at {w.devicePixelRatioF()}, not {scale}"]
    popups = make_popups(app, w, tab, viewer, settle)
    for sw, sh in SCREENS:
        if not offered(sw, scale):
            continue
        lw, lh = round(sw / scale), round((sh - TASKBAR) / scale)
        problems = []
        name = f"{sw}x{sh}_{round(scale * 100)}"
        for window, title in ((w, "main window"), (viewer, "viewer")):
            hint = window.minimumSizeHint()
            if hint.width() > lw or hint.height() > lh:
                problems.append(f"{title}: needs at least {hint.width()}x{hint.height()}, the screen has {lw}x{lh}")
            window.resize(lw, lh)
        settle()
        for index, label in ((0, "vcf"), (w.tabs.indexOf(tab), "cdna")):
            w.tabs.setCurrentIndex(index)
            settle()
            if scroll.horizontalScrollBar().maximum():
                problems.append(f"{label}: scrolls {scroll.horizontalScrollBar().maximum()} px sideways")
            problems += cut_off(w, label)
            if label == "vcf":
                need = w.panes.widget(2).sizeHint().height()
                if w.panes.sizes()[2] > need + 2:
                    problems.append(f"vcf: the Settings pane is {w.panes.sizes()[2]} px, its contents {need}")
            if save:
                w.grab().save(str(OUT / f"{name}_{label}.png"))
            # What is folded away has to fit as well.
            for button in (w.advanced_button, tab.advanced):
                button.setChecked(True)
            settle()
            if scroll.horizontalScrollBar().maximum():
                problems.append(f"{label}, Advanced open: scrolls {scroll.horizontalScrollBar().maximum()} px sideways")
            problems += cut_off(w, label + ", Advanced open")
            for button in (w.advanced_button, tab.advanced):
                button.setChecked(False)
            settle()
            if label == "cdna":
                # The compatibility area, further down the tab.
                scroll.verticalScrollBar().setValue(tab.check_button.mapTo(scroll.widget(), QPoint(0, 0)).y() - 8)
                settle()
                if save:
                    w.grab().save(str(OUT / f"{name}_compatibility.png"))
                scroll.verticalScrollBar().setValue(0)
        for title, d, opens in popups:
            hint = d.minimumSizeHint()
            if hint.width() > lw or hint.height() > lh:
                problems.append(f"{title}: needs at least {hint.width()}x{hint.height()}, the screen has {lw}x{lh}")
                continue
            # Opened on this screen it is fitted to it (theme.fit_to_screen does that for every window).
            d.resize(theme.fitted_size(opens, lw, lh))
            settle(2)
            if d.width() > lw or d.height() > lh:
                problems.append(f"{title}: opens at {d.width()}x{d.height()}, the screen has {lw}x{lh}")
            elif d.size() == opens and (title == "welcome" or title.startswith("help")):
                # A text that fits the screen is shown whole: no scroll bar.
                bars = [a for a in d.findChildren(QAbstractScrollArea) if a.verticalScrollBar().isVisible()]
                if bars:
                    problems.append(f"{title}: scrolls although the screen has room for all of it")
            problems += cut_off(d, title)
        for index, label in ((0, "schematic"), (1, "stacked"), (2, "sequences")):
            viewer.tabs.setCurrentIndex(index)
            settle()
            if index == 1:
                view = viewer.stacked
                for _ in range(400):
                    view.repaint()
                    app.processEvents()
                    if view.tiles_ready():
                        break
                ratio = view.devicePixelRatioF()
                tile = next(iter(view._tiles.values()), None) if view._tiles else None
                if tile is None or tile.devicePixelRatio() != ratio or tile.width() != round(256 * ratio):
                    problems.append(f"stacked: a tile isn't at the screen's ratio {ratio}: "
                                    f"{tile and (tile.width(), tile.devicePixelRatio())}")
                shot = view.grab().toImage().convertToFormat(QImage.Format_RGB32)
                direct = QImage(shot.width(), shot.height(), QImage.Format_RGB32)
                direct.setDevicePixelRatio(ratio)
                direct.fill(QColor(theme.SURFACE))
                p = QPainter(direct)
                view.paint(p, view.width())
                p.end()
                # The rows with text in them, from the screen, against the same rows drawn straight at this ratio.
                y = 2
                for kind, h, payload in view.rows():
                    if kind == TRANSCRIPT:
                        top, bottom = round((y + 12) * ratio), round((y + 26) * ratio)
                        a = np.frombuffer(shot.constBits(), np.uint8).reshape(shot.height(), -1, 4)[top:bottom]
                        b = np.frombuffer(direct.constBits(), np.uint8).reshape(direct.height(), -1, 4)[top:bottom]
                        difference = float(np.abs(a.astype(int) - b.astype(int)).mean())
                        if difference > 6:
                            problems.append(f"stacked: {payload.name} on screen differs from a direct drawing by "
                                            f"{difference:.1f} of 255 (blurred?)")
                    y += h + ROW_GAP
            else:
                drawn = [(k, r) for r, k, _ in viewer.exons._hit if k in ("note", "event")]
                past = [k for k, r in drawn if r.right() > viewer.exons.width() + 1]
                if past:
                    problems.append(f"schematic: {len(past)} painted labels run past its right edge")
                labels = sorted((r.left(), r.right()) for r, k, _ in viewer.exons._hit
                                if k in ("event", "note") and r.top() > 30 and r.height() == 16)
                if any(b[0] < a[1] for a, b in zip(labels, labels[1:])):
                    problems.append(f"schematic: event labels drawn over each other: {labels}")
            problems += cut_off(viewer, label)
            if save:
                viewer.grab().save(str(OUT / f"{name}_{label}.png"))
        report[name] = problems
    for _, d, _ in popups:
        d.close()
    print("REPORT " + json.dumps(report), flush=True)
    viewer.close()
    from PySide6.QtWidgets import QMessageBox
    QMessageBox.question = staticmethod(lambda *a, **k: QMessageBox.Yes)
    w.close()


def make_popups(app, w, tab, viewer, settle):
    """[(title, dialog, fixed)]: every dialog the program opens, and one message box per "(?)" help (its text is
    the button's tooltip, as a click shows it). fixed: a message box, which can't be made smaller than it is."""
    from types import SimpleNamespace

    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QToolButton

    from spliceai_gui.help_button import help_dialog

    from spliceai_gui.cdna_tab import IntroDialog
    from spliceai_gui.download_dialog import ReferenceDownloadDialog
    from spliceai_gui.first_launch_dialog import FirstLaunchDialog
    from spliceai_gui.gene_import_dialog import GeneImportDialog
    from spliceai_gui.license_dialog import LicenseDialog
    from spliceai_gui.mane_download_dialog import ManeDownloadDialog
    from spliceai_gui.refdb_wizard import RefdbDialog, RefdbWizard
    from spliceai_gui.snpeff_download_dialog import SnpEffDownloadDialog
    from spliceai_gui.spliceai_setup_dialog import SpliceAISetupDialog
    from spliceai_gui.transcript_effect_dialog import TranscriptEffectDialog
    from spliceai_pipeline import vcf_bridge

    gene = tab.gene
    row = SimpleNamespace(chrom="chr17", pos=7674221, ref="G", alt="A", variant_id="chr17:7674221:G:A", gene=gene.name,
                          source="", ds_ag=0.0, ds_al=0.0, ds_dg=0.0, ds_dl=0.0)
    vcf = "\n".join(["##fileformat=VCFv4.2", "\t".join("#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT S".split()),
                     "\t".join("chr17 7674221 . G A . . . GT 0/1".split())]) + "\n"
    makers = [
        ("welcome", lambda: FirstLaunchDialog(w)), ("licence", lambda: LicenseDialog(w)),
        ("cDNA intro", lambda: IntroDialog(tab)), ("reference download", lambda: ReferenceDownloadDialog(w)),
        ("MANE download", lambda: ManeDownloadDialog(w)), ("SnpEff download", lambda: SnpEffDownloadDialog(w)),
        ("SpliceAI setup", lambda: SpliceAISetupDialog(w)), ("gene list", lambda: GeneImportDialog(w)),
        ("database", lambda: RefdbDialog(w.refdb, w)), ("database wizard", lambda: RefdbWizard(w.refdb, w)),
        ("transcript effect", lambda: TranscriptEffectDialog(
            w, "hg38", [row], lambda r: ([gene], []), lambda r, g: [],
            lambda r: vcf_bridge.genotype(vcf, r.chrom, r.pos, r.ref, r.alt), lambda g, n: (60, "estimate"), "masked")),
    ]
    out, seen = [], set()           # [(title, dialog, the size it opens at)]
    for title, make in makers:
        try:
            d = make()
        except Exception as exc:        # reported, so a dialog that can't be built here is seen, not skipped
            print(f"{title}: could not be built here: {exc!r}", file=sys.stderr)
            continue
        out.append((title, d))
    for window in (w, viewer):
        for b in window.findChildren(QToolButton):
            tip = b.toolTip()
            if b.text().endswith("(?)") and tip and tip not in seen:
                seen.add(tip)
                out.append((f"help '{b.text()[:30]}'", help_dialog(window, b.text(), tip)))
    for _, d in out:
        d.setAttribute(Qt.WA_DontShowOnScreen)
        d.show()
    settle()
    return [(title, d, d.size()) for title, d in out]


def make_bundle(path):
    helpers.setup()
    from pyfaidx import Fasta

    from spliceai_pipeline import bundle, cdna_analysis as ca, refdb
    from spliceai_pipeline.mane import find_cached_summary
    from spliceai_pipeline.wt_concordance import in_order

    build = "hg38"
    table = ca.GeneTable(build, str(helpers.THIRDPARTY / "refseq"), find_cached_summary(str(helpers.MANE_DIR)))
    gene = table.find("TP53")[0]
    record = Fasta(str(helpers.HG38))["chr" + gene.chrom]
    fetch = lambda a, b: record[a - 1:b].seq
    exons = in_order(gene.canonical)
    at = exons[3][0] - 400          # intron 4 (minus strand: below exon 4), a 193 bp insertion
    inserted = "".join("ACGT"[(i * i + 3 * i) % 4] for i in range(193))
    text = f"insertion chr17:{at} {inserted}\ndeletion chr17:{exons[6][0] - 300}-{exons[6][0] - 180}"
    compat = ca.from_coordinates(text, gene, build, fetch, len(record))
    root = refdb.default_root()
    scores = refdb.ReferenceScores(build, str(helpers.HG38), reader=refdb.Reader(os.path.join(root, build)),
                                   refseq_dir=str(helpers.THIRDPARTY / "refseq"), mane_dir=str(helpers.MANE_DIR))
    result = ca.run_analysis(build, gene, compat, scores, fetch, table)
    return bundle.save(result, path)


def main():
    save = "--no-screens" not in sys.argv
    helpers.watchdog(1500)
    os.makedirs(OUT, exist_ok=True)
    path = make_bundle(os.path.join(tempfile.mkdtemp(prefix="spliceai_screens_"), "tp53"))
    problems = {}
    platform = dict(os.environ, QT_QPA_PLATFORM="windows" if sys.platform == "win32" else "offscreen")
    native = float(subprocess.run(
        [sys.executable, "-c", "from PySide6.QtGui import QGuiApplication as A; a = A([]); "
                               "print(a.primaryScreen().devicePixelRatio())"],
        env=platform, capture_output=True, text=True).stdout.strip())
    print(f"(Qt's own ratio on this monitor: {native:.0%})")
    for scale in SCALES:
        env = dict(platform, QT_SCALE_FACTOR=str(scale / native), QT_SCALE_FACTOR_ROUNDING_POLICY="PassThrough")
        out = subprocess.run([sys.executable, "-W", "ignore", __file__, "--child", path, str(scale)]
                             + ([] if save else ["--no-screens"]), env=env, capture_output=True, text=True,
                             encoding="utf-8", errors="replace")
        line = next((x for x in out.stdout.splitlines() if x.startswith("REPORT ")), None)
        if line is None:
            print(out.stdout[-1500:], out.stderr[-1500:])
            problems[f"scale {scale}"] = ["the child process failed"]
        else:
            problems.update(json.loads(line[7:]))
    for name, found in problems.items():
        check(f"{name}: nothing cut off, the text sharp", not found, "; ".join(found[:6]) + (
            f" ... and {len(found) - 6} more" if len(found) > 6 else ""))
    helpers.finish()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--child":
        child(sys.argv[2], float(sys.argv[3]), "--no-screens" not in sys.argv)
    else:
        main()
